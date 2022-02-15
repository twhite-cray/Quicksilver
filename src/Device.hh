#pragma once
#include "MaterialDatabase.hh"
#include "MC_Cell_State.hh"
#include "MC_Facet_Adjacency.hh"
#include "MC_Particle.hh"
#include "MC_RNG_State.hh"
#include "MC_Tally_Event.hh"
#include "Messages.hh"
#include "NuclearData.hh"

class MC_Base_Particle;
class MonteCarlo;

struct DeviceFacet {
  double4 plane;
  int3 point;
  int adjacentCell;
  int adjacentDomain;
  int adjacentFacet;
  int currentDomain;
  int neighbor;
  MC_Subfacet_Adjacency_Event::Enum event;
};

struct DeviceCell {
  static constexpr int numFacets = 24;
  static constexpr int numQuadPoints = 14;
  DeviceCell &operator=(const MC_Cell_State &that)
  {
    cellNumberDensity = that._cellNumberDensity;
    material = that._material;
    return *this;
  }
  double *__restrict__ totals;
  double *__restrict__ groupTallies;
  DeviceFacet facets[numFacets];
  double cellNumberDensity;
  int quadPoints[numQuadPoints];
  int material;
};

struct DeviceDomain {
  DeviceCell *__restrict__ cells;
  double3 *__restrict__ nodes;
  int *__restrict__ neighbors;
};

struct DeviceIsotope {
  DeviceIsotope &operator=(const Isotope &that)
  {
    atomFraction = that._atomFraction;
    gid = that._gid;
    return *this;
  }
  double atomFraction;
  int gid;
};

struct DeviceMaterial {
  DeviceIsotope *__restrict__ isos;
  double mass;
  int isoSize;
};

struct DeviceReaction {
  double *__restrict__ crossSections;
  NuclearDataReaction::Enum type;
};

struct DeviceNuclearDataIsotope {
  DeviceReaction *__restrict__ reactions;
};

struct DeviceParticle {
  DeviceParticle &operator=(const MC_Base_Particle &that);
  DeviceParticle &operator=(const MC_Particle &that);

  __device__  DeviceParticle &operator=(const MessageParticle &that)
  {
    coordinate = that.coordinate;
    velocity = that.velocity;
    kineticEnergy = that.kineticEnergy;
    weight = that.weight;
    timeToCensus = that.timeToCensus;
    age = that.age;
    numMeanFreePaths = that.numMeanFreePaths;
    numSegments = that.numSegments;
    randomNumberSeed = that.randomNumberSeed;
    identifier = that.identifier;
    lastEvent = MC_Tally_Event::Facet_Crossing_Communication;
    numCollisions = that.numCollisions;
    breed = that.breed;
    domain = that.domain;
    cell = that.cell;
    return *this;
  }

  bool operator==(const MC_Base_Particle &that) const;
  void set(MC_Base_Particle &that) const;

  void set(MC_Particle &that) const
  {
    that.coordinate.x = coordinate.x;
    that.coordinate.y = coordinate.y;
    that.coordinate.z = coordinate.z;
    that.velocity.x = velocity.x;
    that.velocity.y = velocity.y;
    that.velocity.z = velocity.z;
    {
      const double divSpeed = 1.0/sqrt(velocity.x*velocity.x+velocity.y*velocity.y+velocity.z*velocity.z);
      that.direction_cosine.alpha = divSpeed*velocity.x;
      that.direction_cosine.beta = divSpeed*velocity.y;
      that.direction_cosine.gamma = divSpeed*velocity.z;
    }
    that.kinetic_energy = kineticEnergy;
    that.weight = weight;
    that.time_to_census = timeToCensus;
    that.totalCrossSection = 0;
    that.age = age;
    that.num_mean_free_paths = numMeanFreePaths;
    that.mean_free_path = 0;
    that.segment_path_length = 0;
    that.random_number_seed = randomNumberSeed;
    that.identifier = identifier;
    that.last_event = lastEvent;
    that.num_collisions = numCollisions;
    that.num_segments = numSegments;
    that.breed = breed;
    that.energy_group = 0;
    that.domain = domain;
    that.cell = cell;
    that.facet = 0;
    that.normal_dot = 0;
  }

  double3 coordinate;
  double3 velocity;
  double kineticEnergy;
  double weight;
  double timeToCensus;
  double age;
  double numMeanFreePaths;
  double numSegments;
  uint64_t randomNumberSeed;
  uint64_t identifier;
  MC_Tally_Event::Enum lastEvent;
  int numCollisions;
  int breed;
  int domain;
  int cell;
};


struct Device {

  Device():
    domains(nullptr),
    mats(nullptr),
    isotopes(nullptr),
    tallies(nullptr),
    particleSizes(nullptr),
    processing(nullptr),
    processed(nullptr),
    extras(nullptr),
    nuBar(0),
    logLow(0),
    divDelta(0),
    timeStep(0),
    numGroups(0),
    reactionSize(0)
  {}

  Device(const Device &) = default;
  void init(MonteCarlo &mc);
  void cycleInit(MonteCarlo &mc);
  void cycleFinalize(MonteCarlo &mc);

  void collide(const NuclearDataReaction::Enum type, const double energyIn, const double mass, double *__restrict__ const energyOut, double *__restrict__ const angleOut, int &nOut, uint64_t &seed) const
  {
    switch(type) {
      case NuclearDataReaction::Scatter:
        {
          nOut = 1;
          energyOut[0] = energyIn*(1.0-(rngSample(&seed)*(1.0/mass)));
          angleOut[0] = rngSample(&seed)*2.0-1.0;
        }
        break;
      case NuclearDataReaction::Absorption:
        break;
      case NuclearDataReaction::Fission:
        {
          nOut = int(nuBar+rngSample(&seed));
          for (int i = 0; i < nOut; i++) {
            const double ran = rngSample(&seed)/2.0+0.5;
            energyOut[i] = 20.0*ran*ran;
            angleOut[i] = rngSample(&seed)*2.0-1.0;
          }
        }
        break;
      case NuclearDataReaction::Undefined:
        abort();
    }
  }

  int getEnergyGroup(const double e)
  {
    const int i = int((log(e)-logLow)*divDelta);
    return std::max(0,std::min(numGroups,i));
  }

  DeviceDomain *__restrict__ domains;
  DeviceMaterial *__restrict__ mats;
  DeviceNuclearDataIsotope *__restrict__ isotopes;

  enum Tallies { ABSORB, CENSUS, COLLISION, ESCAPE, FISSION, PRODUCE, SEGMENTS, SCATTER, TALLIES_SIZE };
  long *__restrict__ tallies;

  enum ParticleSizes { PROCESSING, PROCESSED, EXTRAS, PARTICLE_SIZES_SIZE };
  int *__restrict__ particleSizes;
  DeviceParticle *__restrict__ processing;
  DeviceParticle *__restrict__ processed;
  DeviceParticle *__restrict__ extras;

  double nuBar;
  double logLow;
  double divDelta;
  double timeStep;
  int numGroups;
  int reactionSize;
};


