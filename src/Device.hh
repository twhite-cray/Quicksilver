#pragma once
#include "MaterialDatabase.hh"
#include "MC_Cell_State.hh"
#include "MC_RNG_State.hh"
#include "NuclearData.hh"

class MC_Base_Particle;
class MC_Particle;
class MonteCarlo;

struct DeviceCell {
  static constexpr int numFacets = 24;
  static constexpr int numQuadPoints = 14;
  DeviceCell &operator=(const MC_Cell_State &that)
  {
    cellNumberDensity = that._cellNumberDensity;
    material = that._material;
    return *this;
  }
  double *totals;
  double *groupTallies;
  double4 facets[numFacets];
  double cellNumberDensity;
  int3 facetPoints[numFacets];
  int quadPoints[numQuadPoints];
  int material;
};

struct DeviceDomain {
  DeviceCell *cells;
  double3 *nodes;
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
  DeviceIsotope *isos;
  double mass;
  int isoSize;
};

struct DeviceReaction {
  double *crossSections;
  NuclearDataReaction::Enum type;
};

struct DeviceNuclearDataIsotope {
  DeviceReaction *reactions;
};

struct DeviceParticle {
  DeviceParticle &operator=(const MC_Base_Particle &that);
  DeviceParticle &operator=(const MC_Particle &that);
  bool operator==(const MC_Base_Particle &that) const;

  uint64_t identifier;
  int species;
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

  DeviceDomain *domains;
  DeviceMaterial *mats;
  DeviceNuclearDataIsotope *isotopes;

  enum Tallies { ABSORB, COLLISION, FISSION, PRODUCE, SEGMENTS, SCATTER, TALLIES_SIZE };
  long *tallies;

  enum ParticleSizes { PROCESSING, PROCESSED, EXTRAS, PARTICLE_SIZES_SIZE };
  int *particleSizes;
  DeviceParticle *processing;
  DeviceParticle *processed;
  DeviceParticle *extras;

  double nuBar;
  double logLow;
  double divDelta;
  int numGroups;
  int reactionSize;
};

