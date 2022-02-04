#pragma once
#include "MaterialDatabase.hh"
#include "MC_Cell_State.hh"

class MC_Base_Particle;
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
};

struct DeviceNuclearDataIsotope {
  DeviceReaction *reactions;
};

struct DeviceParticle {
  DeviceParticle &operator=(const MC_Base_Particle &that);
  bool operator==(const MC_Base_Particle &that) const;

  uint64_t identifier;
  int species;
};

struct Device {
  Device():
    domains(nullptr),
    mats(nullptr),
    isotopes(nullptr),
    processingSize(nullptr),
    processing(nullptr),
    processed(nullptr)
  {}
  Device(const Device &) = default;
  void init(MonteCarlo &mc);
  void cycleInit(MonteCarlo &mc);
  void cycleFinalize(MonteCarlo &mc);

  DeviceDomain *domains;
  DeviceMaterial *mats;
  DeviceNuclearDataIsotope *isotopes;
  long *numSegments;

  int *processingSize;
  DeviceParticle *processing;
  DeviceParticle *processed;
  int reactionSize;
};

