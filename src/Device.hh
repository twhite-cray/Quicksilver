#pragma once
#include "MaterialDatabase.hh"
#include "MC_Cell_State.hh"

class MC_Base_Particle;
class MonteCarlo;

struct DeviceCellState {
  static constexpr int numFacets = 24;
  DeviceCellState &operator=(const MC_Cell_State &that)
  {
    cellNumberDensity = that._cellNumberDensity;
    material = that._material;
    return *this;
  }
  double *totals;
  double4 facets[numFacets];
  double cellNumberDensity;
  int3 points[numFacets];
  int material;
};

struct DeviceDomain {
  DeviceCellState *cellStates;
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
  uint64_t identifier;
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

  DeviceDomain *domains;
  DeviceMaterial *mats;
  DeviceNuclearDataIsotope *isotopes;

  int *processingSize;
  DeviceParticle *processing;
  DeviceParticle *processed;
};

