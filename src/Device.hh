#pragma once
#include "MaterialDatabase.hh"
#include "MC_Cell_State.hh"

class MonteCarlo;

struct DeviceCellState {
  DeviceCellState &operator=(const MC_Cell_State &that)
  {
    cellNumberDensity = that._cellNumberDensity;
    material = that._material;
    return *this;
  }
  double *totals;
  double cellNumberDensity;
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
};

struct DeviceReaction {
  double *crossSections;
};

struct DeviceNuclearDataIsotope {
  DeviceReaction *reactions;
};

struct Device {
  Device(): domains(nullptr), mats(nullptr), isotopes(nullptr) {}
  Device(const Device &) = default;
  void init(MonteCarlo &mc);
  void cycleInit(MonteCarlo &mc);
  DeviceDomain *domains;
  DeviceMaterial *mats;
  DeviceNuclearDataIsotope *isotopes;
};
