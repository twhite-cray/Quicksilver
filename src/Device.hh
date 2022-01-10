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
  double cellNumberDensity;
  int material;
};

struct DeviceDomain {
  DeviceCellState *cellState;
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
  DeviceIsotope *iso;
};

struct DeviceNuclearDataIsotope {
  double *totalCrossSections;
};

struct Device {
  Device(): domain(nullptr), mat(nullptr) {}
  Device(const Device &) = default;
  void init(MonteCarlo &mc);
  DeviceDomain *domain;
  DeviceMaterial *mat;
  DeviceNuclearDataIsotope *isotopes;
};
