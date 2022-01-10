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
    return *this;
  }
  double atomFraction;
};

struct DeviceMaterial {
  DeviceIsotope *iso;
};

struct Device {
  Device(): domain(nullptr), mat(nullptr) {}
  Device(const Device &) = default;
  void init(MonteCarlo &mc);
  DeviceDomain *domain;
  DeviceMaterial *mat;
};
