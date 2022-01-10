#pragma once
#include "MC_Cell_State.hh"

class MonteCarlo;

struct DeviceCellState {
  DeviceCellState &operator=(const MC_Cell_State &that)
  {
    _material = that._material;
    return *this;
  }
  int _material;
};

struct DeviceDomain {
  DeviceCellState *cell_state;
};

struct Device {
  Device(): domain(nullptr) {}
  void init(MonteCarlo &mc);
  DeviceDomain *domain;
};
