#include "Device.hh"

#include "cudaUtils.hh"
#include "MonteCarlo.hh"

void Device::init(MonteCarlo &mc)
{
  assert(domain == nullptr);
  const int domainSize = mc.domain.size();
  CHECK(hipHostMalloc(&domain,domainSize*sizeof(DeviceDomain)));

  int csSizeSum = 0;
  for (int i = 0; i < domainSize; i++) csSizeSum += mc.domain[i].cell_state.size();
  
  DeviceCellState *cs = nullptr;
  CHECK(hipHostMalloc(&cs,csSizeSum*sizeof(DeviceCellState)));
  for (int i = 0; i < domainSize; i++) {
    domain[i].cell_state = cs;
    const int csSize = mc.domain[i].cell_state.size();
    for (int j = 0; j < csSize; j++) cs[j] = mc.domain[i].cell_state[j];
    cs += csSize;
  }
}

