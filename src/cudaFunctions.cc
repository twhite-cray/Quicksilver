#include "cudaFunctions.hh"
#include "cudaUtils.hh"

__global__ static void warmup(const bool stop)
{
  if (stop) abort();
}

void warmup_kernel()
{
  warmup<<<1,1>>>(false);
  CHECK(hipDeviceSynchronize());
}
