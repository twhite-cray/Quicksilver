#ifndef CUDAUTILS_HH
#define CUDAUTILS_HH

#include <hip/hip_runtime.h>
__attribute__((unused))
static void checkHip(const hipError_t err, const char *const file, const int line)
{
  if (err == hipSuccess) return;
  fprintf(stderr,"HIP ERROR AT LINE %d OF FILE '%s': %s %s\n",line,file,hipGetErrorName(err),hipGetErrorString(err));
  fflush(stderr);
  exit(err);
}
#define CHECK(X) checkHip(X,__FILE__,__LINE__)

#define VAR_MEM MemoryControl::AllocationPolicy::HOST_MEM

template <typename T>
__host__ __device__ T atomicFetchAdd(T *const p, const T x)
{
#ifdef __HIP_DEVICE_COMPILE__
  return atomicAdd(p,x);
#else
  return __atomic_fetch_add(p,x,__ATOMIC_RELAXED);
#endif
}

template <typename T>
void hipCalloc(const size_t n, T *__restrict__ &p)
{
  const size_t bytes = n*sizeof(T);
  CHECK(hipMalloc((void **)&p,bytes));
  CHECK(hipMemset(p,0,bytes));
  CHECK(hipDeviceSynchronize());
}

#endif
