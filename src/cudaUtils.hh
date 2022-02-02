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

#endif
