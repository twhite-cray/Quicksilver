#ifndef CUDAUTILS_HH
#define CUDAUTILS_HH

#if defined(HAVE_HIP)
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
#endif

#define VAR_MEM MemoryControl::AllocationPolicy::HOST_MEM

enum ExecutionPolicy{ cpu, gpuWithOpenMP, gpuWithHip };

inline ExecutionPolicy getExecutionPolicy( int useGPU )
{
    ExecutionPolicy execPolicy = ExecutionPolicy::cpu;

    if( useGPU )
    {
        #if defined (HAVE_OPENMP_TARGET)
        execPolicy = ExecutionPolicy::gpuWithOpenMP;
        #elif defined (HAVE_HIP)
        execPolicy = ExecutionPolicy::gpuWithHip;
        #endif
    }
    return execPolicy;
}
#endif
