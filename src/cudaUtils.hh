#ifndef CUDAUTILS_HH
#define CUDAUTILS_HH

#if defined(HAVE_CUDA) || defined(HAVE_OPENMP_TARGET) 
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif

#ifdef HAVE_OPENMP_TARGET
    #ifdef USE_OPENMP_NO_GPU
        #define VAR_MEM MemoryControl::AllocationPolicy::HOST_MEM
    #else
        #define VAR_MEM MemoryControl::AllocationPolicy::UVM_MEM
        #define HAVE_UVM
    #endif
#elif HAVE_CUDA
    #define VAR_MEM MemoryControl::AllocationPolicy::UVM_MEM
    #define HAVE_UVM
#else
    #define VAR_MEM MemoryControl::AllocationPolicy::HOST_MEM
#endif

enum ExecutionPolicy{ cpu, gpuWithCUDA, gpuWithOpenMP };

inline ExecutionPolicy getExecutionPolicy( int useGPU )
{
    ExecutionPolicy execPolicy = ExecutionPolicy::cpu;

    if( useGPU )
    {
        #if defined (HAVE_CUDA)
        execPolicy = ExecutionPolicy::gpuWithCUDA;
        #elif defined (HAVE_OPENMP_TARGET)
        execPolicy = ExecutionPolicy::gpuWithOpenMP;
        #endif
    }
    return execPolicy;
}

#ifdef HAVE_CUDA
#include <cstdio>
#include <cstdlib>
__attribute__((unused))
static void checkCuda(const cudaError_t err, const char *const file, const int line)
{
  if (err == cudaSuccess) return;
  fprintf(stderr,"CUDA ERROR AT LINE %d OF FILE '%s': %s %s\n",line,file,cudaGetErrorName(err),cudaGetErrorString(err));
  fflush(stderr);
  exit(err);
}
#define CHECK(X) checkCuda(X,__FILE__,__LINE__)
#endif

#endif
