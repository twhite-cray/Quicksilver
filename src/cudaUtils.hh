#ifndef CUDAUTILS_HH
#define CUDAUTILS_HH

#if defined(HAVE_CUDA) && defined(HAVE_HIP)
#error "Do not define both HAVE_CUDA and HAVE_HIP"
#endif

#if defined(HAVE_CUDA) || defined(HAVE_OPENMP_TARGET) 
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif

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

#ifdef HAVE_OPENMP_TARGET
    #ifdef USE_OPENMP_NO_GPU
        #define VAR_MEM MemoryControl::AllocationPolicy::HOST_MEM
    #endif
#elif HAVE_CUDA
#else
    #define VAR_MEM MemoryControl::AllocationPolicy::HOST_MEM
#endif

enum ExecutionPolicy{ cpu, gpuWithCUDA, gpuWithOpenMP, gpuWithHip };

inline ExecutionPolicy getExecutionPolicy( int useGPU )
{
    ExecutionPolicy execPolicy = ExecutionPolicy::cpu;

    if( useGPU )
    {
        #if defined (HAVE_CUDA)
        execPolicy = ExecutionPolicy::gpuWithCUDA;
        #elif defined (HAVE_OPENMP_TARGET)
        execPolicy = ExecutionPolicy::gpuWithOpenMP;
        #elif defined (HAVE_HIP)
        execPolicy = ExecutionPolicy::gpuWithHip;
        #endif
    }
    return execPolicy;
}
#endif
