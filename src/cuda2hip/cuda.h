#pragma once

#include <hip/hip_runtime.h>

#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaError_t hipError_t
#define cudaFree hipFree
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaGetErrorName hipGetErrorName
#define cudaGetErrorString hipGetErrorString
#define cudaMemAttachGlobal NULL
#define cudaPeekAtLastError hipPeekAtLastError
#define cudaSetDevice hipSetDevice
#define cudaSuccess hipSuccess

#define cudaMallocManaged(A,B,C) hipMalloc((A),(B))
