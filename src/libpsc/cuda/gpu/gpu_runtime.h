
#ifdef GTENSOR_DEVICE_HIP
#include <hip/hip_runtime.h>
#else
#include <cuda_runtime.h>

#define hipError_t cudaError_t
#define hipSuccess cudaSuccess
#define hipGetLastError cudaGetLastError
#define hipDeviceSynchronize cudaDeviceSynchronize
#define hipGetErrorName cudaGetErrorName

#define hipGetDeviceCount cudaGetDeviceCount
#define hipDeviceProp_t cudaDeviceProp
#define hipGetDeviceProperties cudaGetDeviceProperties
#define hipComputeModeDefault cudaComputeModeDefault
#define hipComputeModeExclusive cudaComputeModeExclusive
#define hipComputeModeProhibited cudaComputeModeProhibited

#endif
