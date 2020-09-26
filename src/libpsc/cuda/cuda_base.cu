
#include "PscConfig.h"

#ifdef PSC_HAVE_RMM
#include <rmm/mr/device/cnmem_memory_resource.hpp>
#include <rmm/mr/device/cuda_memory_resource.hpp>
#include <rmm/mr/device/default_memory_resource.hpp>
#include <rmm/thrust_rmm_allocator.h>
#endif

#include <cstdio>
#include <cassert>
#include <cuda_bits.h>
#include <mrc_common.h>

void cuda_base_init(void)
{
  static bool first_time = true;
  if (!first_time)
    return;

  first_time = false;

#ifdef PSC_HAVE_RMM
  static rmm::mr::cnmem_memory_resource pool_mr;
  rmm::mr::set_default_resource(&pool_mr);
#endif

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);

  // This function call returns 0 if there are no CUDA capable devices.
  if (deviceCount == 0) {
    printf("There is no device supporting CUDA\n");
    return;
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return;
  }

  for (int dev = 0; dev < deviceCount; ++dev) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    if (dev == 0) {
      // This function call returns 9999 for both major & minor fields, if no
      // CUDA capable devices are present
      if (deviceProp.major == 9999 && deviceProp.minor == 9999)
        printf("There is no device supporting CUDA.\n");
      else if (deviceCount == 1)
        printf("There is 1 device supporting CUDA\n");
      else
        printf("There are %d devices supporting CUDA\n", deviceCount);
    }
    printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
    printf("  CUDA Capability Major revision number:         %d\n",
           deviceProp.major);
    printf("  CUDA Capability Minor revision number:         %d\n",
           deviceProp.minor);
    printf("  Total amount of global memory:                 %lu bytes\n",
           deviceProp.totalGlobalMem);
#if CUDART_VERSION >= 2000
    printf("  Number of multiprocessors:                     %d\n",
           deviceProp.multiProcessorCount);
    printf("  Number of cores:                               %d\n",
           8 * deviceProp.multiProcessorCount);
#endif
    printf("  Total amount of constant memory:               %lu bytes\n",
           deviceProp.totalConstMem);
    printf("  Total amount of shared memory per block:       %lu bytes\n",
           deviceProp.sharedMemPerBlock);
    printf("  Total number of registers available per block: %d\n",
           deviceProp.regsPerBlock);
    printf("  Warp size:                                     %d\n",
           deviceProp.warpSize);
    printf("  Maximum number of threads per block:           %d\n",
           deviceProp.maxThreadsPerBlock);
    printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
           deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
           deviceProp.maxThreadsDim[2]);
    printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
           deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
           deviceProp.maxGridSize[2]);
    printf("  Maximum memory pitch:                          %lu bytes\n",
           deviceProp.memPitch);
    printf("  Texture alignment:                             %lu bytes\n",
           deviceProp.textureAlignment);
    printf("  Clock rate:                                    %.2f GHz\n",
           deviceProp.clockRate * 1e-6f);
#if CUDART_VERSION >= 2000
    printf("  Concurrent copy and execution:                 %s\n",
           deviceProp.deviceOverlap ? "Yes" : "No");
#endif
#if CUDART_VERSION >= 2020
    printf("  Run time limit on kernels:                     %s\n",
           deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
    printf("  Integrated:                                    %s\n",
           deviceProp.integrated ? "Yes" : "No");
    printf("  Support host page-locked memory mapping:       %s\n",
           deviceProp.canMapHostMemory ? "Yes" : "No");
    printf(
      "  Compute mode:                                  %s\n",
      deviceProp.computeMode == cudaComputeModeDefault
        ? "Default (multiple host threads can use this device simultaneously)"
        : deviceProp.computeMode == cudaComputeModeExclusive
            ? "Exclusive (only one host thread at a time can use this device)"
            : deviceProp.computeMode == cudaComputeModeProhibited
                ? "Prohibited (no host thread can use this device)"
                : "Unknown");
#endif
  }
}

// ----------------------------------------------------------------------
// myCudaMalloc

void* myCudaMalloc(size_t len)
{
  void* rv;
  cudaError_t ierr;
  ierr = cudaMalloc(&rv, len);
  cudaCheck(ierr);
  // mprintf("myCudaMalloc %ld\n", len);

  return rv;
}

// ----------------------------------------------------------------------
// myCudaFree

void myCudaFree(void* ptr)
{
  cudaError_t ierr;
  ierr = cudaFree(ptr);
  cudaCheck(ierr);
  // mprintf("myCudaFree\n");
}
