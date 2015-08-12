
#ifndef CUDA_WRAP_H
#define CUDA_WRAP_H

#define check(a) do { int ierr = a; if (ierr != cudaSuccess) fprintf(stderr, "IERR = %d (%d)\n", ierr, cudaSuccess); assert(ierr == cudaSuccess); } while(0)

// ======================================================================

#ifdef __CUDACC__
#define CUDA_DEVICE __device__
#define CUDA_CONSTANT __constant__ __device__
#else
#define CUDA_DEVICE
#define CUDA_CONSTANT
#define __forceinline__
#endif

#ifndef __CUDACC__

// ======================================================================
// CUDA emulation

#include <stdlib.h>
#include <string.h>
#include <math.h>

static struct {
  int x, y;
} threadIdx _mrc_unused;

static struct {
  int x, y;
} blockIdx _mrc_unused;

static struct {
  int x, y;
} blockDim _mrc_unused;

#define RUN_KERNEL(dimGrid, dimBlock, func, params) do {		\
    blockDim.x = dimBlock[0];						\
    blockDim.y = dimBlock[1];						\
    for (blockIdx.y = 0; blockIdx.y < dimGrid[1]; blockIdx.y++) {	\
      for (blockIdx.x = 0; blockIdx.x < dimGrid[0]; blockIdx.x++) {	\
	for (threadIdx.y = 0; threadIdx.y < dimBlock[1]; threadIdx.y++) { \
	  for (threadIdx.x = 0; threadIdx.x < dimBlock[0]; threadIdx.x++) { \
	    func params;						\
	  }								\
	}								\
      }									\
    }									\
  } while (0)

//#define __syncthreads() do {} while (0)

#define __device__
#define __global__
#define __constant__
#define __shared__

enum {
  cudaSuccess,
};

enum {
  cudaMemcpyHostToDevice,
  cudaMemcpyDeviceToHost,
};

static inline int
cudaMalloc(void **p, size_t len)
{
  *p = malloc(len);
  return cudaSuccess;
}

static inline int
cudaMemcpy(void *to, void *from, size_t len, int dir)
{
  memcpy(to, from, len);
  return cudaSuccess;
}

static inline int
cudaFree(void *p)
{
  free(p);
  return cudaSuccess;
}

#define cudaMemcpyToSymbol(to, from, len)		\
  cudaMemcpy(&to, from, len, cudaMemcpyHostToDevice) 
 
typedef struct {
  float x, y, z, w;
} float4;

#define EXTERN_C

static inline float
rsqrtf(float x)
{
  return 1.f / sqrtf(x);
}

static inline int
cuda_nint(float x)
{
  // FIXME?
  return (int)(x + (float)(10.5)) - 10;
}

#else

static bool CUDA_SYNC = true;

static inline void
cuda_sync_if_enabled()
{
  if (CUDA_SYNC) {
    check(cudaThreadSynchronize());
  }
}

#define RUN_KERNEL(dimGrid, dimBlock, func, params) do {	\
    dim3 dG(dimGrid[0], dimGrid[1]);				\
    dim3 dB(dimBlock[0], dimBlock[1]);				\
    func<<<dG, dB>>>params;					\
    cuda_sync_if_enabled();					\
  } while (0)

#define EXTERN_C extern "C"

__device__ static inline int
cuda_nint(float x)
{
  return __float2int_rn(x);
}

#endif

#define fabsr fabsf

#define rsqrtr rsqrtf

static inline float
cuda_int_as_float(int i)
{
  union { int i; float f; } u;
  u.i = i;
  return u.f;
};

static inline int
cuda_float_as_int(float f)
{
  union { int i; float f; } u;
  u.f = f;
  return u.i;
};

#endif
