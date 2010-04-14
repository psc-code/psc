
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

#include "psc.h"

typedef struct {
  float x, y, z, w;
} float4;

struct d_part {
  float4 *xi4;    // xi , yi , zi , qni_div_mni
  float4 *pxi4;   // pxi, pyi, pzi, qni_wni (if qni==0, then qni_wni = wni)
};

struct psc_cuda {
  float4 *xi4;
  float4 *pxi4;
  struct d_part d_part; // on device
};

void cuda_push_part_yz_a();

// ======================================================================
// CUDA emulation

#include <stdlib.h>
#include <string.h>

#ifndef __unused
#define __unused __attribute__((unused))
#endif

static struct {
  int x, y;
} threadIdx __unused;

static struct {
  int x, y;
} blockIdx __unused;

static struct {
  int x, y;
} blockDim __unused;

#define RUN_KERNEL(dimBlock, dimGrid, func, params) do {		\
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

#define __device__
#define __global__
#define __constant__

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

#define check(a) do { int ierr = a; if (ierr != cudaSuccess) fprintf(stderr, "IERR = %d (%d)\n", ierr, cudaSuccess); assert(ierr == cudaSuccess); } while(0)

// ======================================================================

#endif
