
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

#include "psc_fields_cuda.h"

#include <assert.h>
#include <math.h>
#include <psc.h>

#define rsqrtr rsqrtf

#ifndef __CUDACC__

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
nint(real x)
{
  // FIXME?
  return (int)(x + real(10.5)) - 10;
}

#else

#define RUN_KERNEL(dimGrid, dimBlock, func, params) do {	\
    dim3 dG(dimGrid[0], dimGrid[1]);				\
    dim3 dB(dimBlock[0], dimBlock[1]);				\
    func<<<dG, dB>>>params;					\
    check(cudaThreadSynchronize()); /* FIXME */			\
  } while (0)

#define EXTERN_C extern "C"

__device__ static inline int
nint(real x)
{
  return __float2int_rn(x);
}

#endif

// ======================================================================

#define check(a) do { int ierr = a; if (ierr != cudaSuccess) fprintf(stderr, "IERR = %d (%d)\n", ierr, cudaSuccess); assert(ierr == cudaSuccess); } while(0)

// ======================================================================

typedef float particle_cuda_real_t;

#define MPI_PARTICLES_CUDA_REAL MPI_FLOAT

struct d_part {
  float4 *xi4;    // xi , yi , zi , qni_div_mni
  float4 *pxi4;   // pxi, pyi, pzi, qni_wni (if qni==0, then qni_wni = wni)
  int *offsets;   // particles per block are
                  // are at indices offsets[block] .. offsets[block+1]-1
};

typedef struct {
  struct d_part h_part; // all particles, on host
  struct d_part d_part; // all particles, on device
} particles_cuda_t;

struct psc_cuda {
  particles_cuda_t p;
  int nr_blocks;        // number of blocks
  int b_mx[3];          // number of blocks by direction
};

EXTERN_C void cuda_push_part_yz_a();
EXTERN_C void cuda_push_part_yz_b();
EXTERN_C void cuda_push_part_yz_b2();
EXTERN_C void __cuda_particles_from_fortran(struct psc_cuda *cuda);
EXTERN_C void __cuda_particles_to_fortran(struct psc_cuda *cuda);
EXTERN_C void __cuda_fields_from_fortran(fields_cuda_t *pf);
EXTERN_C void __cuda_fields_to_fortran(fields_cuda_t *pf);

struct d_particle {
  real xi[3];
  real qni_div_mni;
  real pxi[3];
  real qni_wni;
};

#define THREADS_PER_BLOCK 128
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 1
#define BLOCKSIZE_Z 1

#define LOAD_PARTICLE(pp, d_p, n) do {					\
    (pp).xi[0]       = d_p.xi4[n].x;					\
    (pp).xi[1]       = d_p.xi4[n].y;					\
    (pp).xi[2]       = d_p.xi4[n].z;					\
    (pp).qni_div_mni = d_p.xi4[n].w;					\
    (pp).pxi[0]      = d_p.pxi4[n].x;					\
    (pp).pxi[1]      = d_p.pxi4[n].y;					\
    (pp).pxi[2]      = d_p.pxi4[n].z;					\
    (pp).qni_wni     = d_p.pxi4[n].w;					\
} while (0)

#define STORE_PARTICLE_POS(pp, d_p, n) do {				\
    d_p.xi4[n].x = (pp).xi[0];						\
    d_p.xi4[n].y = (pp).xi[1];						\
    d_p.xi4[n].z = (pp).xi[2];						\
    d_p.xi4[n].w = (pp).qni_div_mni;					\
} while (0)

#define STORE_PARTICLE_MOM(pp, d_p, n) do {				\
    d_p.pxi4[n].x = (pp).pxi[0];					\
    d_p.pxi4[n].y = (pp).pxi[1];					\
    d_p.pxi4[n].z = (pp).pxi[2];					\
    d_p.pxi4[n].w = (pp).qni_wni;					\
} while (0)

#endif
