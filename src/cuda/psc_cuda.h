
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

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


struct d_part {
  float4 *xi4;    // xi , yi , zi , qni_div_mni
  float4 *pxi4;   // pxi, pyi, pzi, qni_wni (if qni==0, then qni_wni = wni)
};

struct psc_cuda {
  float4 *xi4;
  float4 *pxi4;
  float *flds;
  struct d_part d_part; // all particles, on device
  float *d_flds;        // all fields, on device
};

EXTERN_C void cuda_push_part_yz_a();
EXTERN_C void cuda_push_part_yz_b();
EXTERN_C void __cuda_particles_from_fortran(struct psc_cuda *cuda);
EXTERN_C void __cuda_particles_to_fortran(struct psc_cuda *cuda);
EXTERN_C void __cuda_fields_from_fortran(struct psc_cuda *cuda);
EXTERN_C void __cuda_fields_to_fortran(struct psc_cuda *cuda);

struct d_particle {
  real xi[3];
  real qni_div_mni;
  real pxi[3];
  real qni_wni;
};

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

// ----------------------------------------------------------------------
// macros to access fields from CUDA

#define F3_OFF(fldnr, jx,jy,jz)						\
  ((((fldnr)								\
     *d_mx[2] + ((jz)-d_iglo[2]))					\
    *d_mx[1] + ((jy)-d_iglo[1]))					\
   *d_mx[0] + ((jx)-d_iglo[0]))

#if 1

#define F3(fldnr, jx,jy,jz) \
  (d_flds)[F3_OFF(fldnr, jx,jy,jz)]

#else

#define F3(fldnr, jx,jy,jz)						\
  (*({int off = F3_OFF(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS * d_mx[0] * d_mx[1] * d_mx[2]);		\
      &(d_flds[off]);							\
    }))

#endif

// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define CF3_OFF(fldnr, jx,jy,jz)					\
  ((((fldnr)								\
     *psc.img[2] + ((jz)-psc.ilg[2]))					\
    *psc.img[1] + ((jy)-psc.ilg[1]))					\
   *psc.img[0] + ((jx)-psc.ilg[0]))

#define CF3(fldnr, jx,jy,jz) \
  (cuda->flds)[CF3_OFF(fldnr, jx,jy,jz)]

// ----------------------------------------------------------------------
// macros to access fields from CUDA

#define F3_OFF(fldnr, jx,jy,jz)						\
  ((((fldnr)								\
     *d_mx[2] + ((jz)-d_iglo[2]))					\
    *d_mx[1] + ((jy)-d_iglo[1]))					\
   *d_mx[0] + ((jx)-d_iglo[0]))

#if 1

#define F3(fldnr, jx,jy,jz) \
  (d_flds)[F3_OFF(fldnr, jx,jy,jz)]

#else

#define F3(fldnr, jx,jy,jz)						\
  (*({int off = F3_OFF(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS * d_mx[0] * d_mx[1] * d_mx[2]);		\
      &(d_flds[off]);							\
    }))

#endif

// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define CF3_OFF(fldnr, jx,jy,jz)					\
  ((((fldnr)								\
     *psc.img[2] + ((jz)-psc.ilg[2]))					\
    *psc.img[1] + ((jy)-psc.ilg[1]))					\
   *psc.img[0] + ((jx)-psc.ilg[0]))

#define CF3(fldnr, jx,jy,jz) \
  (cuda->flds)[CF3_OFF(fldnr, jx,jy,jz)]

#endif
