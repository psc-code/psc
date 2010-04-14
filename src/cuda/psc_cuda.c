
#include "psc_cuda.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>

// ======================================================================

struct {
  int x, y;
} threadIdx;

struct {
  int x, y;
} blockIdx;

struct {
  int x, y;
} blockDim;

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

// ======================================================================

static void
cuda_create()
{
  struct psc_cuda *cuda = malloc(sizeof(*cuda));
  psc.c_ctx = cuda;
}

static void
cuda_destroy()
{
  struct psc_cuda *cuda = psc.c_ctx;
  free(cuda);
}

static void
cuda_particles_from_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;

  float4 *xi4  = calloc(psc.n_part, sizeof(float4));
  float4 *pxi4 = calloc(psc.n_part, sizeof(float4));

  for (int i = 0; i < psc.n_part; i++) {
    real qni = psc.f_part[i].qni;
    real wni = psc.f_part[i].wni;
    real qni_div_mni = qni / psc.f_part[i].mni;
    real qni_wni;
    if (qni != 0.) {
      qni_wni = qni * wni;
    } else {
      qni_wni = wni;
    }

    xi4[i].x  = psc.f_part[i].xi;
    xi4[i].y  = psc.f_part[i].yi;
    xi4[i].z  = psc.f_part[i].zi;
    xi4[i].w  = qni_div_mni;
    pxi4[i].x = psc.f_part[i].pxi;
    pxi4[i].y = psc.f_part[i].pyi;
    pxi4[i].z = psc.f_part[i].pzi;
    pxi4[i].w = qni_wni;
  }

  cuda->xi4 = xi4;
  cuda->pxi4 = pxi4;
}

static void
cuda_particles_to_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;
  float4 *xi4  = cuda->xi4;
  float4 *pxi4 = cuda->pxi4;

  for (int i = 0; i < psc.n_part; i++) {
    f_real qni_div_mni = xi4[i].w;
    f_real qni_wni = pxi4[i].w;
    f_real qni, mni, wni;
    if (qni_div_mni == 0.) {
      qni = 0.;
      wni = qni_wni;
      mni = -1.;
      assert(0); // can't recover the mass of a neutral particle
    } else {
      qni = qni_div_mni > 0 ? 1. : -1.;
      mni = qni / qni_div_mni;
      wni = qni_wni / qni;
    }

    psc.f_part[i].xi  = xi4[i].x;
    psc.f_part[i].yi  = xi4[i].y;
    psc.f_part[i].zi  = xi4[i].z;
    psc.f_part[i].pxi = pxi4[i].x;
    psc.f_part[i].pyi = pxi4[i].y;
    psc.f_part[i].pzi = pxi4[i].z;
    psc.f_part[i].qni = qni;
    psc.f_part[i].mni = mni;
    psc.f_part[i].wni = wni;
  }

  free(xi4);
  free(pxi4);
}

static float _dt;

static void set_constants()
{
  _dt = psc.dt;
}

static void
push_part_yz_a_one(int n, float4 *d_xi4, float4 *d_pxi4)
{
  float4 xi4  = d_xi4[n];
  float4 pxi4 = d_pxi4[n];
  
  float root = 1. / sqrt(1. + sqr(pxi4.x) + sqr(pxi4.y) + sqr(pxi4.z));
  float vyi = pxi4.y * root;
  float vzi = pxi4.z * root;
  
  xi4.y += vyi * .5 * _dt;
  xi4.z += vzi * .5 * _dt;

  d_xi4[n] = xi4;
}

static void
push_part_yz_a(int n_part, float4 *d_xi4, float4 *d_pxi4, int stride)
{
  int n = threadIdx.x + blockDim.x * blockIdx.x;

  while (n < n_part) {
    push_part_yz_a_one(n, d_xi4, d_pxi4);
    n += stride;
  }
}

static void
cuda_push_part_yz_a()
{
  struct psc_cuda *cuda = psc.c_ctx;

  set_constants();

  const int threadsPerBlock = 128;
  const int gridSize = 256;
  int dimGrid[2]  = { gridSize, 1 };
  int dimBlock[2] = { threadsPerBlock, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_yz_a, (psc.n_part, cuda->xi4, cuda->pxi4,
			      gridSize * threadsPerBlock));

  for (int n = 0; n < psc.n_part; n++) {
  }
}

struct psc_ops psc_ops_cuda = {
  .name = "cuda",
  .create                 = cuda_create,
  .destroy                = cuda_destroy,
  .particles_from_fortran = cuda_particles_from_fortran,
  .particles_to_fortran   = cuda_particles_to_fortran,
  .push_part_yz_a         = cuda_push_part_yz_a,
};
