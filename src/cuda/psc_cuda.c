
#include "psc_cuda.h"

#include <stdlib.h>
#include <assert.h>

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

  check(cudaMalloc((void **) &cuda->d_part.xi4, psc.n_part * sizeof(float4)));
  check(cudaMalloc((void **) &cuda->d_part.pxi4, psc.n_part * sizeof(float4)));

  check(cudaMemcpy(cuda->d_part.xi4, xi4, psc.n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(cuda->d_part.pxi4, pxi4, psc.n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
}

static void
cuda_particles_to_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;
  float4 *xi4  = cuda->xi4;
  float4 *pxi4 = cuda->pxi4;

  check(cudaMemcpy(xi4, cuda->d_part.xi4, psc.n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(pxi4, cuda->d_part.pxi4, psc.n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaFree(cuda->d_part.xi4));
  check(cudaFree(cuda->d_part.pxi4));

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

struct psc_ops psc_ops_cuda = {
  .name = "cuda",
  .create                 = cuda_create,
  .destroy                = cuda_destroy,
  .particles_from_fortran = cuda_particles_from_fortran,
  .particles_to_fortran   = cuda_particles_to_fortran,
  .push_part_yz_a         = cuda_push_part_yz_a,
};
