
#include "cuda_mparticles.h"

#include <cstdio>
#include <cassert>

#define cudaCheck(ierr) do {						\
    if (ierr != cudaSuccess)						\
      fprintf(stderr, "IERR = %d (%s)\n", ierr, cudaGetErrorName(ierr)); \
    assert(ierr == cudaSuccess);					\
  } while(0)

// ----------------------------------------------------------------------
// cuda_mparticles_create

struct cuda_mparticles *
cuda_mparticles_create()
{
  struct cuda_mparticles *cmprts = 
    (struct cuda_mparticles *) calloc(1, sizeof(*cmprts));

  return cmprts;
}

// ----------------------------------------------------------------------
// cuda_mparticles_destroy

void
cuda_mparticles_destroy(struct cuda_mparticles *cmprts)
{
  free(cmprts);
}

// ----------------------------------------------------------------------
// cuda_mparticles_alloc

void
cuda_mparticles_alloc(struct cuda_mparticles *cmprts, int nr_alloced)
{
  cudaError_t ierr;

  ierr = cudaMalloc((void **) &cmprts->d_xi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_pxi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_alt_xi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_alt_pxi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_bidx, nr_alloced * sizeof(unsigned int)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_id, nr_alloced * sizeof(unsigned int)); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_dealloc

void
cuda_mparticles_dealloc(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  ierr = cudaFree(cmprts->d_xi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_alt_xi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_alt_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_bidx); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_id); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_swap_alt

void
cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts)
{
  float4 *tmp_xi4 = cmprts->d_alt_xi4;
  float4 *tmp_pxi4 = cmprts->d_alt_pxi4;
  cmprts->d_alt_xi4 = cmprts->d_xi4;
  cmprts->d_alt_pxi4 = cmprts->d_pxi4;
  cmprts->d_xi4 = tmp_xi4;
  cmprts->d_pxi4 = tmp_pxi4;
}

