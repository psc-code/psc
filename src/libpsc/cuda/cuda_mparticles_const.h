
#ifndef CUDA_MPARTICLES_CONST_H
#define CUDA_MPARTICLES_CONST_H

#include <cuda_bits.h>
#include <cstdio>

// ----------------------------------------------------------------------
// cuda_mparticles_const
//
// cuda_mparticles parameters in CUDA constant memory
// (to avoid having to pass them in)

#define MAX_KINDS (4)

struct cuda_mparticles_const {
  float dxi[3];
  float b_dxi[3];
  int b_mx[3];
  float dt;
  float dqs;
  float dq[MAX_KINDS];
  float fnqs;
  float fnqxs, fnqys, fnqzs;
};

__constant__ __device__ struct cuda_mparticles_const d_cmprts_const;

static void
cuda_mparticles_const_set(struct cuda_mparticles *cmprts)
{
  struct cuda_mparticles_const c;
  for (int d = 0; d < 3; d++) {
    c.dxi[d] = 1.f / cmprts->dx[d];
    c.b_mx[d] = cmprts->b_mx[d];
    c.b_dxi[d] = cmprts->b_dxi[d];
  }

  c.dt = cmprts->dt;
  c.fnqs  = cmprts->fnqs;
  c.fnqxs = cmprts->dx[0] * c.fnqs / cmprts->dt;
  c.fnqys = cmprts->dx[1] * c.fnqs / cmprts->dt;
  c.fnqzs = cmprts->dx[2] * c.fnqs / cmprts->dt;
  c.dqs = .5f * cmprts->eta * cmprts->dt;

  assert(cmprts->n_kinds <= MAX_KINDS);
  for (int k = 0; k < cmprts->n_kinds; k++) {
    c.dq[k] = c.dqs * cmprts->kind_q[k] / cmprts->kind_m[k];
  }

  cudaError_t ierr = cudaMemcpyToSymbol(d_cmprts_const, &c, sizeof(c)); cudaCheck(ierr);
}

#endif

