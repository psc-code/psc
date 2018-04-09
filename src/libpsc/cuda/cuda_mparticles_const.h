
#ifndef CUDA_MPARTICLES_CONST_H
#define CUDA_MPARTICLES_CONST_H

#include "cuda_bits.h"
#include <cstdio>

// ----------------------------------------------------------------------
// cuda_mparticles_const
//
// cuda_mparticles parameters in CUDA constant memory
// (to avoid having to pass them in)

#define MAX_KINDS (4)

struct cuda_mparticles_const {
  DParticleIndexer dpi;
  float dxi[3];
  float dt;
  float dqs;
  float dq[MAX_KINDS];
  float q_inv[MAX_KINDS];
  float fnqs;
  float fnqxs, fnqys, fnqzs;
};

__constant__ __device__ struct cuda_mparticles_const d_cmprts_const;

static void
cuda_mparticles_const_set(struct cuda_mparticles *cmprts)
{
  const Grid_t& grid = cmprts->grid_;
  
  struct cuda_mparticles_const c;
  c.dpi = DParticleIndexer(*cmprts);
  for (int d = 0; d < 3; d++) {
    c.dxi[d] = 1.f / grid.domain.dx[d];
  }

  c.dt = grid.dt;
  c.fnqs  = grid.fnqs;
  c.fnqxs = grid.domain.dx[0] * c.fnqs / grid.dt;
  c.fnqys = grid.domain.dx[1] * c.fnqs / grid.dt;
  c.fnqzs = grid.domain.dx[2] * c.fnqs / grid.dt;
  c.dqs = .5f * grid.eta * grid.dt;

  int n_kinds = grid.kinds.size();
  assert(n_kinds <= MAX_KINDS);
  for (int k = 0; k < n_kinds; k++) {
    c.dq[k] = c.dqs * grid.kinds[k].q / grid.kinds[k].m;
    c.q_inv[k] = 1.f / grid.kinds[k].q;
  }

  cudaError_t ierr = cudaMemcpyToSymbol(d_cmprts_const, &c, sizeof(c)); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// find_idx_off_1st

__device__ static void
find_idx_off_1st(const float xi[3], int j[3], float h[3], float shift)
{
  for (int d = 0; d < 3; d++) {
    float pos = xi[d] * d_cmprts_const.dxi[d] + shift;
    j[d] = __float2int_rd(pos);
    h[d] = pos - j[d];
  }
}

// ----------------------------------------------------------------------
// find_idx_off_pos_1st

__device__ static void
find_idx_off_pos_1st(const float xi[3], int j[3], float h[3], float pos[3], float shift)
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * d_cmprts_const.dxi[d] + shift;
    j[d] = __float2int_rd(pos[d]);
    h[d] = pos[d] - j[d];
  }
}

#endif

