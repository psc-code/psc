
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
// struct d_particle

struct d_particle {
  float xi[3];
  float kind_as_float;
  float pxi[3];
  float qni_wni;
};

#define LOAD_PARTICLE_POS(pp, d_xi4, n) do {				\
    float4 _xi4 = d_xi4[n];						\
    (pp).xi[0]         = _xi4.x;					\
    (pp).xi[1]         = _xi4.y;					\
    (pp).xi[2]         = _xi4.z;					\
    (pp).kind_as_float = _xi4.w;					\
} while (0)

#define LOAD_PARTICLE_MOM(pp, d_pxi4, n) do {				\
    float4 _pxi4 = d_pxi4[n];						\
    (pp).pxi[0]        = _pxi4.x;					\
    (pp).pxi[1]        = _pxi4.y;					\
    (pp).pxi[2]        = _pxi4.z;					\
    (pp).qni_wni       = _pxi4.w;					\
} while (0)

#define STORE_PARTICLE_POS(pp, d_xi4, n) do {				\
    float4 xi4 = { (pp).xi[0], (pp).xi[1], (pp).xi[2], (pp).kind_as_float }; \
    d_xi4[n] = xi4;							\
} while (0)

#define STORE_PARTICLE_MOM(pp, d_pxi4, n) do {				\
    float4 pxi4 = { (pp).pxi[0], (pp).pxi[1], (pp).pxi[2], (pp).qni_wni }; \
    d_pxi4[n] = pxi4;							\
} while (0)

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

