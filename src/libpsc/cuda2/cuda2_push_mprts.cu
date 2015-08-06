
#include "psc_cuda2.h"

#include "psc_fields_cuda2.h"
#include "psc_particles_as_cuda2.h"

#define EM_CACHE EM_CACHE_NONE
#define CALC_J CALC_J_1VB_VAR1
#define F3_CURR(flds, m, ix,iy,iz) ((float *) flds->data)[0]

#if DIM == DIM_YZ
#define BLOCKBND_X 0
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#elif DIM == DIM_XYZ
#define BLOCKBND_X 2
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#endif

#define BLOCKGSIZE_X (BLOCKSIZE_X + 2 * BLOCKBND_X)
#define BLOCKGSIZE_Y (BLOCKSIZE_Y + 2 * BLOCKBND_Y)
#define BLOCKGSIZE_Z (BLOCKSIZE_Z + 2 * BLOCKBND_Z)

#include "../psc_push_particles/inc_params.c"
#include "../psc_push_particles/inc_push.c"
#include "../psc_push_particles/inc_interpolate.c"

// ======================================================================

#if DIM == DIM_YZ

#define F3_DEV_OFF(fldnr, jx,jy,jz)					\
  ((((fldnr)								\
     *prm.mx[2] + ((jz)-prm.ilg[2]))					\
    *prm.mx[1] + ((jy)-prm.ilg[1])))

#else

#define F3_DEV_OFF(fldnr, jx,jy,jz)					\
  ((((fldnr)								\
     *prm.mx[2] + ((jz)-prm.ilg[2]))					\
    *prm.mx[1] + ((jy)-prm.ilg[1]))					\
   *prm.mx[0] + ((jx)-prm.ilg[0]))

#endif

#define F3_DEV(d_flds, fldnr, jx,jy,jz)		\
  ((d_flds)[F3_DEV_OFF(fldnr, jx,jy,jz)])


__device__ static void
curr_add(real *d_flds, int m, int jx, int jy, int jz, real val, int *ci0)
{
  float *addr = &F3_DEV(d_flds, JXI+m, jx,jy,jz);
  atomicAdd(addr, val);
}

#include "../psc_push_particles/inc_curr.c"

#undef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK (512)

// OPT: precalc offsets into fld_cache (including ci[])
// OPT: use more shmem?

#define LOAD_PARTICLE_POS_(pp, d_xi4, n) do {				\
    float4 _xi4 = d_xi4[n];						\
    (pp).xi[0]         = _xi4.x;					\
    (pp).xi[1]         = _xi4.y;					\
    (pp).xi[2]         = _xi4.z;					\
    (pp).kind_as_float = _xi4.w;					\
} while (0)

#define LOAD_PARTICLE_MOM_(pp, d_pxi4, n) do {				\
    float4 _pxi4 = d_pxi4[n];						\
    (pp).pxi[0]        = _pxi4.x;					\
    (pp).pxi[1]        = _pxi4.y;					\
    (pp).pxi[2]        = _pxi4.z;					\
    (pp).qni_wni       = _pxi4.w;					\
} while (0)

#define STORE_PARTICLE_POS_(pp, d_xi4, n) do {				\
    float4 xi4 = { (pp).xi[0], (pp).xi[1], (pp).xi[2], (pp).kind_as_float }; \
    d_xi4[n] = xi4;							\
} while (0)

#define STORE_PARTICLE_MOM_(pp, d_pxi4, n) do {				\
    float4 pxi4 = { (pp).pxi[0], (pp).pxi[1], (pp).pxi[2], (pp).qni_wni }; \
    d_pxi4[n] = pxi4;							\
} while (0)

// ======================================================================
// field caching

#if EM_CACHE == EM_CACHE_NONE

#define F3_CACHE(flds_em, m, jx, jy, jz)	\
  (F3_DEV(flds_em, m, jx,jy,jz))

#define DECLARE_EM_CACHE(flds_em, d_flds, size, ci0)	\
  real *flds_em = d_flds

#elif EM_CACHE == EM_CACHE_CUDA

#if DIM == DIM_YZ
#define F3_CACHE(flds_em, m, jx, jy, jz)				\
  ((flds_em)[(((m-EX)							\
	       *BLOCKGSIZE_Z + ((jz-ci0[2])-(-BLOCKBND_Z)))		\
	      *BLOCKGSIZE_Y + ((jy-ci0[1])-(-BLOCKBND_Y)))])
#elif DIM == DIM_XYZ
#define F3_CACHE(flds_em, m, jx, jy, jz)				\
  ((flds_em)[((((m-EX)							\
		*BLOCKGSIZE_Z + ((jz-ci0[2])-(-BLOCKBND_Z)))		\
	       *BLOCKGSIZE_Y + ((jy-ci0[1])-(-BLOCKBND_Y)))		\
	      *BLOCKGSIZE_X + ((jx-ci0[0])-(-BLOCKBND_X)))])
#endif

__device__ static void
cache_fields(float *flds_em, float *d_flds, int size, int *ci0)
{
  int ti = threadIdx.x;
  while (ti < n) {
    int tmp = ti;
    int jx = tmp % BLOCKGSIZE_X - BLOCKBND_X;
    tmp /= dims[0];
    int jy = tmp % BLOCKGSIZE_Y - BLOCKBND_Y;
    tmp /= dims[1];
    int jz = tmp % BLOCKGSIZE_Z - BLOCKBND_Z;
    // OPT? currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3_CACHE(flds_em, m, jx+ci0[0],jy+ci0[1] jz+ci0[2]) = 
	F3_DEV(d_flds, m, jx+ci0[0],jy+ci0[1],jz+ci0[2]);
    }
    ti += THREADS_PER_BLOCK;
  }
}

#define DECLARE_EM_CACHE(flds_em, d_flds, size, ci0)	\
  __shared__ real flds_em[6 * BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z];\
  cache_fields(flds_em, d_flds, size, ci0)

#endif

// ----------------------------------------------------------------------
// push_part_one

__device__ static void
push_part_one(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	      real *d_flds, real *flds_em, int ci0[3])
{
  LOAD_PARTICLE_POS_(*prt, d_xi4, n);

  // here we have x^{n+.5}, p^n

  // field interpolation
  real exq, eyq, ezq, hxq, hyq, hzq;
  int lg[3];
  real og[3];
  find_idx_off_1st_rel(prt->xi, lg, og, real(0.));
  INTERPOLATE_1ST_EC(flds_em, exq, eyq, ezq, hxq, hyq, hzq);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  LOAD_PARTICLE_MOM_(*prt, d_pxi4, n);
  int kind = particle_kind(prt);
  real dq = prm.dq_kind[kind];
  push_pxi(prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  STORE_PARTICLE_MOM_(*prt, d_pxi4, n);
}

__device__ static int
find_block_pos_patch(int *block_pos, int *ci0)
{
  block_pos[0] = blockIdx.x;
  block_pos[1] = blockIdx.y;
  block_pos[2] = blockIdx.z % prm.b_mx[2];

#if EM_CACHE == EM_CACHE_CUDA
  ci0[0] = block_pos[0] * BLOCKSIZE_X;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;
#endif

  return blockIdx.z / prm.b_mx[2];
}

__device__ static int
find_bid()
{
  return (blockIdx.z * prm.b_mx[1] + blockIdx.y) * prm.b_mx[0] + blockIdx.x;
}

// ======================================================================
// depositing current

// ----------------------------------------------------------------------
// calc_j

__device__ static void
calc_j(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
       real *d_flds, int p_nr, int bid, int *ci0)
{
  real vxi[3];
  calc_vxi(vxi, prt);

  // position xm at x^(n+.5)
  real h0[3], h1[3];
  real xm[3], xp[3];
  int j[3], k[3];
  
  find_idx_off_pos_1st_rel(prt->xi, j, h0, xm, real(0.));

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, prm.dt);
  STORE_PARTICLE_POS_(*prt, d_xi4, n);

#if 0
  // save block_idx for new particle position at x^(n+1.5)
  unsigned int block_pos_y = __float2int_rd(prt->xi[1] * prm.b_dxi[1]);
  unsigned int block_pos_z = __float2int_rd(prt->xi[2] * prm.b_dxi[2]);
  int nr_blocks = prm.b_mx[1] * prm.b_mx[2];

  int block_idx;
  if (block_pos_y >= prm.b_mx[1] || block_pos_z >= prm.b_mx[2]) {
    block_idx = CUDA_BND_S_OOB;
  } else {
    int bidx = block_pos_z * prm.b_mx[1] + block_pos_y + p_nr * nr_blocks;
    int b_diff = bid - bidx + prm.b_mx[1] + 1;
    int d1 = b_diff % prm.b_mx[1];
    int d2 = b_diff / prm.b_mx[1];
    block_idx = d2 * 3 + d1;
  }
  d_bidx[n] = block_idx;
#endif

  // position xp at x^(n+.5)
  find_idx_off_pos_1st_rel(prt->xi, k, h1, xp, real(0.));

#if DIM == DIM_YZ

  // deposit xm -> xp
  int idiff[3] = { 0, k[1] - j[1], k[2] - j[2] };
  int i[3] = { 0, j[1], j[2] };
  real x[3] = { 0.f, xm[1] - j[1] - real(.5), xm[2] - j[2] - real(.5) };
  //real dx[3] = { 0.f, xp[1] - xm[1], xp[2] - xm[2] };
  real dx[3] = { vxi[0] * prm.dt * prm.dxi[0], xp[1] - xm[1], xp[2] - xm[2] };
  
  real x1 = x[1] * idiff[1];
  real x2 = x[2] * idiff[2];
  int d_first = (abs(dx[2]) * (.5f - x1) >= abs(dx[1]) * (.5f - x2));

  int off[3];
  if (d_first == 0) {
    off[1] = idiff[1];
    off[2] = 0;
  } else {
    off[1] = 0;
    off[2] = idiff[2];
  }

  real dx1[3];
  calc_3d_dx1(dx1, x, dx, off);
  curr_3d_vb_cell(d_flds, i, x, dx1, prt->qni_wni, ci0);
  curr_3d_vb_cell_upd(i, x, dx1, dx, off);
  
  off[1] = idiff[1] - off[1];
  off[2] = idiff[2] - off[2];
  calc_3d_dx1(dx1, x, dx, off);
  curr_3d_vb_cell(d_flds, i, x, dx1, prt->qni_wni, ci0);
  curr_3d_vb_cell_upd(i, x, dx1, dx, off);
    
  curr_3d_vb_cell(d_flds, i, x, dx, prt->qni_wni, ci0);
#endif
}


// ======================================================================

__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(float4 *d_xi4, float4 *d_pxi4,
	      unsigned int *d_off,
	      float *d_flds0, unsigned int size)
{
  int block_pos[3], ci0[3];
  int p, bid;
  p = find_block_pos_patch(block_pos, ci0);
  real *d_flds = d_flds0 + p * size;

  DECLARE_EM_CACHE(flds_em, d_flds, size, ci0);

  bid = find_bid();
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    particle_t prt;
    push_part_one(&prt, n, d_xi4, d_pxi4, d_flds0 + p * size, flds_em, ci0);
    calc_j(&prt, n, d_xi4, d_pxi4, d_flds, p, bid, ci0);
  }

}

// ----------------------------------------------------------------------
// zero_currents

static void
zero_currents(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  unsigned int size = mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_cuda2_real_t *d_flds = mflds_sub->d_flds + p * size * mflds->nr_fields;
    check(cudaMemset(d_flds + JXI * size, 0, 3 * size * sizeof(*d_flds)));
  }
}

// ----------------------------------------------------------------------
// cuda_push_mprts_ab

static void
cuda_push_mprts_ab(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  params_1vb_set(ppsc, mprts, mflds);

  unsigned int fld_size = mflds->nr_fields *
    mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

  zero_currents(mflds);

  int gx, gy, gz;
#if DIM == DIM_YZ
  assert(mprts_sub->b_mx[0] == 1);
#endif
  gx = mprts_sub->b_mx[0];
  gy = mprts_sub->b_mx[1];
  gz = mprts_sub->b_mx[2] * mprts->nr_patches;
  dim3 dimGrid(gx, gy, gz);

  push_mprts_ab<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_sub->d_xi4, mprts_sub->d_pxi4, mprts_sub->d_b_off,
     mflds_sub->d_flds, fld_size);

  cuda_sync_if_enabled();
}

