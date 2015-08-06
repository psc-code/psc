
#include "psc_cuda2.h"

#include "psc_fields_cuda2.h"
#include "psc_particles_as_cuda2.h"

#include "../psc_push_particles/inc_params.c"
#include "../psc_push_particles/inc_push.c"

#define NO_CACHE

#define MAX_KINDS (4)

#define BND (2)

#define F3_DEV_OFF_XYZ(fldnr, jx,jy,jz)					\
  ((((fldnr)								\
     *prm.mx[2] + ((jz)-prm.ilg[2]))					\
    *prm.mx[1] + ((jy)-prm.ilg[1]))					\
   *prm.mx[0] + ((jx)-prm.ilg[0]))

#define F3_DEV_XYZ(fldnr,jx,jy,jz)		\
  ((d_flds)[F3_DEV_OFF_XYZ(fldnr, jx,jy,jz)])

#define F3_DEV_YZ(fldnr,jy,jz) F3_DEV_XYZ(fldnr, 0,jy,jz)


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

#ifdef NO_CACHE

#define F3_CACHE_YZ(fld_cache, m, jy, jz)	\
  (F3_DEV_YZ(m, jy+ci0[1],jz+ci0[2]))

#define F3_CACHE_XYZ(fld_cache, m, jx, jy, jz)	\
  (F3_DEV_XYZ(m, jx+ci0[0],jy+ci0[1],jz+ci0[2]))

#else

#define F3_CACHE_YZ(fld_cache, m, jy, jz)				\
  ((fld_cache)[(((m-EX)							\
		 *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))			\
		*(BLOCKSIZE_Y + 4) + ((jy)-(-2)))])

#define F3_CACHE_XYZ(fld_cache, m, jx, jy, jz)				\
  ((fld_cache)[((((m-EX)						\
		  *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))			\
		 *(BLOCKSIZE_Y + 4) + ((jy)-(-2)))			\
		*(BLOCKSIZE_Z + 4) + ((jx)-(-2)))])

#endif

// ----------------------------------------------------------------------
// push_part_one

__device__ static void
push_part_one(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	      real *d_flds, real *fld_cache, int ci0[3])
{
  LOAD_PARTICLE_POS_(*prt, d_xi4, n);

  // here we have x^{n+.5}, p^n

  // field interpolation
  real exq, eyq, ezq, hxq, hyq, hzq;
  int lg[3];
  real og[3];
  find_idx_off_1st_rel(prt->xi, lg, og, real(0.));
#if DIM == DIM_YZ
  lg[1] -= ci0[1];
  lg[2] -= ci0[2];
  
  exq = ((1.f - og[1]) * (1.f - og[2]) * F3_CACHE_YZ(fld_cache, EX, lg[1]+0, lg[2]+0) +
	 (      og[1]) * (1.f - og[2]) * F3_CACHE_YZ(fld_cache, EX, lg[1]+1, lg[2]+0) +
	 (1.f - og[1]) * (      og[2]) * F3_CACHE_YZ(fld_cache, EX, lg[1]+0, lg[2]+1) +
	 (      og[1]) * (      og[2]) * F3_CACHE_YZ(fld_cache, EX, lg[1]+1, lg[2]+1));
  eyq = ((1.f - og[2]) * F3_CACHE_YZ(fld_cache, EY, lg[1]  , lg[2]+0) +
	 (      og[2]) * F3_CACHE_YZ(fld_cache, EY, lg[1]  , lg[2]+1));
  ezq = ((1.f - og[1]) * F3_CACHE_YZ(fld_cache, EZ, lg[1]+0, lg[2]  ) +
	 (      og[1]) * F3_CACHE_YZ(fld_cache, EZ, lg[1]+1, lg[2]  ));
  hxq = (F3_CACHE_YZ(fld_cache, HX, lg[1]  , lg[2]  ));
  hyq = ((1.f - og[1]) * F3_CACHE_YZ(fld_cache, HY, lg[1]+0, lg[2]  ) +
	 (      og[1]) * F3_CACHE_YZ(fld_cache, HY, lg[1]+1, lg[2]  ));
  hzq = ((1.f - og[2]) * F3_CACHE_YZ(fld_cache, HZ, lg[1]  , lg[2]+0) +
	 (      og[2]) * F3_CACHE_YZ(fld_cache, HZ, lg[1]  , lg[2]+1));
#elif DIM == DIM_XYZ
  lg[0] -= ci0[0];
  lg[1] -= ci0[1];
  lg[2] -= ci0[2];
  
  exq = ((1.f - og[1]) * (1.f - og[2]) * F3_CACHE_XYZ(fld_cache, EX, lg[0]  ,lg[1]  ,lg[2]  ) +
	 (      og[1]) * (1.f - og[2]) * F3_CACHE_XYZ(fld_cache, EX, lg[0]  ,lg[1]+1,lg[2]  ) +
	 (1.f - og[1]) * (      og[2]) * F3_CACHE_XYZ(fld_cache, EX, lg[0]  ,lg[1]  ,lg[2]+1) +
	 (      og[1]) * (      og[2]) * F3_CACHE_XYZ(fld_cache, EX, lg[0]  ,lg[1]+1,lg[2]+1));
  eyq = ((1.f - og[2]) * (1.f - og[0]) * F3_CACHE_XYZ(fld_cache, EY, lg[0]  ,lg[1]  ,lg[2]  ) +
	 (      og[2]) * (1.f - og[0]) * F3_CACHE_XYZ(fld_cache, EY, lg[0]  ,lg[1]  ,lg[2]+1) +
	 (1.f - og[2]) * (      og[0]) * F3_CACHE_XYZ(fld_cache, EY, lg[0]+1,lg[1]  ,lg[2]  ) +
	 (      og[2]) * (      og[0]) * F3_CACHE_XYZ(fld_cache, EY, lg[0]+1,lg[1]  ,lg[2]+1));
  ezq = ((1.f - og[0]) * (1.f - og[1]) * F3_CACHE_XYZ(fld_cache, EZ, lg[0]  ,lg[1]  ,lg[2]  ) +
	 (      og[0]) * (1.f - og[1]) * F3_CACHE_XYZ(fld_cache, EZ, lg[0]+1,lg[1]  ,lg[2]  ) +
	 (1.f - og[0]) * (      og[1]) * F3_CACHE_XYZ(fld_cache, EZ, lg[0]  ,lg[1]+1,lg[2]  ) +
	 (      og[0]) * (      og[1]) * F3_CACHE_XYZ(fld_cache, EZ, lg[0]+1,lg[1]+1,lg[2]  ));
  hxq = ((1.f - og[0]) * F3_CACHE_XYZ(fld_cache, HX, lg[0]  ,lg[1]  , lg[2]  ) +
	 (      og[0]) * F3_CACHE_XYZ(fld_cache, HX, lg[0]+1,lg[1]  , lg[2]  ));
  hyq = ((1.f - og[1]) * F3_CACHE_XYZ(fld_cache, HY, lg[0]  ,lg[1]  , lg[2]  ) +
	 (      og[1]) * F3_CACHE_XYZ(fld_cache, HY, lg[0]  ,lg[1]+1, lg[2]  ));
  hzq = ((1.f - og[2]) * F3_CACHE_XYZ(fld_cache, HZ, lg[0]  ,lg[1]  , lg[2]  ) +
	 (      og[2]) * F3_CACHE_XYZ(fld_cache, HZ, lg[0]  ,lg[1]  , lg[2]+1));
#endif

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
#if DIM == DIM_YZ
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % prm.b_mx[2];

  ci0[0] = 0;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.y / prm.b_mx[2];

#elif DIM == DIM_XYZ
  block_pos[0] = blockIdx.x;
  block_pos[1] = blockIdx.y;
  block_pos[2] = blockIdx.z % prm.b_mx[2];

  ci0[0] = block_pos[0] * BLOCKSIZE_X;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.z / prm.b_mx[2];
#endif
}

__device__ static int
find_bid()
{
#if DIM == DIM_YZ
  return blockIdx.y * prm.b_mx[1] + blockIdx.x;
#elif DIM == DIM_XYZ
  return (blockIdx.z * prm.b_mx[1] + blockIdx.y) * prm.b_mx[0] + blockIdx.x;
#endif
}

__device__ static void
cache_fields(float *fld_cache, float *d_flds0, int size, int *ci0, int p)
{
  real *d_flds = d_flds0 + p * size;

#if DIM == DIM_YZ
  int n = (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
  int ti = threadIdx.x;
  while (ti < n) {
    int tmp = ti;
    int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
    tmp /= BLOCKSIZE_Y + 4;
    int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
    // OPT? currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3_CACHE_YZ(fld_cache, m, jy, jz) = F3_DEV_YZ(m, jy+ci0[1],jz+ci0[2]);
    }
    ti += THREADS_PER_BLOCK;
  }
#elif DIM == DIM_XYZ
  int n = (BLOCKSIZE_X + 4) * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
  int ti = threadIdx.x;
  while (ti < n) {
    int tmp = ti;
    int jx = tmp % (BLOCKSIZE_X + 4) - 2;
    tmp /= BLOCKSIZE_X + 4;
    int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
    tmp /= BLOCKSIZE_Y + 4;
    int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
    // OPT? currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3_CACHE_XYZ(fld_cache, m, jx, jy, jz) = F3_DEV_XYZ(m, jx+ci0[0],jy+ci0[1],jz+ci0[2]);
    }
    ti += THREADS_PER_BLOCK;
  }
#endif
}

class GCurr {
public:
  real *d_flds;

  __device__ GCurr(real *_d_flds) :
    d_flds(_d_flds)
  {
  }

  __device__ void add(int m, int jy, int jz, float val, int *ci0)
  {
    float *addr = &F3_DEV_YZ(JXI+m, jy+ci0[1],jz+ci0[2]);
    atomicAdd(addr, val);
  }
};

// ======================================================================
// depositing current

// ----------------------------------------------------------------------
// calc_dx1

__device__ static void
calc_dx1(real dx1[3], real x[3], real dx[3], int off[3])
{
  real o1, x1, dx_0, dx_1, dx_2, v0, v1, v2;
  if (off[1] == 0) {
    o1 = off[2];
    x1 = x[2];
    dx_0 = dx[0];
    dx_1 = dx[2];
    dx_2 = dx[1];
  } else {
    o1 = off[1];
    x1 = x[1];
    dx_0 = dx[0];
    dx_1 = dx[1];
    dx_2 = dx[2];
  }
  if ((off[1] == 0 && off[2] == 0) || dx_1 == 0.f) {
    v0 = 0.f;
    v1 = 0.f;
    v2 = 0.f;
  } else {
    v1 = .5f * o1 - x1;
    v2 = dx_2 / dx_1 * v1;
    v0 = dx_0 / dx_1 * v1;
  }
  if (off[1] == 0) {
    dx1[0] = v0;
    dx1[1] = v2;
    dx1[2] = v1;
  } else {
    dx1[0] = v0;
    dx1[1] = v1;
    dx1[2] = v2;
  }
}

// ----------------------------------------------------------------------
// curr_vb_cell

__device__ static void
curr_vb_cell(int i[3], real x[3], real dx[3], real qni_wni,
	     GCurr &scurr, int *ci0)
{
  real xa[3] = { 0.,
		 x[1] + .5f * dx[1],
		 x[2] + .5f * dx[2], };
  if (dx[0] != 0.f) {
    real fnqx = qni_wni * prm.fnqxs;
    real h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
    scurr.add(0, i[1]  , i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), ci0);
    scurr.add(0, i[1]+1, i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), ci0);
    scurr.add(0, i[1]  , i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), ci0);
    scurr.add(0, i[1]+1, i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), ci0);
  }
  if (dx[1] != 0.f) {
    real fnqy = qni_wni * prm.fnqys;
    scurr.add(1, i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]), ci0);
    scurr.add(1, i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]), ci0);
  }
  if (dx[2] != 0.f) {
    real fnqz = qni_wni * prm.fnqzs;
    scurr.add(2, i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]), ci0);
    scurr.add(2, i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]), ci0);
  }
}

// ----------------------------------------------------------------------
// curr_vb_cell_upd

__device__ static void
curr_vb_cell_upd(int i[3], real x[3], real dx1[3], real dx[3], int off[3])
{
  dx[0] -= dx1[0];
  dx[1] -= dx1[1];
  dx[2] -= dx1[2];
  x[1] += dx1[1] - off[1];
  x[2] += dx1[2] - off[2];
  i[1] += off[1];
  i[2] += off[2];
}

// ----------------------------------------------------------------------
// calc_j

__device__ static void
calc_j(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
       GCurr &scurr, int p_nr, int bid, int *ci0)
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
  int i[3] = { 0, j[1] - ci0[1], j[2] - ci0[2] };
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
  calc_dx1(dx1, x, dx, off);
  curr_vb_cell(i, x, dx1, prt->qni_wni, scurr, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
  
  off[1] = idiff[1] - off[1];
  off[2] = idiff[2] - off[2];
  calc_dx1(dx1, x, dx, off);
  curr_vb_cell(i, x, dx1, prt->qni_wni, scurr, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
    
  curr_vb_cell(i, x, dx, prt->qni_wni, scurr, ci0);
#endif
}


// ======================================================================

__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(float4 *d_xi4, float4 *d_pxi4,
	      unsigned int *d_off,
	      float *d_flds0, unsigned int size)
{
  int block_pos[3], ci0[3];
  int p, bid;
  p = find_block_pos_patch(block_pos, ci0);
  bid = find_bid();

  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

#ifdef NO_CACHE
  real *fld_cache = NULL;
#else

#if DIM == DIM_YZ
  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
#elif DIM == DIM_XYZ
  __shared__ real fld_cache[6 * (BLOCKSIZE_X + 4) * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
#endif
  cache_fields(fld_cache, d_flds0, size, ci0, p);

#endif

  GCurr scurr(d_flds0 + p * size);
  
  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    particle_t prt;
    push_part_one(&prt, n, d_xi4, d_pxi4, d_flds0 + p * size, fld_cache, ci0);
    calc_j(&prt, n, d_xi4, d_pxi4, scurr, p, bid, ci0);
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

#if DIM == DIM_YZ
  int gx, gy;
  gx = mprts_sub->b_mx[1];
  gy = mprts_sub->b_mx[2] * mprts->nr_patches;
  dim3 dimGrid(gx, gy);
#elif DIM == DIM_XYZ
  int gx, gy, gz;
  gx = mprts_sub->b_mx[0];
  gy = mprts_sub->b_mx[1];
  gz = mprts_sub->b_mx[2] * mprts->nr_patches;
  dim3 dimGrid(gx, gy, gz);
#endif

  push_mprts_ab<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_sub->d_xi4, mprts_sub->d_pxi4, mprts_sub->d_b_off,
     mflds_sub->d_flds, fld_size);

  cuda_sync_if_enabled();
}

