
#include "psc_cuda2.h"

#include "psc_fields_cuda2.h"
#include "psc_particles_cuda2.h"

struct d_particle {
  real xi[3];
  real kind_as_float;
  real pxi[3];
  real qni_wni;
};

EXTERN_C void psc_mparticles_cuda_copy_to_dev(struct psc_mparticles *mprts);

#define MAX_KINDS (4)

struct cuda_params {
  real dt;
  real dxi[3];
  real b_dxi[3];
  real dqs;
  real fnqs;
  real fnqxs, fnqys, fnqzs;
  int mx[3];
  int ilg[3];
  int b_mx[3];
  real dq[MAX_KINDS];
};

#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#undef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK (512)

// OPT: precalc offsets into fld_cache (including ci[])
// OPT: use more shmem?

#define LOAD_PARTICLE_(pp, d_xi4, d_pxi4, n) do {			\
    float4 xi4 = d_xi4[n];						\
    (pp).xi[0]         = xi4.x;						\
    (pp).xi[1]         = xi4.y;						\
    (pp).xi[2]         = xi4.z;						\
    (pp).kind_as_float = xi4.w;						\
    float4 pxi4 = d_pxi4[n];						\
    (pp).pxi[0]        = pxi4.x;					\
    (pp).pxi[1]        = pxi4.y;					\
    (pp).pxi[2]        = pxi4.z;					\
    (pp).qni_wni       = pxi4.w;					\
} while (0)

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

#if 0
#define STORE_PARTICLE_POS_(pp, d_xi4, n) do {				\
    d_xi4[n].x = (pp).xi[0];						\
    d_xi4[n].y = (pp).xi[1];						\
    d_xi4[n].z = (pp).xi[2];						\
    d_xi4[n].w = (pp).kind_as_float;					\
} while (0)
#else
#define STORE_PARTICLE_POS_(pp, d_xi4, n) do {				\
    float4 xi4 = { (pp).xi[0], (pp).xi[1], (pp).xi[2], (pp).kind_as_float }; \
    d_xi4[n] = xi4;							\
} while (0)
#endif

#if 0
#define STORE_PARTICLE_MOM_(pp, d_pxi4, n) do {				\
    d_pxi4[n].x = (pp).pxi[0];						\
    d_pxi4[n].y = (pp).pxi[1];						\
    d_pxi4[n].z = (pp).pxi[2];						\
    d_pxi4[n].w = (pp).qni_wni;						\
} while (0)
#else
#define STORE_PARTICLE_MOM_(pp, d_pxi4, n) do {				\
    float4 pxi4 = { (pp).pxi[0], (pp).pxi[1], (pp).pxi[2], (pp).qni_wni }; \
    d_pxi4[n] = pxi4;							\
} while (0)
#endif

// ======================================================================

// FIXME -> common.c

__device__ static void
find_idx_off_1st(const real xi[3], int j[3], real h[3], real shift,
		 struct cuda_params prm)
{
  for (int d = 0; d < 3; d++) {
    real pos = xi[d] * prm.dxi[d] + shift;
    j[d] = __float2int_rd(pos);
    h[d] = pos - j[d];
  }
}

__device__ static void
find_idx_off_pos_1st(const real xi[3], int j[3], real h[3], real pos[3], real shift,
		     struct cuda_params prm)
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * prm.dxi[d] + shift;
    j[d] = __float2int_rd(pos[d]);
    h[d] = pos[d] - j[d];
  }
}

static __constant__ __device__ float c_dqs[4]; // FIXME hardcoded

static void
set_consts(struct cuda_params *prm)
{
  check(cudaMemcpyToSymbol(c_dqs, prm->dq, sizeof(c_dqs)));
}

static void
_set_params(struct cuda_params *prm, struct psc *psc,
	    struct psc_mfields *mflds,
	    struct psc_mparticles *mprts_cuda)
{
  prm->dt = psc->dt;
  for (int d = 0; d < 3; d++) {
    prm->dxi[d] = 1.f / ppsc->patch[0].dx[d];
  }

  prm->dqs    = .5f * psc->coeff.eta * psc->dt;
  prm->fnqs   = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
  prm->fnqxs  = psc->patch[0].dx[0] * prm->fnqs / psc->dt;
  prm->fnqys  = psc->patch[0].dx[1] * prm->fnqs / psc->dt;
  prm->fnqzs  = psc->patch[0].dx[2] * prm->fnqs / psc->dt;
  assert(psc->nr_kinds <= MAX_KINDS);
  for (int k = 0; k < psc->nr_kinds; k++) {
    prm->dq[k] = prm->dqs * psc->kinds[k].q / psc->kinds[k].m;
  }

  if (mprts_cuda && mprts_cuda->nr_patches > 0) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts_cuda, 0);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);
    for (int d = 0; d < 3; d++) {
      prm->b_mx[d] = prts_cuda->b_mx[d];
      prm->b_dxi[d] = prts_cuda->b_dxi[d];
    }
  }

  if (mflds) {
    struct psc_mfields_cuda2 * mflds_sub = psc_mfields_cuda2(mflds);
    for (int d = 0; d < 3; d++) {
      prm->mx[d] = mflds_sub->im[d];
      prm->ilg[d] = mflds_sub->ib[d];
      if (d > 0) {
	assert(mflds_sub->ib[d] == -BND);
      } else {
	assert(mflds_sub->im[d] == 1);
      }
    }
  }
}

static void
_free_params(struct cuda_params *prm)
{
}

// ======================================================================
// field caching

#define F3_CACHE(fld_cache, m, jy, jz)					\
  ((fld_cache)[(((m-EX)							\
		 *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))			\
		*(BLOCKSIZE_Y + 4) + ((jy)-(-2)))])

// ----------------------------------------------------------------------
// push_xi
//
// advance position using velocity

__device__ static void
push_xi(struct d_particle *p, const real vxi[3], real dt)
{
  int d;
  for (d = 1; d < 3; d++) {
    p->xi[d] += dt * vxi[d];
  }
}

// ----------------------------------------------------------------------
// calc_vxi
//
// calculate velocity from momentum

__device__ static void
calc_vxi(real vxi[3], struct d_particle p)
{
  real root = rsqrtr(real(1.) + sqr(p.pxi[0]) + sqr(p.pxi[1]) + sqr(p.pxi[2]));

  int d;
  for (d = 0; d < 3; d++) {
    vxi[d] = p.pxi[d] * root;
  }
}

// ----------------------------------------------------------------------
// push_pxi_dt
//
// advance moments according to EM fields

__device__ static void
push_pxi_dt(struct d_particle *p,
	    real exq, real eyq, real ezq, real hxq, real hyq, real hzq)
{
  int kind = __float_as_int(p->kind_as_float);
  real dq = c_dqs[kind];
  real pxm = p->pxi[0] + dq*exq;
  real pym = p->pxi[1] + dq*eyq;
  real pzm = p->pxi[2] + dq*ezq;
  
  real root = dq * rsqrtr(real(1.) + sqr(pxm) + sqr(pym) + sqr(pzm));
  real taux = hxq * root, tauy = hyq * root, tauz = hzq * root;
  
  real tau = real(1.) / (real(1.) + sqr(taux) + sqr(tauy) + sqr(tauz));
  real pxp = ( (real(1.) + sqr(taux) - sqr(tauy) - sqr(tauz)) * pxm
	       +(real(2.)*taux*tauy + real(2.)*tauz)*pym
	       +(real(2.)*taux*tauz - real(2.)*tauy)*pzm)*tau;
  real pyp = ( (real(2.)*taux*tauy - real(2.)*tauz)*pxm
	       +(real(1.) - sqr(taux) + sqr(tauy) - sqr(tauz)) * pym
	       +(real(2.)*tauy*tauz + real(2.)*taux)*pzm)*tau;
  real pzp = ( (real(2.)*taux*tauz + real(2.)*tauy)*pxm
	       +(real(2.)*tauy*tauz - real(2.)*taux)*pym
	       +(real(1.) - sqr(taux) - sqr(tauy) + sqr(tauz))*pzm)*tau;
  
  p->pxi[0] = pxp + dq * exq;
  p->pxi[1] = pyp + dq * eyq;
  p->pxi[2] = pzp + dq * ezq;
}

// ----------------------------------------------------------------------
// push_part_one

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
push_part_one(struct d_particle *prt, int n, unsigned int *d_ids, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      real *fld_cache, int ci0[3], struct cuda_params prm)
{
  LOAD_PARTICLE_POS_(*prt, d_xi4, n);

  // here we have x^{n+.5}, p^n

  // field interpolation
  real exq, eyq, ezq, hxq, hyq, hzq;
  int lg[3];
  real og[3];
  find_idx_off_1st(prt->xi, lg, og, real(0.), prm);
  lg[1] -= ci0[1];
  lg[2] -= ci0[2];
  
  exq = ((1.f - og[1]) * (1.f - og[2]) * F3_CACHE(fld_cache, EX, lg[1]+0, lg[2]+0) +
	 (      og[1]) * (1.f - og[2]) * F3_CACHE(fld_cache, EX, lg[1]+1, lg[2]+0) +
	 (1.f - og[1]) * (      og[2]) * F3_CACHE(fld_cache, EX, lg[1]+0, lg[2]+1) +
	 (      og[1]) * (      og[2]) * F3_CACHE(fld_cache, EX, lg[1]+1, lg[2]+1));
  eyq = ((1.f - og[2]) * F3_CACHE(fld_cache, EY, lg[1]  , lg[2]+0) +
	 (      og[2]) * F3_CACHE(fld_cache, EY, lg[1]  , lg[2]+1));
  ezq = ((1.f - og[1]) * F3_CACHE(fld_cache, EZ, lg[1]+0, lg[2]  ) +
	 (      og[1]) * F3_CACHE(fld_cache, EZ, lg[1]+1, lg[2]  ));
  hxq = (F3_CACHE(fld_cache, HX, lg[1]  , lg[2]  ));
  hyq = ((1.f - og[1]) * F3_CACHE(fld_cache, HY, lg[1]+0, lg[2]  ) +
	 (      og[1]) * F3_CACHE(fld_cache, HY, lg[1]+1, lg[2]  ));
  hzq = ((1.f - og[2]) * F3_CACHE(fld_cache, HZ, lg[1]  , lg[2]+0) +
	 (      og[2]) * F3_CACHE(fld_cache, HZ, lg[1]  , lg[2]+1));

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  LOAD_PARTICLE_MOM_(*prt, d_pxi4, n);
  push_pxi_dt(prt, exq, eyq, ezq, hxq, hyq, hzq);
  STORE_PARTICLE_MOM_(*prt, d_pxi4, n);
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch(struct cuda_params prm, int *block_pos, int *ci0)
{
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % prm.b_mx[2];

  ci0[0] = 0;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.y / prm.b_mx[2];
}

__device__ static int
find_bid(struct cuda_params prm)
{
  return blockIdx.y * prm.b_mx[1] + blockIdx.x;
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
cache_fields(struct cuda_params prm, float *fld_cache, float *d_flds0, int size, int *ci0, int p)
{
  real *d_flds = d_flds0 + p * size;

  int ti = threadIdx.x;
  int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
  while (ti < n) {
    int tmp = ti;
    int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
    tmp /= BLOCKSIZE_Y + 4;
    int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
    // OPT? currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3_CACHE(fld_cache, m, jy, jz) = F3_DEV_YZ(m, jy+ci0[1],jz+ci0[2]);
    }
    ti += THREADS_PER_BLOCK;
  }
}

class GCurr {
public:
  real *d_flds;

  __device__ GCurr(real *_d_flds) :
    d_flds(_d_flds)
  {
  }

  __device__ void add(int m, int jy, int jz, float val, struct cuda_params prm, int *ci0)
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
	     GCurr &scurr, struct cuda_params prm, int *ci0)
{
  real xa[3] = { 0.,
		 x[1] + .5f * dx[1],
		 x[2] + .5f * dx[2], };
  if (dx[0] != 0.f) {
    real fnqx = qni_wni * prm.fnqxs;
    real h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
    scurr.add(0, i[1]  , i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), prm, ci0);
    scurr.add(0, i[1]+1, i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), prm, ci0);
    scurr.add(0, i[1]  , i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), prm, ci0);
    scurr.add(0, i[1]+1, i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), prm, ci0);
  }
  if (dx[1] != 0.f) {
    real fnqy = qni_wni * prm.fnqys;
    scurr.add(1, i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]), prm, ci0);
    scurr.add(1, i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]), prm, ci0);
  }
  if (dx[2] != 0.f) {
    real fnqz = qni_wni * prm.fnqzs;
    scurr.add(2, i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]), prm, ci0);
    scurr.add(2, i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]), prm, ci0);
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
// yz_calc_j

__device__ static void
yz_calc_j(struct d_particle *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	  GCurr &scurr,
	  struct cuda_params prm, int nr_total_blocks, int p_nr,
	  unsigned int *d_bidx, int bid, int *ci0)
{
  real vxi[3];
  calc_vxi(vxi, *prt);

  // position xm at x^(n+.5)
  real h0[3], h1[3];
  real xm[3], xp[3];
  int j[3], k[3];
  
  find_idx_off_pos_1st(prt->xi, j, h0, xm, real(0.), prm);

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

  // position xm at x^(n+.5)
  find_idx_off_pos_1st(prt->xi, k, h1, xp, real(0.), prm);

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
  curr_vb_cell(i, x, dx1, prt->qni_wni, scurr, prm, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
  
  off[1] = idiff[1] - off[1];
  off[2] = idiff[2] - off[2];
  calc_dx1(dx1, x, dx, off);
  curr_vb_cell(i, x, dx1, prt->qni_wni, scurr, prm, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
    
  curr_vb_cell(i, x, dx, prt->qni_wni, scurr, prm, ci0);
}

// ======================================================================

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      unsigned int *d_off, int nr_total_blocks, unsigned int *d_ids, unsigned int *d_bidx,
	      float *d_flds0, unsigned int size)
{
  int block_pos[3], ci0[3];
  int p, bid;
  p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, block_pos, ci0);
  bid = find_bid(prm);

  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, fld_cache, d_flds0, size, ci0, p);

  GCurr scurr(d_flds0 + p * size);
  
  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      (&prt, n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0, prm);

    yz_calc_j
      (&prt, n, d_xi4, d_pxi4, scurr, prm, nr_total_blocks, p, d_bidx, bid, ci0);
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
    fields_cuda_real_t *d_flds = mflds_sub->d_flds + p * size * mflds->nr_fields;
    check(cudaMemset(d_flds + JXI * size, 0, 3 * size * sizeof(*d_flds)));
  }
}

// ----------------------------------------------------------------------
// cuda_push_mprts_ab

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_ab(struct psc_mparticles *mprts, struct psc_mfields *mflds,
		   struct psc_mparticles *mprts_cuda)
{
  struct psc_mparticles_cuda *mprts_cuda_sub = psc_mparticles_cuda(mprts_cuda);
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  struct cuda_params prm;
  _set_params(&prm, ppsc, mflds, mprts_cuda);
  set_consts(&prm);

  unsigned int fld_size = mflds->nr_fields *
    mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

  zero_currents(mflds);

  int gx, gy;
  gx = prm.b_mx[1];
  gy = prm.b_mx[2] * mprts_cuda->nr_patches;

  dim3 dimGrid(gx, gy);

  push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (0, prm, mprts_cuda_sub->d_xi4, mprts_cuda_sub->d_pxi4,
     mprts_cuda_sub->d_alt_xi4, mprts_cuda_sub->d_alt_pxi4, mprts_cuda_sub->d_off,
     mprts_cuda_sub->nr_total_blocks, mprts_cuda_sub->d_ids, mprts_cuda_sub->d_bidx,
     mflds_sub->d_flds, fld_size);
  cuda_sync_if_enabled();

  _free_params(&prm);
}

// ----------------------------------------------------------------------
// yz_cuda_push_mprts

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
yz_cuda_push_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds,
		   struct psc_mparticles *mprts_cuda)
{
  struct psc_mparticles_cuda *mprts_cuda_sub = psc_mparticles_cuda(mprts_cuda);
    
  psc_mparticles_cuda_copy_to_dev(mprts_cuda);
  
  assert(!mprts_cuda_sub->need_reorder);
  cuda_push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(mprts, mflds, mprts_cuda);
}

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts_yz

void
cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds,
			  struct psc_mparticles *mprts_cuda)
{
  yz_cuda_push_mprts<1, 4, 4>(mprts, mflds, mprts_cuda);
}

