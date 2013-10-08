
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <mrc_profile.h>

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


#define NO_CHECKERBOARD
//#define DEBUG

#define SW (2)

#include "cuda_common.h"

__device__ int *__d_error_count;

static __constant__ __device__ float c_dqs[4];

static void
set_consts(struct cuda_params *prm)
{
  check(cudaMemcpyToSymbol(c_dqs, prm->dq, sizeof(c_dqs)));
}

void
set_params(struct cuda_params *prm, struct psc *psc,
	   struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  prm->dt = psc->dt;
  for (int d = 0; d < 3; d++) {
    prm->dxi[d] = 1.f / ppsc->patch[0].dx[d];
  }

  prm->dqs    = .5f * psc->coeff.eta * psc->dt;
  prm->fnqs   = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
  prm->fnqys  = psc->patch[0].dx[1] * prm->fnqs / psc->dt;
  prm->fnqzs  = psc->patch[0].dx[2] * prm->fnqs / psc->dt;
  assert(psc->nr_kinds <= MAX_KINDS);
  for (int k = 0; k < psc->nr_kinds; k++) {
    prm->dq[k] = prm->dqs * psc->kinds[k].q / psc->kinds[k].m;
  }

  if (mprts && mprts->nr_patches > 0) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, 0);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);
    for (int d = 0; d < 3; d++) {
      prm->b_mx[d] = prts_cuda->b_mx[d];
      prm->b_dxi[d] = prts_cuda->b_dxi[d];
    }
  }

  if (mflds) {
    struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
    for (int d = 0; d < 3; d++) {
      prm->mx[d] = mflds_cuda->im[d];
      prm->ilg[d] = mflds_cuda->ib[d];
      if (d > 0) {
	assert(mflds_cuda->ib[d] == -BND);
      } else {
	assert(mflds_cuda->im[d] == 1);// + 2*BND);
      }
    }
  }

  //  check(cudaMalloc(&prm->d_error_count, 1 * sizeof(int)));
  //  check(cudaMemset(prm->d_error_count, 0, 1 * sizeof(int)));
}

void
free_params(struct cuda_params *prm)
{
  //  int h_error_count[1];
  //  check(cudaMemcpy(h_error_count, prm->d_error_count, 1 * sizeof(int),
  //		   cudaMemcpyDeviceToHost));
  //  check(cudaFree(prm->d_error_count));
  //  if (h_error_count[0] != 0) {
  //    printf("err cnt %d\n", h_error_count[0]);
  //  }
  //  assert(h_error_count[0] == 0);
}

// ======================================================================

void
psc_mparticles_cuda_copy_to_dev(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  check(cudaMemcpy(mprts_cuda->d_dev, mprts_cuda->h_dev,
		   mprts->nr_patches * sizeof(*mprts_cuda->d_dev),
		   cudaMemcpyHostToDevice));
}

// ======================================================================
// field caching

#if 0
template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
class F3cache {
  real *fld_cache;

public:
  __device__ F3cache(real *_fld_cache, real *d_flds, int l[3],
		     struct cuda_params prm) :
    fld_cache(_fld_cache)
  {
    int ti = threadIdx.x;
    int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
    while (ti < n) {
      int tmp = ti;
      int jx = tmp % BLOCKSIZE_X;
      tmp /= BLOCKSIZE_X;
      int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
      tmp /= BLOCKSIZE_Y + 4;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      //    tmp /= BLOCKSIZE_Z + 4;
      //    int m = tmp + EX;
      //    printf("n %d ti %d m %d, jx %d,%d,%d\n", n, ti, m, jx, jy, jz);
      // currently it seems faster to do the loop rather than do m by threadidx
      for (int m = EX; m <= HZ; m++) {
	(*this)(m, jx,jy,jz) = F3_DEV_YZ(m, jy+l[1],jz+l[2]);
      }
      ti += blockDim.x;
    }
    __syncthreads();
  }

  __host__ __device__ real operator()(int fldnr, int jx, int jy, int jz) const
  {
    int off = ((((fldnr-EX)
		 *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))
		*(BLOCKSIZE_Y + 4) + ((jy)-(-2)))
	       *1 + ((jx)));
    return fld_cache[off];
  }
  __host__ __device__ real& operator()(int fldnr, int jx, int jy, int jz)
  {
    int off = ((((fldnr-EX)
		 *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))
		*(BLOCKSIZE_Y + 4) + ((jy)-(-2)))
	       *1 + ((jx)));
    return fld_cache[off];
  }
};
#endif

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

#define OFF(g, d) o##g[d]
  
__device__ static real
ip1_to_grid_0(real h)
{
  return real(1.) - h;
}

__device__ static real
ip1_to_grid_p(real h)
{
  return h;
}

#define INTERP_FIELD_1ST(cache, exq, fldnr, g1, g2)			\
  do {									\
    int ddy = l##g1[1], ddz = l##g2[2];			\
    /* printf("C %g [%d,%d,%d]\n", F3C(fldnr, 0, ddy, ddz), 0, ddy, ddz); */ \
    exq =								\
      ip1_to_grid_0(OFF(g1, 1)) * ip1_to_grid_0(OFF(g2, 2)) *		\
      F3_CACHE(fld_cache, fldnr, ddy+0, ddz+0) +			\
      ip1_to_grid_p(OFF(g1, 1)) * ip1_to_grid_0(OFF(g2, 2)) *		\
      F3_CACHE(fld_cache, fldnr, ddy+1, ddz+0) +			\
      ip1_to_grid_0(OFF(g1, 1)) * ip1_to_grid_p(OFF(g2, 2)) *		\
      F3_CACHE(fld_cache, fldnr, ddy+0, ddz+1) +			\
      ip1_to_grid_p(OFF(g1, 1)) * ip1_to_grid_p(OFF(g2, 2)) *		\
      F3_CACHE(fld_cache, fldnr, ddy+1, ddz+1);				\
  } while(0)

// ----------------------------------------------------------------------
// push_part_one
//
// push one particle

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
push_part_one(struct d_particle *prt, int n, float4 *d_xi4, float4 *d_pxi4, real *fld_cache, int ci0[3],
	      struct cuda_params prm, bool do_read, bool do_write, bool do_push_pxi)
{
  // here we have x^{n+.5}, p^n
  
  // field interpolation

  int lh[3], lg[3];
  real oh[3], og[3];

  real exq, eyq, ezq, hxq, hyq, hzq;

  if (do_push_pxi) {
    find_idx_off_1st(prt->xi, lh, oh, real(-.5), prm);
    find_idx_off_1st(prt->xi, lg, og, real(0.), prm);
    lg[1] -= ci0[1];
    lh[1] -= ci0[1];
    lg[2] -= ci0[2];
    lh[2] -= ci0[2];

    INTERP_FIELD_1ST(cached_flds, exq, EX, g, g);
    INTERP_FIELD_1ST(cached_flds, eyq, EY, h, g);
    INTERP_FIELD_1ST(cached_flds, ezq, EZ, g, h);
    INTERP_FIELD_1ST(cached_flds, hxq, HX, h, h);
    INTERP_FIELD_1ST(cached_flds, hyq, HY, g, h);
    INTERP_FIELD_1ST(cached_flds, hzq, HZ, h, g);
  }

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  if (do_read) {
    LOAD_PARTICLE_MOM_(*prt, d_pxi4, n);
  }
  if (do_push_pxi) {
    push_pxi_dt(prt, exq, eyq, ezq, hxq, hyq, hzq);
  }
  if (do_write) {
    STORE_PARTICLE_MOM_(*prt, d_pxi4, n);
  }
}

// ----------------------------------------------------------------------
// push_part_one_reorder
//
// push one particle

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
push_part_one_reorder(struct d_particle *prt, int n, unsigned int *d_ids, float4 *d_xi4, float4 *d_pxi4,
		      float4 *d_alt_xi4, float4 *d_alt_pxi4,
		      real *fld_cache, int ci0[3], struct cuda_params prm, bool write_later,
		      bool do_read, bool do_write, bool do_push_pxi)
{
  unsigned int id = d_ids[n];
  if (do_read) {
    LOAD_PARTICLE_POS_(*prt, d_xi4, id);
  }
  if (!write_later) {
    if (do_write) {
      STORE_PARTICLE_POS_(*prt, d_alt_xi4, n);
    }
  }

  // here we have x^{n+.5}, p^n
  
  // field interpolation

  int lh[3], lg[3];
  real oh[3], og[3];
  real exq, eyq, ezq, hxq, hyq, hzq;

  if (do_push_pxi) {
    find_idx_off_1st(prt->xi, lh, oh, real(-.5), prm);
    find_idx_off_1st(prt->xi, lg, og, real(0.), prm);
    lg[1] -= ci0[1];
    lh[1] -= ci0[1];
    lg[2] -= ci0[2];
    lh[2] -= ci0[2];
    INTERP_FIELD_1ST(cached_flds, exq, EX, g, g);
    INTERP_FIELD_1ST(cached_flds, eyq, EY, h, g);
    INTERP_FIELD_1ST(cached_flds, ezq, EZ, g, h);
    INTERP_FIELD_1ST(cached_flds, hxq, HX, h, h);
    INTERP_FIELD_1ST(cached_flds, hyq, HY, g, h);
    INTERP_FIELD_1ST(cached_flds, hzq, HZ, h, g);
  }

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  if (do_read) {
    LOAD_PARTICLE_MOM_(*prt, d_pxi4, id);
  }
  if (do_push_pxi) {
    push_pxi_dt(prt, exq, eyq, ezq, hxq, hyq, hzq);
  }
  if (do_read) {
    STORE_PARTICLE_MOM_(*prt, d_alt_pxi4, n);
  }
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

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch_q(struct cuda_params prm, int *block_pos, int *ci0, int block_start)
{
  int grid_dim_y = (prm.b_mx[2] + 1) / 2;
  block_pos[1] = blockIdx.x * 2;
  block_pos[2] = (blockIdx.y % grid_dim_y) * 2;
  block_pos[1] += block_start & 1;
  block_pos[2] += block_start >> 1;
  if (block_pos[1] >= prm.b_mx[1] ||
      block_pos[2] >= prm.b_mx[2])
    return -1;

  ci0[0] = 0;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.y / grid_dim_y;
}

__device__ static int
find_bid(struct cuda_params prm)
{
  return blockIdx.y * prm.b_mx[1] + blockIdx.x;
}

__device__ static int
find_bid_q(struct cuda_params prm, int p, int *block_pos)
{
  // FIXME won't work if b_mx[1,2] not even (?)
  return block_pos_to_block_idx(block_pos, prm.b_mx) + p * prm.b_mx[1] * prm.b_mx[2];
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

// ----------------------------------------------------------------------
// push_mprts_p1
//
// push particles

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_p1(struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	      unsigned int *d_off, float *d_flds0, unsigned int size,
	      bool do_read, bool do_write, bool do_push_pxi)
{
  int block_pos[3], ci0[3];
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(prm, block_pos, ci0);
  int bid = find_bid(prm);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(prm, fld_cache, d_flds0, size, ci0, p);
  __syncthreads();

  float4 *xi4_begin = d_xi4 + block_begin;
  float4 *xi4 = d_xi4 + (block_begin & ~31) + threadIdx.x;
  float4 *pxi4 = d_pxi4 + (block_begin & ~31) + threadIdx.x;
  float4 *xi4_end = d_xi4 + block_end;

  for (; xi4 < xi4_end; xi4 += THREADS_PER_BLOCK, pxi4 += THREADS_PER_BLOCK) {
    if (xi4 >= xi4_begin) {
      struct d_particle prt;
      if (1||do_read) {
	LOAD_PARTICLE_POS_(prt, xi4, 0);
      }

      push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(&prt, 0, xi4, pxi4, fld_cache, ci0,
							   prm, 1||do_read, 1||do_write, 1||do_push_pxi);
    }
  }
}

// ----------------------------------------------------------------------
// push_mprts_p1q
//
// push particles, quarter at a time

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_p1q(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	       unsigned int *d_off, float *d_flds0, unsigned int size,
	       bool do_read, bool do_write, bool do_push_pxi)
{
  int block_pos[3], ci0[3];
  int p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, block_pos, ci0, block_start);
  if (p < 0)
    return;
  int bid = find_bid_q(prm, p, block_pos);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(prm, fld_cache, d_flds0, size, ci0, p);
  __syncthreads();

  float4 *xi4_begin = d_xi4 + block_begin;
  float4 *xi4 = d_xi4 + (block_begin & ~31) + threadIdx.x;
  float4 *pxi4 = d_pxi4 + (block_begin & ~31) + threadIdx.x;
  float4 *xi4_end = d_xi4 + block_end;

  for (; xi4 < xi4_end; xi4 += THREADS_PER_BLOCK, pxi4 += THREADS_PER_BLOCK) {
    if (xi4 >= xi4_begin) {
      struct d_particle prt;
      if (1||do_read) {
	LOAD_PARTICLE_POS_(prt, xi4, 0);
      }
      push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(&prt, 0, xi4, pxi4, fld_cache, ci0,
							   prm, 1||do_read, 1||do_write, 1||do_push_pxi);
    }
  }
}

// ----------------------------------------------------------------------
// push_mprts_p1_reorder
//
// push particles

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_p1_reorder(struct cuda_params prm, unsigned int *d_ids, float4 *d_xi4, float4 *d_pxi4,
		      float4 *d_alt_xi4, float4 *d_alt_pxi4,
		      unsigned int *d_off, float *d_flds0, unsigned int size)
{
  int block_pos[3], ci0[3];
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, block_pos, ci0);
  int bid = find_bid(prm);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(prm, fld_cache, d_flds0, size, ci0, p);
  __syncthreads();

  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      (&prt, n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0, prm, false,
       true, true, true);
  }
}

// ----------------------------------------------------------------------
// push_mprts_p1q_reorder
//
// push particles

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_p1q_reorder(int block_start, 
		       struct cuda_params prm, unsigned int *d_ids, float4 *d_xi4, float4 *d_pxi4,
		       float4 *d_alt_xi4, float4 *d_alt_pxi4,
		       unsigned int *d_off, float *d_flds0, unsigned int size)
{
  int block_pos[3], ci0[3];
  int p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, block_pos, ci0, block_start);
  if (p < 0)
    return;
  int bid = find_bid_q(prm, p, block_pos);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(prm, fld_cache, d_flds0, size, ci0, p);
  __syncthreads();

  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    push_part_one_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      (n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0, prm,
       true, true, true);
  }
}

// ----------------------------------------------------------------------
// cuda_push_mprts_a

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_a(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  struct cuda_params prm;
  set_params(&prm, ppsc, mprts, mflds);
  set_consts(&prm);

  unsigned int size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  
  dim3 dimGrid(prm.b_mx[1], prm.b_mx[2] * mprts->nr_patches);
  
  push_mprts_p1<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4, mprts_cuda->d_off,
     mflds_cuda->d_flds, size, prm.b_mx[1], prm.b_mx[2],
     true, true, true);
  cuda_sync_if_enabled();
  
  free_params(&prm);
}

// ----------------------------------------------------------------------
// cuda_push_mprts_aq

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_aq(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  struct cuda_params prm;
  set_params(&prm, ppsc, mprts, mflds);
  set_consts(&prm);

  unsigned int size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  
  dim3 dimGrid((prm.b_mx[1] + 1) / 2, ((prm.b_mx[2] + 1) / 2) * mprts->nr_patches);
  
  for (int block_start = 0; block_start < 4; block_start++) {
    push_mprts_p1q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4, mprts_cuda->d_off,
       mflds_cuda->d_flds, size,
       true, true, true);
    cuda_sync_if_enabled();
  }
  free_params(&prm);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_swap_alt
// FIXME, duplicated

static void
psc_mparticles_cuda_swap_alt(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  float4 *tmp_xi4 = mprts_cuda->d_alt_xi4;
  float4 *tmp_pxi4 = mprts_cuda->d_alt_pxi4;
  mprts_cuda->d_alt_xi4 = mprts_cuda->d_xi4;
  mprts_cuda->d_alt_pxi4 = mprts_cuda->d_pxi4;
  mprts_cuda->d_xi4 = tmp_xi4;
  mprts_cuda->d_pxi4 = tmp_pxi4;
}

// ----------------------------------------------------------------------
// cuda_push_mprts_a_reorder

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_a_reorder(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  struct cuda_params prm;
  set_params(&prm, ppsc, mprts, mflds);
  set_consts(&prm);

  psc_mparticles_cuda_copy_to_dev(mprts);

  unsigned int size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  
  dim3 dimGrid(prm.b_mx[1], prm.b_mx[2] * mprts->nr_patches);

  push_mprts_p1_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (prm, mprts_cuda->d_ids, mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
     mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4, mprts_cuda->d_off,
     mflds_cuda->d_flds, size);
  cuda_sync_if_enabled();
  
  psc_mparticles_cuda_swap_alt(mprts);
    
  free_params(&prm);
}

// OPT: take i < cell_end condition out of load
// OPT: reduce two at a time
// OPT: try splitting current calc / measuring by itself

// OPT: don't need as many ghost points for current and EM fields (?)

#define NR_CBLOCKS 16
#define CBLOCK_ID (threadIdx.x & (NR_CBLOCKS - 1))
#define CBLOCK_SIZE_Y (BLOCKSIZE_Y + 2*SW)
#define CBLOCK_SIZE_Z (BLOCKSIZE_Z + 2*SW)
#define CBLOCK_SIZE (CBLOCK_SIZE_Y * CBLOCK_SIZE_Z * (NR_CBLOCKS))

#define CBLOCK_OFF(jy, jz, wid) ((((jz) + SW) * CBLOCK_SIZE_Y + ((jy) + SW)) * (NR_CBLOCKS) + wid)

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
class SCurr {
  real *scurr;

public:
  __device__ SCurr(real *_scurr) :
    scurr(_scurr)
  {
  }

  __device__ void zero()
  {
    int i = threadIdx.x;
    while (i < CBLOCK_SIZE) {
      scurr[i] = real(0.);
      i += THREADS_PER_BLOCK;
    }
  }

  __device__ void add_to_fld(real *d_flds, int m, struct cuda_params prm, int *ci0)
  {
    int i = threadIdx.x;
    int stride = (BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW);
    while (i < stride) {
      int rem = i;
      int jz = rem / (BLOCKSIZE_Y + 2*SW);
      rem -= jz * (BLOCKSIZE_Y + 2*SW);
      int jy = rem;
      jz -= SW;
      jy -= SW;
      real val = real(0.);
      // FIXME, OPT
      for (int wid = 0; wid < NR_CBLOCKS; wid++) {
	val += (*this)(wid, jy, jz);
      }
      F3_DEV_YZ(JXI+m, jy+ci0[1],jz+ci0[2]) += val;
      i += THREADS_PER_BLOCK;
    }
  }

  __device__ real operator()(int wid, int jy, int jz) const
  {
    unsigned int off = CBLOCK_OFF(jy, jz, wid);
 #ifdef DEBUG
    if (off >= NR_CBLOCKS * BLOCK_STRIDE) {
      (*__d_error_count)++;
      off = 0;
    }
#endif

    return scurr[off];
  }
  __device__ real& operator()(int wid, int jy, int jz)
  {
    unsigned int off = CBLOCK_OFF(jy, jz, wid);
#ifdef DEBUG
    if (off >= NR_CBLOCKS * CBLOCK_STRIDE) {
      (*__d_error_count)++;
      off = 0;
    }
#endif

    return scurr[off];
  }

  __device__ real operator()(int jy, int jz) const
  {
    return (*this)(CBLOCK_ID, jy, jz);
  }

  __device__ real& operator()(int jy, int jz)
  {
    return (*this)(CBLOCK_ID, jy, jz);
  }
};

// ======================================================================

// ----------------------------------------------------------------------
// current_add

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
current_add(SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr, int jy, int jz, real val,
	    bool do_reduce)
{
  float *addr = &scurr(jy, jz);
  if (do_reduce) {
#if __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
    atomicAdd(addr, val);
#else
#if 0
    while ((val = atomicExch(addr, atomicExch(addr, 0.0f)+val))!=0.0f);
#else
    int lid = threadIdx.x & 31;
    for (int i = 0; i < 32; i++) {
      if (lid == i) {
	*addr += val;
      }
    }
#endif
#endif
  } else {
    *addr += val;
  }
}

// ----------------------------------------------------------------------
// yz_calc_jy

__device__ static void
calc_dx1(real dx1[2], real x[2], real dx[2], int off[2])
{
  real o0, x0, dx_0, dx_1, v0, v1;
  if (off[0] == 0) {
    o0 = off[1];
    x0 = x[1];
    dx_0 = dx[1];
    dx_1 = dx[0];
  } else {
    o0 = off[0];
    x0 = x[0];
    dx_0 = dx[0];
    dx_1 = dx[1];
  }
  if ((off[0] == 0 && off[1] == 0) || dx_0 == 0.f) {
    v0 = 0.f;
    v1 = 0.f;
  } else {
    v0 = .5f * o0 - x0;
    v1 = dx_1 / dx_0 * v0;
  }
  if (off[0] == 0) {
    dx1[0] = v1;
    dx1[1] = v0;
  } else {
    dx1[0] = v0;
    dx1[1] = v1;
  }
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
curr_2d_vb_cell(int i[2], real x[2], real dx[2], real qni_wni,
		SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr_y,
		SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr_z,
		struct cuda_params prm, bool do_reduce)
{
  if (dx[0] != 0.f) {
    real fnqy = qni_wni * prm.fnqys;
    current_add(scurr_y, i[0],i[1]  , fnqy * dx[0] * (.5f - x[1] - .5f * dx[1]), do_reduce);
    current_add(scurr_y, i[0],i[1]+1, fnqy * dx[0] * (.5f + x[1] + .5f * dx[1]), do_reduce);
  }
  if (dx[1] != 0.f) {
    real fnqz = qni_wni * prm.fnqzs;
    current_add(scurr_z, i[0],i[1]  , fnqz * dx[1] * (.5f - x[0] - .5f * dx[0]), do_reduce);
    current_add(scurr_z, i[0]+1,i[1], fnqz * dx[1] * (.5f + x[0] + .5f * dx[0]), do_reduce);
  }
}

__device__ static void
curr_2d_vb_cell_upd(int i[2], real x[2], real dx1[2], real dx[2], int off[2])
{
  dx[0] -= dx1[0];
  dx[1] -= dx1[1];
  x[0] += dx1[0] - off[0];
  x[1] += dx1[1] - off[1];
  i[0] += off[0];
  i[1] += off[1];
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
yz_calc_j(struct d_particle *prt, int i, float4 *d_xi4, float4 *d_pxi4,
	  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr_x,
	  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr_y,
	  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr_z,
	  struct cuda_params prm, int nr_total_blocks, int p_nr,
	  unsigned int *d_bidx, int bid, int *ci0,
	  bool do_read, bool do_write, bool do_reduce,
	  bool do_calc_jx, bool do_calc_jyjz)
{
  real vxi[3];
  calc_vxi(vxi, *prt);

  real h0[3], h1[3];
  real xm[3], xp[3];
  int j[3], k[3];
  int lf[3];
  real of[3];
  real fnqx;

  if (do_calc_jyjz) {
    find_idx_off_pos_1st(prt->xi, j, h0, xm, real(0.), prm);
  }

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
  push_xi(prt, vxi, .5f * prm.dt);
  
  if (do_calc_jx) {
    fnqx = vxi[0] * prt->qni_wni * prm.fnqs;
      
    find_idx_off_1st(prt->xi, lf, of, real(0.), prm);
    lf[1] -= ci0[1];
    lf[2] -= ci0[2];
  }

  // x^(n+1.0), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, .5f * prm.dt);
  if (do_write) {
    STORE_PARTICLE_POS_(*prt, d_xi4, i);
  }
  
  if (do_calc_jx) {
    current_add(scurr_x, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnqx, do_reduce);
    current_add(scurr_x, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnqx, do_reduce);
    current_add(scurr_x, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnqx, do_reduce);
    current_add(scurr_x, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnqx, do_reduce);
  }

  if (do_write) {
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
    d_bidx[i] = block_idx;
  }

  if (do_calc_jyjz) {
    find_idx_off_pos_1st(prt->xi, k, h1, xp, real(0.), prm);
    
    int idiff[2] = { k[1] - j[1], k[2] - j[2] };
    real dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };
    real x[2] = { xm[1] - j[1] - real(.5), xm[2] - j[2] - real(.5) };
    int i[2] = { j[1] - ci0[1], j[2] - ci0[2] };

    real x0 = x[0] * idiff[0];
    real x1 = x[1] * idiff[1];
    int d_first = (abs(dx[1]) * (.5f - x0) >= abs(dx[0]) * (.5f - x1));

    int off[2];
    if (d_first == 0) {
      off[0] = idiff[0];
      off[1] = 0;
    } else {
      off[0] = 0;
      off[1] = idiff[1];
    }
    real dx1[2];
    calc_dx1(dx1, x, dx, off);
    curr_2d_vb_cell(i, x, dx1, prt->qni_wni, scurr_y, scurr_z, prm, do_reduce);
    curr_2d_vb_cell_upd(i, x, dx1, dx, off);
    
    off[0] = idiff[0] - off[0];
    off[1] = idiff[1] - off[1];
    calc_dx1(dx1, x, dx, off);
    curr_2d_vb_cell(i, x, dx1, prt->qni_wni, scurr_y, scurr_z, prm, do_reduce);
    curr_2d_vb_cell_upd(i, x, dx1, dx, off);
    
    curr_2d_vb_cell(i, x, dx, prt->qni_wni, scurr_y, scurr_z, prm, do_reduce);
  }
}

// ======================================================================

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_p3(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	      unsigned int *d_off, int nr_total_blocks, unsigned int *d_bidx,
	      float *d_flds0, unsigned int size,
	      bool do_read, bool do_write, bool do_reduce, bool do_calc_jx,
	      bool do_calc_jyjz)
{
  __d_error_count = prm.d_error_count;

  __shared__ real _scurrx[CBLOCK_SIZE];
  __shared__ real _scurry[CBLOCK_SIZE];
  __shared__ real _scurrz[CBLOCK_SIZE];

  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_x(_scurrx);
  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_y(_scurry);
  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_z(_scurrz);

  if (do_write) {
    scurr_x.zero();
    scurr_y.zero();
    scurr_z.zero();
  }

  int block_pos[3], ci0[3];
  int p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, block_pos, ci0, block_start);
  if (p < 0)
    return;
  int bid = find_bid_q(prm, p, block_pos);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (1||do_read) {
      LOAD_PARTICLE_(prt, d_xi4, d_pxi4, n);
    }
    yz_calc_j(&prt, n, d_xi4, d_pxi4, scurr_x, scurr_y, scurr_z, prm, nr_total_blocks, p, d_bidx, bid, ci0,
	      1||do_read, 1||do_write, 1||do_reduce, 1||do_calc_jx, 1||do_calc_jyjz);
  }
  
  if (do_write) {
    __syncthreads();
    real *d_flds = d_flds0 + p * size;
    scurr_x.add_to_fld(d_flds, 0, prm, ci0);
    scurr_y.add_to_fld(d_flds, 1, prm, ci0);
    scurr_z.add_to_fld(d_flds, 2, prm, ci0);
  }
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_p13(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	       unsigned int *d_off, int nr_total_blocks, unsigned int *d_bidx,
	       float *d_flds0, unsigned int size,
	       bool do_read, bool do_write, bool do_reduce, bool do_push_pxi,
	       bool do_calc_jx, bool do_calc_jyjz)
{
  int block_pos[3], ci0[3];
  int p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, block_pos, ci0, block_start);
  if (p < 0)
    return;

  int bid = find_bid_q(prm, p, block_pos);
  int block_begin = d_off[bid];
  __shared__ int block_end;
  block_end = d_off[bid + 1];

  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(prm, fld_cache, d_flds0, size, ci0, p);

  __shared__ real _scurrx[CBLOCK_SIZE];
  __shared__ real _scurry[CBLOCK_SIZE];
  __shared__ real _scurrz[CBLOCK_SIZE];

  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_x(_scurrx);
  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_y(_scurry);
  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_z(_scurrz);

  if (do_write) {
    scurr_x.zero();
    scurr_y.zero();
    scurr_z.zero();
  }

  __syncthreads();

  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (1||do_read) {
      LOAD_PARTICLE_POS_(prt, d_xi4, n);
    }
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(&prt, n, d_xi4, d_pxi4, fld_cache, ci0,
							 prm, 1||do_read, 1||do_write, 1||do_push_pxi);
    yz_calc_j(&prt, n, d_xi4, d_pxi4, scurr_x, scurr_y, scurr_z, prm, nr_total_blocks, p, d_bidx, bid, ci0,
	      1||do_read, 1||do_write, 1||do_reduce, 1||do_calc_jx, 1||do_calc_jyjz);
  }
  
  if (do_write) {
    __syncthreads();
    real *d_flds = d_flds0 + p * size;
    scurr_x.add_to_fld(d_flds, 0, prm, ci0);
    scurr_y.add_to_fld(d_flds, 1, prm, ci0);
    scurr_z.add_to_fld(d_flds, 2, prm, ci0);
  }
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_p13_reorder(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
		       float4 *d_alt_xi4, float4 *d_alt_pxi4,
		       unsigned int *d_off, int nr_total_blocks, unsigned int *d_ids, unsigned int *d_bidx,
		       float *d_flds0, unsigned int size,
		       bool do_read, bool do_write, bool do_reduce, bool do_calc_jx,
		       bool do_calc_jyjz)
{
  int block_pos[3], ci0[3];
  int p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (prm, block_pos, ci0, block_start);
  if (p < 0)
    return;
  int bid = find_bid_q(prm, p, block_pos);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(prm, fld_cache, d_flds0, size, ci0, p);

  __shared__ real _scurrx[CBLOCK_SIZE];
  __shared__ real _scurry[CBLOCK_SIZE];
  __shared__ real _scurrz[CBLOCK_SIZE];

  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_x(_scurrx);
  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_y(_scurry);
  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> scurr_z(_scurrz);

  if (do_write) {
    scurr_x.zero();
    scurr_y.zero();
    scurr_z.zero();
  }
  __syncthreads();

  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      (&prt, n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0, prm, true,
       true, true, true);

    yz_calc_j(&prt, n, d_alt_xi4, d_alt_pxi4, scurr_x, scurr_y, scurr_z, prm, nr_total_blocks, p, d_bidx, bid, ci0,
	      1||do_read, 1||do_write, 1||do_reduce, 1||do_calc_jx, 1||do_calc_jyjz);
  }
  
  if (do_write) {
    __syncthreads();
    real *d_flds = d_flds0 + p * size;
    scurr_x.add_to_fld(d_flds, 0, prm, ci0);
    scurr_y.add_to_fld(d_flds, 1, prm, ci0);
    scurr_z.add_to_fld(d_flds, 2, prm, ci0);
  }
}


static void
zero_currents(struct psc_mfields *mflds)
{
  // FIXME, one memset should do OPT
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    check(cudaMemset(flds_cuda->d_flds + JXI * size, 0,
		     3 * size * sizeof(*flds_cuda->d_flds)));
  }
}

// ----------------------------------------------------------------------
// cuda_push_mprts_b

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_b(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  if (mprts->nr_patches == 0)
    return;

  struct cuda_params prm;
  set_params(&prm, ppsc, mprts, mflds);

  unsigned int fld_size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];

  zero_currents(mflds);
  
  dim3 dimGrid((prm.b_mx[1] + 1) / 2, ((prm.b_mx[2] + 1) / 2) * mprts->nr_patches);
  
  bool do_reduce = !(ppsc->timestep == -100);
  bool do_calc_jx = !(ppsc->timestep == -100);
  bool do_calc_jyjz = !(ppsc->timestep == -100);
  bool do_write = !(ppsc->timestep == -100);
  bool do_read = !(ppsc->timestep == -100);
  
  for (int block_start = 0; block_start < 4; block_start++) {
    push_mprts_p3<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4, mprts_cuda->d_off,
       mprts_cuda->nr_total_blocks, mprts_cuda->d_bidx,
       mflds_cuda->d_flds, fld_size,
       do_read, do_write, do_reduce, do_calc_jx, do_calc_jyjz);
    cuda_sync_if_enabled();
  }
  
  free_params(&prm);
}

// ----------------------------------------------------------------------
// cuda_push_mprts_ab

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_ab(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  struct cuda_params prm;
  set_params(&prm, ppsc, mprts, mflds);
  set_consts(&prm);

  unsigned int fld_size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];

  zero_currents(mflds);
  
  bool do_reduce = !(ppsc->timestep == -100);
  bool do_calc_jx = !(ppsc->timestep == -100);
  bool do_calc_jyjz = !(ppsc->timestep == -100);
  bool do_write = !(ppsc->timestep == -100);
  bool do_read = !(ppsc->timestep == -100);

  dim3 dimGrid((prm.b_mx[1] + 1) / 2, ((prm.b_mx[2] + 1) / 2) * mprts->nr_patches);
  
  for (int block_start = 0; block_start < 4; block_start++) {
    push_mprts_p13<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4, mprts_cuda->d_off,
       mprts_cuda->nr_total_blocks, mprts_cuda->d_bidx,
       mflds_cuda->d_flds, fld_size,
       do_read, do_write, do_reduce, true, do_calc_jx, do_calc_jyjz);
    cuda_sync_if_enabled();
  }

  free_params(&prm);
}

// ----------------------------------------------------------------------
// cuda_push_mprts_ab_reorder

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_ab_reorder(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  struct cuda_params prm;
  set_params(&prm, ppsc, mprts, mflds);
  set_consts(&prm);

  psc_mparticles_cuda_copy_to_dev(mprts);

  zero_currents(mflds);
  
  unsigned int fld_size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  
  bool do_reduce = !(ppsc->timestep == -100);
  bool do_calc_jx = !(ppsc->timestep == -100);
  bool do_calc_jyjz = !(ppsc->timestep == -100);
  bool do_write = !(ppsc->timestep == -100);
  bool do_read = !(ppsc->timestep == -100);
    
  dim3 dimGrid((prm.b_mx[1] + 1) / 2, ((prm.b_mx[2] + 1) / 2) * mprts->nr_patches);

  for (int block_start = 0; block_start < 4; block_start++) {
    // push_mprts_p1q_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    //   <<<dimGrid, THREADS_PER_BLOCK>>>
    //   (block_start, prm, mprts_cuda->d_ids, mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
    //    mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4, mprts_cuda->d_off,
    //    mflds_cuda->d_flds, fld_size);
    // cuda_sync_if_enabled();
    
    push_mprts_p13_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
       mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4, mprts_cuda->d_off,
       mprts_cuda->nr_total_blocks, mprts_cuda->d_ids, mprts_cuda->d_bidx,
       mflds_cuda->d_flds, fld_size,
       do_read, do_write, do_reduce, do_calc_jx, do_calc_jyjz);
    cuda_sync_if_enabled();
  }
    
  psc_mparticles_cuda_swap_alt(mprts);

  free_params(&prm);
}

// ======================================================================

EXTERN_C void
yz2x2_1vb_cuda_push_part_p2(struct psc_particles *prts, struct psc_fields *pf)
{
  assert(0);
  //  cuda_push_part_p2<1, 2, 2>(prts, pf);
}

EXTERN_C void
yz2x2_1vb_cuda_push_part_p3(struct psc_particles *prts, struct psc_fields *pf, real *dummy,
			    int block_stride)
{
  assert(0);
  //  cuda_push_part_p3<1, 2, 2>(prts, pf);
}

EXTERN_C void
yz8x8_1vb_cuda_push_part_p2(struct psc_particles *prts, struct psc_fields *pf)
{
  assert(0);
  //  cuda_push_part_p2<1, 8, 8>(prts, pf);
}

EXTERN_C void
yz8x8_1vb_cuda_push_part_p3(struct psc_particles *prts, struct psc_fields *pf, real *dummy,
			    int block_stride)
{
  assert(0);
  //  cuda_push_part_p3<1, 8, 8>(prts, pf);
}

static void
yz4x4_1vb_cuda_push_mprts_a(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_mprts_a", 1., 0, 0);
  }

  if (mprts->nr_patches == 0) {
    return;
  }

  prof_start(pr);
  psc_mparticles_cuda_copy_to_dev(mprts);
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (!mprts_cuda->need_reorder) {
    MHERE;
    cuda_push_mprts_aq<1, 4, 4>(mprts, mflds);
  } else {
    cuda_push_mprts_a_reorder<1, 4, 4>(mprts, mflds);
    mprts_cuda->need_reorder = false;
  }
  prof_stop(pr);
}

static void
yz4x4_1vb_cuda_push_mprts_b(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_mprts_b", 1., 0, 0);
  }

  prof_start(pr);
  cuda_push_mprts_b<1, 4, 4>(mprts, mflds);
  prof_stop(pr);
}

EXTERN_C void
yz4x4_1vb_cuda_push_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_mprts", 1., 0, 0);
  }

  prof_start(pr);
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
    
  if (0) {
    yz4x4_1vb_cuda_push_mprts_a(mprts, mflds);
    yz4x4_1vb_cuda_push_mprts_b(mprts, mflds);
  } else {
    psc_mparticles_cuda_copy_to_dev(mprts);
    
    if (!mprts_cuda->need_reorder) {
      MHERE;
      cuda_push_mprts_ab<1, 4, 4>(mprts, mflds);
    } else {
      cuda_push_mprts_ab_reorder<1, 4, 4>(mprts, mflds);
      mprts_cuda->need_reorder = false;
    }
  }

  prof_stop(pr);
}
