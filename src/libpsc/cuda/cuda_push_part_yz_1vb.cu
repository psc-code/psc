
#include "psc_cuda.h"
#include "particles_cuda.h"

#undef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK (512)

enum IP {
  IP_STD, // standard interpolation
  IP_EC,  // energy-conserving interpolation
};

enum DEPOSIT {
  DEPOSIT_VB_2D,
  DEPOSIT_VB_3D,
};

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

#include "cuda_common.h"

static __constant__ __device__ float c_dqs[4]; // FIXME hardcoded

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
  prm->fnqxs  = psc->patch[0].dx[0] * prm->fnqs / psc->dt;
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
}

void
free_params(struct cuda_params *prm)
{
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

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER, enum IP IP>
__device__ static void
push_part_one(struct d_particle *prt, int n, unsigned int *d_ids, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      real *fld_cache, int ci0[3], struct cuda_params prm)
{
  unsigned int id;
  if (REORDER) {
    id = d_ids[n];
    LOAD_PARTICLE_POS_(*prt, d_xi4, id);
  } else {
    LOAD_PARTICLE_POS_(*prt, d_xi4, n);
  }
  // here we have x^{n+.5}, p^n

  // field interpolation
  real exq, eyq, ezq, hxq, hyq, hzq;
  int lg[3];
  real og[3];
  find_idx_off_1st(prt->xi, lg, og, real(0.), prm);
  lg[1] -= ci0[1];
  lg[2] -= ci0[2];
  
  if (IP == IP_STD) {
    int lh[3];
    real oh[3];
    
    find_idx_off_1st(prt->xi, lh, oh, real(-.5), prm);
    lh[1] -= ci0[1];
    lh[2] -= ci0[2];
    INTERP_FIELD_1ST(cached_flds, exq, EX, g, g);
    INTERP_FIELD_1ST(cached_flds, eyq, EY, h, g);
    INTERP_FIELD_1ST(cached_flds, ezq, EZ, g, h);
    INTERP_FIELD_1ST(cached_flds, hxq, HX, h, h);
    INTERP_FIELD_1ST(cached_flds, hyq, HY, g, h);
    INTERP_FIELD_1ST(cached_flds, hzq, HZ, h, g);
  } else if (IP == IP_EC) {
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
  } else {
    assert(0);
  }

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  if (REORDER) {
    LOAD_PARTICLE_MOM_(*prt, d_pxi4, id);
    push_pxi_dt(prt, exq, eyq, ezq, hxq, hyq, hzq);
    STORE_PARTICLE_MOM_(*prt, d_alt_pxi4, n);
  } else {
    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    LOAD_PARTICLE_MOM_(*prt, d_pxi4, n);
    push_pxi_dt(prt, exq, eyq, ezq, hxq, hyq, hzq);
    STORE_PARTICLE_MOM_(*prt, d_pxi4, n);
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

// OPT: take i < cell_end condition out of load
// OPT: reduce two at a time
// OPT: try splitting current calc / measuring by itself

// OPT: don't need as many ghost points for current and EM fields (?)

#define BND_CURR_L (1)
#define BND_CURR_R (2)

#define NR_CBLOCKS 16
#define CBLOCK_ID (threadIdx.x & (NR_CBLOCKS - 1))
#define CBLOCK_SIZE_Y (BLOCKSIZE_Y + BND_CURR_L + BND_CURR_R)
#define CBLOCK_SIZE_Z (BLOCKSIZE_Z + BND_CURR_L + BND_CURR_R)
#define CBLOCK_SIZE (CBLOCK_SIZE_Y * CBLOCK_SIZE_Z * (NR_CBLOCKS))

#define CBLOCK_OFF(jy, jz, m, wid) ((((m) * CBLOCK_SIZE_Z + ((jz) + BND_CURR_L)) * CBLOCK_SIZE_Y + ((jy) + BND_CURR_L)) * (NR_CBLOCKS) + wid)

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
class SCurr {
public:
  static const int shared_size = 3 * CBLOCK_SIZE;

  real *scurr;
  real *d_flds;

  __device__ SCurr(real *_scurr, real *_d_flds) :
    scurr(_scurr), d_flds(_d_flds)
  {
    int i = threadIdx.x;
    while (i < shared_size) {
      scurr[i] = real(0.);
      i += THREADS_PER_BLOCK;
    }
  }

  __device__ void add_to_fld(struct cuda_params prm, int *ci0)
  {
    __syncthreads();				\
    int i = threadIdx.x;
    int stride = (BLOCKSIZE_Y + BND_CURR_L + BND_CURR_R) * (BLOCKSIZE_Z + BND_CURR_L + BND_CURR_R);
    while (i < stride) {
      int rem = i;
      int jz = rem / (BLOCKSIZE_Y + BND_CURR_L + BND_CURR_R);
      rem -= jz * (BLOCKSIZE_Y + BND_CURR_L + BND_CURR_R);
      int jy = rem;
      jz -= BND_CURR_L;
      jy -= BND_CURR_L;
      for (int m = 0; m < 3; m++) {
	real val = real(0.);
	// FIXME, OPT
	for (int wid = 0; wid < NR_CBLOCKS; wid++) {
	  val += (*this)(wid, jy, jz, m);
	}
	F3_DEV_YZ(JXI+m, jy+ci0[1],jz+ci0[2]) += val;
	i += THREADS_PER_BLOCK;
      }
    }
  }

  __device__ real operator()(int wid, int jy, int jz, int m) const
  {
    unsigned int off = CBLOCK_OFF(jy, jz, m, wid);
    return scurr[off];
  }
  __device__ real& operator()(int wid, int jy, int jz, int m)
  {
    unsigned int off = CBLOCK_OFF(jy, jz, m, wid);
    return scurr[off];
  }

  __device__ void add(int m, int jy, int jz, float val, struct cuda_params prm, int *ci0)
  {
    float *addr = &(*this)(CBLOCK_ID, jy, jz, m);
    atomicAdd(addr, val);
#if 0
    float *d_flds = scurr.d_flds;
    float *addr = &F3_DEV_YZ(JXI+m, jy+ci0[1],jz+ci0[2]);
    atomicAdd(addr, val);
#endif
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

template<enum DEPOSIT DEPOSIT, int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
curr_vb_cell(int i[3], real x[3], real dx[3], real qni_wni,
	     SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr,
	     struct cuda_params prm, int *ci0)
{
  real xa[3] = { 0.,
		 x[1] + .5f * dx[1],
		 x[2] + .5f * dx[2], };
  if (DEPOSIT == DEPOSIT_VB_3D) {
    if (dx[0] != 0.f) {
      real fnqx = qni_wni * prm.fnqxs;
      real h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
      scurr.add(0, i[1]  , i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), prm, ci0);
      scurr.add(0, i[1]+1, i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), prm, ci0);
      scurr.add(0, i[1]  , i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), prm, ci0);
      scurr.add(0, i[1]+1, i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), prm, ci0);
    }
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

template<enum DEPOSIT DEPOSIT, int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
yz_calc_j(struct d_particle *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> &scurr,
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

  if (DEPOSIT == DEPOSIT_VB_2D) {
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    push_xi(prt, vxi, .5f * prm.dt);

    real fnqx = vxi[0] * prt->qni_wni * prm.fnqs;

    // out-of-plane currents at intermediate time
    int lf[3];
    real of[3];
    find_idx_off_1st(prt->xi, lf, of, real(0.), prm);
    lf[1] -= ci0[1];
    lf[2] -= ci0[2];

    scurr.add(0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnqx, prm, ci0);
    scurr.add(0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnqx, prm, ci0);
    scurr.add(0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnqx, prm, ci0);
    scurr.add(0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnqx, prm, ci0);

    // x^(n+1.0), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(prt, vxi, .5f * prm.dt);
    STORE_PARTICLE_POS_(*prt, d_xi4, n);
  } else if (DEPOSIT == DEPOSIT_VB_3D) {
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(prt, vxi, prm.dt);
    STORE_PARTICLE_POS_(*prt, d_xi4, n);
  }

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
  curr_vb_cell<DEPOSIT>(i, x, dx1, prt->qni_wni, scurr, prm, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
  
  off[1] = idiff[1] - off[1];
  off[2] = idiff[2] - off[2];
  calc_dx1(dx1, x, dx, off);
  curr_vb_cell<DEPOSIT>(i, x, dx1, prt->qni_wni, scurr, prm, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
    
  curr_vb_cell<DEPOSIT>(i, x, dx, prt->qni_wni, scurr, prm, ci0);
}

// ======================================================================

#define DECLARE_AND_ZERO_SCURR						\
  __shared__ real _scurr[SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>::shared_size]; \
  SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>				\
  scurr(_scurr, d_flds0 + p * size)					\

#define DECLARE_AND_CACHE_FIELDS					\
  __shared__ real fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)]; \
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>			\
  (prm, fld_cache, d_flds0, size, ci0, p)

#define FIND_BLOCK_RANGE						\
  int block_pos[3], ci0[3];						\
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
    (prm, block_pos, ci0);						\
  int bid = find_bid(prm);						\
  int block_begin = d_off[bid];						\
  int block_end = d_off[bid + 1]

#define FIND_BLOCK_RANGE_Q						\
  int block_pos[3], ci0[3];						\
  int p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
    (prm, block_pos, ci0, block_start);					\
  if (p < 0)								\
    return;								\
									\
  int bid = find_bid_q(prm, p, block_pos);				\
  int block_begin = d_off[bid];						\
  int block_end = d_off[bid + 1]					\

#define FIND_BLOCK_RANGE_QS						\
  int block_pos[3], ci0[3];						\
  int p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
    (prm, block_pos, ci0, block_start);					\
  if (p < 0)								\
    return;								\
									\
  int bid = find_bid_q(prm, p, block_pos);				\
  int block_begin = d_off[bid];						\
  __shared__ int block_end;						\
  block_end = d_off[bid + 1];						\


// ----------------------------------------------------------------------
// push_mprts_a

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_a(struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	     unsigned int *d_off, float *d_flds0, unsigned int size)
{
  FIND_BLOCK_RANGE;
  DECLARE_AND_CACHE_FIELDS;

  __syncthreads();
  float4 *xi4_begin = d_xi4 + block_begin;
  float4 *xi4 = d_xi4 + (block_begin & ~31) + threadIdx.x;
  float4 *pxi4 = d_pxi4 + (block_begin & ~31) + threadIdx.x;
  float4 *xi4_end = d_xi4 + block_end;

  for (; xi4 < xi4_end; xi4 += THREADS_PER_BLOCK, pxi4 += THREADS_PER_BLOCK) {
    if (xi4 >= xi4_begin) {
      struct d_particle prt;
      LOAD_PARTICLE_POS_(prt, xi4, 0);
      push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false, 0>
	(&prt, 0, NULL, d_xi4, d_pxi4, NULL, NULL, fld_cache, ci0, prm);
    }
  }
}

// ----------------------------------------------------------------------
// push_mprts_aq
//
// push particles, quarter at a time

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_aq(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	      unsigned int *d_off, float *d_flds0, unsigned int size)
{
  FIND_BLOCK_RANGE_Q;
  DECLARE_AND_CACHE_FIELDS;

  __syncthreads();
  float4 *xi4_begin = d_xi4 + block_begin;
  float4 *xi4 = d_xi4 + (block_begin & ~31) + threadIdx.x;
  float4 *pxi4 = d_pxi4 + (block_begin & ~31) + threadIdx.x;
  float4 *xi4_end = d_xi4 + block_end;

  for (; xi4 < xi4_end; xi4 += THREADS_PER_BLOCK, pxi4 += THREADS_PER_BLOCK) {
    if (xi4 >= xi4_begin) {
      struct d_particle prt;
      LOAD_PARTICLE_POS_(prt, xi4, 0);
      push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false, IP_STD>
	(&prt, 0, NULL, d_xi4, d_pxi4, NULL, NULL, fld_cache, ci0, prm);
    }
  }
}

// ----------------------------------------------------------------------
// push_mprts_a_reorder
//
// push particles

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_a_reorder(struct cuda_params prm, unsigned int *d_ids, float4 *d_xi4, float4 *d_pxi4,
		     float4 *d_alt_xi4, float4 *d_alt_pxi4,
		     unsigned int *d_off, float *d_flds0, unsigned int size)
{
  FIND_BLOCK_RANGE;
  DECLARE_AND_CACHE_FIELDS;

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, true, IP_STD>
      (&prt, n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0, prm);
    STORE_PARTICLE_POS_(prt, d_alt_xi4, 0);
  }
}

// ----------------------------------------------------------------------
// push_mprts_aq_reorder

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_aq_reorder(int block_start, 
		      struct cuda_params prm, unsigned int *d_ids, float4 *d_xi4, float4 *d_pxi4,
		      float4 *d_alt_xi4, float4 *d_alt_pxi4,
		      unsigned int *d_off, float *d_flds0, unsigned int size)
{
  FIND_BLOCK_RANGE_Q;
  DECLARE_AND_CACHE_FIELDS;

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, true, 0>
      (&prt, n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0, prm);
    STORE_PARTICLE_POS_(prt, d_alt_xi4, 0);
  }
}

// ----------------------------------------------------------------------
// push_mprts_b

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_b(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	     unsigned int *d_off, int nr_total_blocks, unsigned int *d_bidx,
	     float *d_flds0, unsigned int size)
{
  FIND_BLOCK_RANGE_Q;
  DECLARE_AND_ZERO_SCURR;

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    LOAD_PARTICLE_(prt, d_xi4, d_pxi4, n);
    yz_calc_j<DEPOSIT_VB_2D>(&prt, n, d_xi4, d_pxi4, scurr, prm, nr_total_blocks, p, d_bidx, bid, ci0);
  }
  
  scurr.add_to_fld(prm, ci0);
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER, enum IP IP, enum DEPOSIT DEPOSIT>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      unsigned int *d_off, int nr_total_blocks, unsigned int *d_ids, unsigned int *d_bidx,
	      float *d_flds0, unsigned int size)
{
  FIND_BLOCK_RANGE_Q;
  DECLARE_AND_CACHE_FIELDS;
  DECLARE_AND_ZERO_SCURR;
  
  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER, IP>
      (&prt, n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0, prm);

    if (REORDER) {
      yz_calc_j<DEPOSIT>(&prt, n, d_alt_xi4, d_alt_pxi4, scurr, prm, 
			 nr_total_blocks, p, d_bidx, bid, ci0);
    } else {
      yz_calc_j<DEPOSIT>(&prt, n, d_xi4, d_pxi4, scurr, prm, 
			 nr_total_blocks, p, d_bidx, bid, ci0);
    }
  }
  
  scurr.add_to_fld(prm, ci0);
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
// zero_currents

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

#define CUDA_PUSH_MPRTS_TOP						\
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);	\
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);	\
									\
  struct cuda_params prm;						\
  set_params(&prm, ppsc, mprts, mflds);					\
  set_consts(&prm);							\
									\
  unsigned int fld_size = mflds->nr_fields *				\
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];		\
									\
  zero_currents(mflds);							\
									\
  dim3 dimGrid((prm.b_mx[1] + 1) / 2, ((prm.b_mx[2] + 1) / 2) * mprts->nr_patches) \


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
  
  push_mprts_a<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4, mprts_cuda->d_off,
     mflds_cuda->d_flds, size, prm.b_mx[1], prm.b_mx[2]);
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
    push_mprts_aq<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4, mprts_cuda->d_off,
       mflds_cuda->d_flds, size);
    cuda_sync_if_enabled();
  }
  free_params(&prm);
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

  push_mprts_a_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (prm, mprts_cuda->d_ids, mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
     mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4, mprts_cuda->d_off,
     mflds_cuda->d_flds, size);
  cuda_sync_if_enabled();
  
  psc_mparticles_cuda_swap_alt(mprts);
    
  free_params(&prm);
}

// ----------------------------------------------------------------------
// cuda_push_mprts_b

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
cuda_push_mprts_b(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  CUDA_PUSH_MPRTS_TOP;

  for (int block_start = 0; block_start < 4; block_start++) {
    push_mprts_b<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4, mprts_cuda->d_off,
       mprts_cuda->nr_total_blocks, mprts_cuda->d_bidx,
       mflds_cuda->d_flds, fld_size);
    cuda_sync_if_enabled();
  }
  
  free_params(&prm);
}

// ----------------------------------------------------------------------
// cuda_push_mprts_ab

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER, enum IP IP, enum DEPOSIT DEPOSIT>
static void
cuda_push_mprts_ab(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  CUDA_PUSH_MPRTS_TOP;

  for (int block_start = 0; block_start < 4; block_start++) {
    push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER, IP, DEPOSIT>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
       mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4, mprts_cuda->d_off,
       mprts_cuda->nr_total_blocks, mprts_cuda->d_ids, mprts_cuda->d_bidx,
       mflds_cuda->d_flds, fld_size);
    cuda_sync_if_enabled();
  }

  if (REORDER) {
    psc_mparticles_cuda_swap_alt(mprts);
  }

  free_params(&prm);
}

// ----------------------------------------------------------------------
// yz4x4_1vb_cuda_push_mprts_separate
//
// superseded by the combined pusher

EXTERN_C void
yz4x4_1vb_cuda_push_mprts_separate(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  psc_mparticles_cuda_copy_to_dev(mprts);

  if (!mprts_cuda->need_reorder) {
    MHERE;
    cuda_push_mprts_aq<1, 4, 4>(mprts, mflds);
  } else {
    cuda_push_mprts_a_reorder<1, 4, 4>(mprts, mflds);
    mprts_cuda->need_reorder = false;
  }
  cuda_push_mprts_b<1, 4, 4>(mprts, mflds);
}

// ----------------------------------------------------------------------
// yz_cuda_push_mprts

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, enum IP IP, enum DEPOSIT DEPOSIT>
static void
yz_cuda_push_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
    
  psc_mparticles_cuda_copy_to_dev(mprts);
  
  if (!mprts_cuda->need_reorder) {
    MHERE;
    cuda_push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false, IP, DEPOSIT>(mprts, mflds);
  } else {
    cuda_push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, true, IP, DEPOSIT>(mprts, mflds);
    mprts_cuda->need_reorder = false;
  }
}

// ----------------------------------------------------------------------
// yz4x4_1vb_cuda_push_mprts

EXTERN_C void
yz4x4_1vb_cuda_push_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  yz_cuda_push_mprts<1, 4, 4, IP_STD, DEPOSIT_VB_2D>(mprts, mflds);
}

// ----------------------------------------------------------------------
// yz4x4_1vbec3d_cuda_push_mprts

EXTERN_C void
yz4x4_1vbec3d_cuda_push_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  yz_cuda_push_mprts<1, 4, 4, IP_EC, DEPOSIT_VB_3D>(mprts, mflds);
}
