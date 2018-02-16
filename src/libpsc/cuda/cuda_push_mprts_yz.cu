
#include "cuda_iface.h"
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_mparticles_const.h"
#include "cuda_mfields_const.h"

#include "../psc_push_particles/inc_defs.h"

#include "psc.h" // FIXME

#include "interpolate.hxx"
#include "pushp.hxx"

#define BND (2) // FIXME

#define THREADS_PER_BLOCK (512)

using real_t = float;

enum DEPOSIT {
  DEPOSIT_VB_2D,
  DEPOSIT_VB_3D,
};

enum CURRMEM {
  CURRMEM_SHARED,
  CURRMEM_GLOBAL,
};

// FIXME
#define CUDA_BND_S_OOB (10)

// ----------------------------------------------------------------------

// OPT: use more shmem?

// OPT: passing shared memory cache etc around is probably sub-optimal
// OPT: fld cache is much bigger than needed
// OPT: precalculating IP coeffs could be a gain, too

// ======================================================================
// FldCache

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
struct FldCache
{
  using real_t = float;
  
  __device__ FldCache() = default;
  __device__ FldCache(const FldCache&) = delete;
  
  __device__ void load(DMFields d_flds0, int size, int *ci0, int p)
  {
    off_ = (-(ci0[2] - 2) * (BLOCKSIZE_Y + 4) +
	    -(ci0[1] - 2));
    DFields d_flds = d_flds0[p];
    
    int ti = threadIdx.x;
    int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
    while (ti < n) {
      int tmp = ti;
      int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
      tmp /= BLOCKSIZE_Y + 4;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      // OPT? currently it seems faster to do the loop rather than do m by threadidx
      for (int m = EX; m <= HZ; m++) {
	(*this)(m, 0, jy+ci0[1], jz+ci0[2]) = d_flds(m, 0,jy+ci0[1],jz+ci0[2]);
      }
      ti += THREADS_PER_BLOCK;
    }
  }

  __host__ __device__ float operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i,j,k)];
  }

private: // it's supposed to be a (read-only) cache, after all
  __host__ __device__ float& operator()(int m, int i, int j, int k)
  {
    return data_[index(m, i,j,k)];
  }

private:
  __host__ __device__ int index(int m, int i, int j, int k) const
  {
    return (((m - EX) * (BLOCKSIZE_Z + 4)
	     + k) * (BLOCKSIZE_Y + 4)
	    + j) + off_;
  }

  float data_[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  int off_;
};

// ----------------------------------------------------------------------
// push_xi
//
// advance position using velocity

__device__ static void
push_xi(struct d_particle *p, const float vxi[3], float dt)
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
calc_vxi(float vxi[3], struct d_particle p)
{
  float root = rsqrtf(float(1.) + sqr(p.pxi[0]) + sqr(p.pxi[1]) + sqr(p.pxi[2]));

  int d;
  for (d = 0; d < 3; d++) {
    vxi[d] = p.pxi[d] * root;
  }
}

// ----------------------------------------------------------------------
// push_part_one

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER,
	 typename OPT_IP, typename FldCache_t>
__device__ static void
push_part_one(struct d_particle& prt, int n, uint *d_ids, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      const FldCache_t& fld_cache, int ci0[3])
{
  uint id;
  if (REORDER) {
    id = d_ids[n];
    LOAD_PARTICLE_POS(prt, d_xi4, id);
  } else {
    LOAD_PARTICLE_POS(prt, d_xi4, n);
  }
  // here we have x^{n+.5}, p^n

  // field interpolation
  float xm[3];
  xm[1] = prt.xi[1] * d_cmprts_const.dxi[1];
  xm[2] = prt.xi[2] * d_cmprts_const.dxi[2];
  InterpolateEM<FldCache_t, OPT_IP, dim_yz> ip;
  ip.set_coeffs(xm);
  
  real_t E[3] = { ip.ex(fld_cache), ip.ey(fld_cache), ip.ez(fld_cache) };
  real_t H[3] = { ip.hx(fld_cache), ip.hy(fld_cache), ip.hz(fld_cache) };

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  int kind = __float_as_int(prt.kind_as_float);
  real_t dq = d_cmprts_const.dq[kind];
  if (REORDER) {
    LOAD_PARTICLE_MOM(prt, d_pxi4, id);
    push_p(prt.pxi, E, H, dq);
    STORE_PARTICLE_MOM(prt, d_alt_pxi4, n);
  } else {
    LOAD_PARTICLE_MOM(prt, d_pxi4, n);
    push_p(prt.pxi, E, H, dq);
    STORE_PARTICLE_MOM(prt, d_pxi4, n);
  }
}

// ----------------------------------------------------------------------
// SCurr

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

  float *scurr;
  DFields d_flds;

  __device__ SCurr(float *_scurr, DFields _d_flds) :
    scurr(_scurr), d_flds(_d_flds)
  {
    int i = threadIdx.x;
    while (i < shared_size) {
      scurr[i] = float(0.);
      i += THREADS_PER_BLOCK;
    }
  }

  __device__ void add_to_fld(int *ci0)
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
	float val = float(0.);
	// FIXME, OPT
	for (int wid = 0; wid < NR_CBLOCKS; wid++) {
	  val += (*this)(wid, jy, jz, m);
	}
	d_flds(JXI+m, 0,jy+ci0[1],jz+ci0[2]) += val;
	i += THREADS_PER_BLOCK;
      }
    }
  }

  __device__ float operator()(int wid, int jy, int jz, int m) const
  {
    uint off = CBLOCK_OFF(jy, jz, m, wid);
    return scurr[off];
  }
  __device__ float& operator()(int wid, int jy, int jz, int m)
  {
    uint off = CBLOCK_OFF(jy, jz, m, wid);
    return scurr[off];
  }

  __device__ void add(int m, int jy, int jz, float val, int *ci0)
  {
    float *addr = &(*this)(CBLOCK_ID, jy, jz, m);
    atomicAdd(addr, val);
  }
};

// ----------------------------------------------------------------------
// GCurr

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
class GCurr {
public:
  static const int shared_size = 1;

  float *scurr;
  DFields d_flds;

  __device__ GCurr(float *_scurr, DFields _d_flds) :
    scurr(_scurr), d_flds(_d_flds)
  {
  }

  __device__ void add_to_fld(int *ci0)
  {
  }

  __device__ void add(int m, int jy, int jz, float val, int *ci0)
  {
    float *addr = &d_flds(JXI+m, 0,jy+ci0[1],jz+ci0[2]);
    atomicAdd(addr, val);
  }
};

// ======================================================================
// depositing current

// ----------------------------------------------------------------------
// calc_dx1

__device__ static void
calc_dx1(float dx1[3], float x[3], float dx[3], int off[3])
{
  float o1, x1, dx_0, dx_1, dx_2, v0, v1, v2;
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

template<enum DEPOSIT DEPOSIT, class CURR>
__device__ static void
curr_vb_cell(int i[3], float x[3], float dx[3], float qni_wni,
	     CURR &scurr, int *ci0)
{
  float xa[3] = { 0.,
		 x[1] + .5f * dx[1],
		 x[2] + .5f * dx[2], };
  if (DEPOSIT == DEPOSIT_VB_3D) {
    if (dx[0] != 0.f) {
      float fnqx = qni_wni * d_cmprts_const.fnqxs;
      float h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
      scurr.add(0, i[1]  , i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), ci0);
      scurr.add(0, i[1]+1, i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), ci0);
      scurr.add(0, i[1]  , i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), ci0);
      scurr.add(0, i[1]+1, i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), ci0);
    }
  }
  if (dx[1] != 0.f) {
    float fnqy = qni_wni * d_cmprts_const.fnqys;
    scurr.add(1, i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]), ci0);
    scurr.add(1, i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]), ci0);
  }
  if (dx[2] != 0.f) {
    float fnqz = qni_wni * d_cmprts_const.fnqzs;
    scurr.add(2, i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]), ci0);
    scurr.add(2, i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]), ci0);
  }
}

// ----------------------------------------------------------------------
// curr_vb_cell_upd

__device__ static void
curr_vb_cell_upd(int i[3], float x[3], float dx1[3], float dx[3], int off[3])
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

template<enum DEPOSIT DEPOSIT, class CURR>
__device__ static void
yz_calc_j(struct d_particle *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	  CURR &scurr,
	  int nr_total_blocks, int p_nr,
	  uint *d_bidx, int bid, int *ci0)
{
  float vxi[3];
  calc_vxi(vxi, *prt);

  // position xm at x^(n+.5)
  float h0[3], h1[3];
  float xm[3], xp[3];
  int j[3], k[3];
  
  find_idx_off_pos_1st(prt->xi, j, h0, xm, float(0.));

  if (DEPOSIT == DEPOSIT_VB_2D) {
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    push_xi(prt, vxi, .5f * d_cmprts_const.dt);

    float fnqx = vxi[0] * prt->qni_wni * d_cmprts_const.fnqs;

    // out-of-plane currents at intermediate time
    int lf[3];
    float of[3];
    find_idx_off_1st(prt->xi, lf, of, float(0.));
    lf[1] -= ci0[1];
    lf[2] -= ci0[2];

    scurr.add(0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnqx, ci0);
    scurr.add(0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnqx, ci0);
    scurr.add(0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnqx, ci0);
    scurr.add(0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnqx, ci0);

    // x^(n+1.0), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(prt, vxi, .5f * d_cmprts_const.dt);
    STORE_PARTICLE_POS(*prt, d_xi4, n);
  } else if (DEPOSIT == DEPOSIT_VB_3D) {
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(prt, vxi, d_cmprts_const.dt);
    STORE_PARTICLE_POS(*prt, d_xi4, n);
  }

  // save block_idx for new particle position at x^(n+1.5)
  uint block_pos_y = __float2int_rd(prt->xi[1] * d_cmprts_const.b_dxi[1]);
  uint block_pos_z = __float2int_rd(prt->xi[2] * d_cmprts_const.b_dxi[2]);
  int nr_blocks = d_cmprts_const.b_mx[1] * d_cmprts_const.b_mx[2];

  int block_idx;
  if (block_pos_y >= d_cmprts_const.b_mx[1] || block_pos_z >= d_cmprts_const.b_mx[2]) {
    block_idx = CUDA_BND_S_OOB;
  } else {
    int bidx = block_pos_z * d_cmprts_const.b_mx[1] + block_pos_y + p_nr * nr_blocks;
    int b_diff = bid - bidx + d_cmprts_const.b_mx[1] + 1;
    int d1 = b_diff % d_cmprts_const.b_mx[1];
    int d2 = b_diff / d_cmprts_const.b_mx[1];
    block_idx = d2 * 3 + d1;
  }
  d_bidx[n] = block_idx;

  // position xm at x^(n+.5)
  find_idx_off_pos_1st(prt->xi, k, h1, xp, float(0.));

  // deposit xm -> xp
  int idiff[3] = { 0, k[1] - j[1], k[2] - j[2] };
  int i[3] = { 0, j[1] - ci0[1], j[2] - ci0[2] };
  float x[3] = { 0.f, xm[1] - j[1] - float(.5), xm[2] - j[2] - float(.5) };
  //float dx[3] = { 0.f, xp[1] - xm[1], xp[2] - xm[2] };
  float dx[3] = { vxi[0] * d_cmprts_const.dt * d_cmprts_const.dxi[0], xp[1] - xm[1], xp[2] - xm[2] };
  
  float x1 = x[1] * idiff[1];
  float x2 = x[2] * idiff[2];
  int d_first = (fabsf(dx[2]) * (.5f - x1) >= fabsf(dx[1]) * (.5f - x2));

  int off[3];
  if (d_first == 0) {
    off[1] = idiff[1];
    off[2] = 0;
  } else {
    off[1] = 0;
    off[2] = idiff[2];
  }

  float dx1[3];
  calc_dx1(dx1, x, dx, off);
  curr_vb_cell<DEPOSIT, CURR>(i, x, dx1, prt->qni_wni, scurr, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
  
  off[1] = idiff[1] - off[1];
  off[2] = idiff[2] - off[2];
  calc_dx1(dx1, x, dx, off);
  curr_vb_cell<DEPOSIT, CURR>(i, x, dx1, prt->qni_wni, scurr, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
    
  curr_vb_cell<DEPOSIT, CURR>(i, x, dx, prt->qni_wni, scurr, ci0);
}

// ======================================================================

#define DECLARE_AND_ZERO_SCURR						\
  __shared__ float _scurr[CURR::shared_size];				\
  CURR scurr(_scurr, d_flds0[p])					\

#define DECLARE_AND_CACHE_FIELDS					\

#define FIND_BLOCK_RANGE_CURRMEM(CURRMEM)				\
  int block_pos[3], ci0[3];						\
  int p, bid;								\
  if (CURRMEM == CURRMEM_SHARED) {					\
    p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
      (block_pos, ci0, block_start);					\
    if (p < 0)								\
      return;								\
    									\
    bid = find_bid_q(p, block_pos);					\
  } else if (CURRMEM == CURRMEM_GLOBAL) {				\
    p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
      (block_pos, ci0);							\
    bid = find_bid();							\
  }									\
  int block_begin = d_mprts.off_[bid];					\
  int block_end = d_mprts.off_[bid + 1]					\
    
// ----------------------------------------------------------------------
// push_mprts_ab

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER,
	 typename OPT_IP, enum DEPOSIT DEPOSIT, enum CURRMEM CURRMEM, class CURR>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(int block_start, DMParticles d_mprts, DMFields d_flds0, uint size)
{
  using FldCache_t = FldCache<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>;

  FIND_BLOCK_RANGE_CURRMEM(CURRMEM);
  __shared__ FldCache_t fld_cache;
  fld_cache.load(d_flds0, size, ci0, p);
  DECLARE_AND_ZERO_SCURR;
  
  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER, OPT_IP>
      (prt, n, d_mprts.id_, d_mprts.xi4_, d_mprts.pxi4_, d_mprts.alt_xi4_, d_mprts.alt_pxi4_, fld_cache, ci0);

    if (REORDER) {
      yz_calc_j<DEPOSIT_VB_2D, CURR>
	(&prt, n, d_mprts.alt_xi4_, d_mprts.alt_pxi4_, scurr, d_mprts.n_blocks_, p, d_mprts.bidx_, bid, ci0);
    } else {
      yz_calc_j<DEPOSIT_VB_2D, CURR>
	(&prt, n, d_mprts.xi4_, d_mprts.pxi4_, scurr, d_mprts.n_blocks_, p, d_mprts.bidx_, bid, ci0);
    }
  }
  
  scurr.add_to_fld(ci0);
}

// ----------------------------------------------------------------------
// zero_currents

static void
zero_currents(struct cuda_mfields *cmflds)
{
  // OPT: j as separate field, so we can use a single memset?
  for (int p = 0; p < cmflds->n_patches; p++) {
    uint size = cmflds->n_cells_per_patch;
    cudaError ierr = cudaMemset((*cmflds)[p].data() + JXI * size, 0,
				3 * size * sizeof(fields_cuda_real_t));
    cudaCheck(ierr);
  }
}

// ----------------------------------------------------------------------
// cuda_push_mprts_ab

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER,
	 typename OPT_IP, enum DEPOSIT DEPOSIT, enum CURRMEM CURRMEM>
static void
cuda_push_mprts_ab(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds)
{
  cuda_mfields_const_set(cmflds);
  cuda_mparticles_const_set(cmprts);

  uint fld_size = cmflds->n_fields * cmflds->n_cells_per_patch;

  zero_currents(cmflds);

  int gx, gy;
  if (CURRMEM == CURRMEM_SHARED) {
    gx =  (cmprts->b_mx_[1] + 1) / 2;
    gy = ((cmprts->b_mx_[2] + 1) / 2) * cmprts->n_patches;
  } else if (CURRMEM == CURRMEM_GLOBAL) {
    gx = cmprts->b_mx_[1];
    gy = cmprts->b_mx_[2] * cmprts->n_patches;
  }
  dim3 dimGrid(gx, gy);

  if (REORDER) {
    cmprts->d_alt_xi4.resize(cmprts->n_prts);
    cmprts->d_alt_pxi4.resize(cmprts->n_prts);
  }

  if (CURRMEM == CURRMEM_SHARED) {
    for (int block_start = 0; block_start < 4; block_start++) {
      push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER, OPT_IP, DEPOSIT, CURRMEM,
		    SCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> >
	<<<dimGrid, THREADS_PER_BLOCK>>>
	(block_start, *cmprts, *cmflds, fld_size);
      cuda_sync_if_enabled();
    }
  } else if (CURRMEM == CURRMEM_GLOBAL) {
    push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER, OPT_IP, DEPOSIT, CURRMEM,
    		  GCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> >
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (0, *cmprts, *cmflds, fld_size);
    cuda_sync_if_enabled();
  } else {
    assert(0);
  }

  if (REORDER) {
    cmprts->swap_alt();
    cmprts->need_reorder = false;
  }
}

// ----------------------------------------------------------------------
// yz_cuda_push_mprts

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, typename IP, enum DEPOSIT DEPOSIT,
	 enum CURRMEM CURRMEM>
static void
yz_cuda_push_mprts(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds)
{
  if (!cmprts->need_reorder) {
    //    printf("INFO: yz_cuda_push_mprts: need_reorder == false\n");
    cuda_push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false, IP, DEPOSIT, CURRMEM>(cmprts, cmflds);
  } else {
    cuda_push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, true, IP, DEPOSIT, CURRMEM>(cmprts, cmflds);
  }
}

// ----------------------------------------------------------------------
// cuda_push_mprts_yz

void
cuda_push_mprts_yz(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds,
		   const int bs[3], bool ip_ec, bool deposit_vb_3d, bool currmem_global)
{
  if (!ip_ec && !deposit_vb_3d && !currmem_global) {
    if (bs[0] == 1 && bs[1] == 4 && bs[2] == 4) {
      //return yz_cuda_push_mprts<1, 4, 4, opt_ip_1st, DEPOSIT_VB_2D, CURRMEM_SHARED>(cmprts, cmflds);
    }
  }
  
  if (ip_ec && deposit_vb_3d && !currmem_global) {
    if (bs[0] == 1 && bs[1] == 4 && bs[2] == 4) {
      return yz_cuda_push_mprts<1, 4, 4, opt_ip_1st_ec, DEPOSIT_VB_3D, CURRMEM_SHARED>(cmprts, cmflds);
    }
    if (bs[0] == 1 && bs[1] == 8 && bs[2] == 8) {
      //return yz_cuda_push_mprts<1, 8, 8, opt_ip_1st_ec, DEPOSIT_VB_3D, CURRMEM_SHARED>(cmprts, cmflds);
    }
  }

  if (ip_ec && deposit_vb_3d && currmem_global) {
    if (bs[0] == 1 && bs[1] == 4 && bs[2] == 4) {
      //return yz_cuda_push_mprts<1, 4, 4, opt_ip_1st_ec, DEPOSIT_VB_3D, CURRMEM_GLOBAL>(cmprts, cmflds);
    }
    if (bs[0] == 1 && bs[1] == 1 && bs[2] == 1) {
      return yz_cuda_push_mprts<1, 1, 1, opt_ip_1st_ec, DEPOSIT_VB_3D, CURRMEM_GLOBAL>(cmprts, cmflds);
    }
  }

  printf("bs %d:%d:%d ip_ec %d deposit_vb_3d %d currmem_global %d is not supported!\n",
	 bs[0], bs[1], bs[2], ip_ec, deposit_vb_3d, currmem_global);
  assert(0);
}

