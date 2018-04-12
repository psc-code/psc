
#include "cuda_iface.h"
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_push_particles.cuh"
#include "push_particles_cuda_impl.hxx"
#include "range.hxx"

#define DIM DIM_YZ

#include "../psc_push_particles/inc_defs.h"

#include "psc.h" // FIXME

#include "dim.hxx"

using dim = dim_yz;

#include "interpolate.hxx"
#include "pushp.hxx"

#define BND (2) // FIXME

#define THREADS_PER_BLOCK (512)

using real_t = float;

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
  
  __device__ void load(DFields d_flds, int *ci0)
  {
    off_ = (-(ci0[2] - 2) * (BLOCKSIZE_Y + 4) +
	    -(ci0[1] - 2));
    
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
// push_part_one

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER,
	 typename OPT_IP, typename FldCache_t>
__device__ static void
push_part_one(DMparticlesCuda& dmprts, struct d_particle& prt, int n, const FldCache_t& fld_cache)
{
  uint id;
  if (REORDER) {
    id = dmprts.id_[n];
    LOAD_PARTICLE_POS(prt, dmprts.xi4_, id);
  } else {
    LOAD_PARTICLE_POS(prt, dmprts.xi4_, n);
  }
  // here we have x^{n+.5}, p^n

  // field interpolation
  real_t xm[3];
  dmprts.scalePos(xm, prt.xi);
  InterpolateEM<FldCache_t, OPT_IP, dim_yz> ip;
  AdvanceParticle<real_t, dim> advance{dmprts.dt()};

  ip.set_coeffs(xm);
  
  real_t E[3] = { ip.ex(fld_cache), ip.ey(fld_cache), ip.ez(fld_cache) };
  real_t H[3] = { ip.hx(fld_cache), ip.hy(fld_cache), ip.hz(fld_cache) };

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  int kind = __float_as_int(prt.kind_as_float);
  real_t dq = dmprts.dq(kind);
  if (REORDER) {
    LOAD_PARTICLE_MOM(prt, dmprts.pxi4_, id);
    advance.push_p(prt.pxi, E, H, dq);
    STORE_PARTICLE_MOM(prt, dmprts.alt_pxi4_, n);
  } else {
    LOAD_PARTICLE_MOM(prt, dmprts.pxi4_, n);
    advance.push_p(prt.pxi, E, H, dq);
    STORE_PARTICLE_MOM(prt, dmprts.pxi4_, n);
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
#define CBLOCK_SIZE_Y (BS_Y + BND_CURR_L + BND_CURR_R)
#define CBLOCK_SIZE_Z (BS_Z + BND_CURR_L + BND_CURR_R)
#define CBLOCK_SIZE (CBLOCK_SIZE_Y * CBLOCK_SIZE_Z * (NR_CBLOCKS))

#define CBLOCK_OFF(jy, jz, m, wid) ((((m) * CBLOCK_SIZE_Z + ((jz) + BND_CURR_L)) * CBLOCK_SIZE_Y + ((jy) + BND_CURR_L)) * (NR_CBLOCKS) + wid)

template<typename BS>
class SCurr
{
  static const int BS_X = BS::x::value, BS_Y = BS::y::value, BS_Z = BS::z::value;

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
    int stride = (BS_Y + BND_CURR_L + BND_CURR_R) * (BS_Z + BND_CURR_L + BND_CURR_R);
    while (i < stride) {
      int rem = i;
      int jz = rem / (BS_Y + BND_CURR_L + BND_CURR_R);
      rem -= jz * (BS_Y + BND_CURR_L + BND_CURR_R);
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

template<typename BS>
class GCurr
{
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
// CurrmemShared

struct CurrmemShared
{
  template<typename BS>
  using Curr = SCurr<BS>;
  
  static Range<int> block_starts() { return range(4);  }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx =  (cmprts.b_mx()[1] + 1) / 2;
    int gy = ((cmprts.b_mx()[2] + 1) / 2) * cmprts.n_patches;
    return dim3(gx, gy);
  }

  template<typename BS>
  __device__ static int find_block_pos_patch(const DMparticlesCuda& dmprts, int *block_pos, int *ci0, int block_start)
  {
    return dmprts.find_block_pos_patch_q<BS::x::value, BS::y::value, BS::z::value>(block_pos, ci0, block_start);
  }

  __device__ static int find_bid(DMparticlesCuda& dmprts, int p, int *block_pos)
  {
    return dmprts.find_bid_q(p, block_pos);
  }
};

// ======================================================================
// CurrmemGlobal

struct CurrmemGlobal
{
  template<typename BS>
  using Curr = GCurr<BS>;

  static Range<int> block_starts() { return range(1);  }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx = cmprts.b_mx()[1];
    int gy = cmprts.b_mx()[2] * cmprts.n_patches;
    return dim3(gx, gy);
  }

  template<typename BS>
  __device__ static int find_block_pos_patch(const DMparticlesCuda& dmprts, int *block_pos, int *ci0, int block_start)
  {
    return dmprts.find_block_pos_patch<BS::x::value, BS::y::value, BS::z::value>(block_pos, ci0);
  }

  __device__ static int find_bid(DMparticlesCuda& dmprts, int p, int *block_pos)
  {
    return dmprts.find_bid();
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
curr_vb_cell(DMparticlesCuda& dmprts, int i[3], float x[3], float dx[3], float qni_wni,
	     CURR &scurr, int *ci0)
{
  float xa[3] = { 0.,
		 x[1] + .5f * dx[1],
		 x[2] + .5f * dx[2], };
  if (DEPOSIT == DEPOSIT_VB_3D) {
    if (dx[0] != 0.f) {
      float fnqx = qni_wni * dmprts.fnqxs();
      float h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
      scurr.add(0, i[1]  , i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), ci0);
      scurr.add(0, i[1]+1, i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), ci0);
      scurr.add(0, i[1]  , i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), ci0);
      scurr.add(0, i[1]+1, i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), ci0);
    }
  }
  if (dx[1] != 0.f) {
    float fnqy = qni_wni * dmprts.fnqys();
    scurr.add(1, i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]), ci0);
    scurr.add(1, i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]), ci0);
  }
  if (dx[2] != 0.f) {
    float fnqz = qni_wni * dmprts.fnqzs();
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
yz_calc_j(DMparticlesCuda& dmprts, struct d_particle& prt, int n, float4 *d_xi4, float4 *d_pxi4,
	  CURR &scurr,
	  int nr_total_blocks, int p_nr,
	  uint *d_bidx, int bid, int *ci0)
{
  AdvanceParticle<real_t, dim> advance{dmprts.dt()};

  float vxi[3];
  advance.calc_v(vxi, prt.pxi);

  // position xm at x^(n+.5)
  float h0[3], h1[3];
  float xm[3], xp[3];
  int j[3], k[3];
  
  dmprts.find_idx_off_pos_1st(prt.xi, j, h0, xm, float(0.));

  if (DEPOSIT == DEPOSIT_VB_2D) {
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    advance.push_x(prt.xi, vxi, .5f);

    float fnqx = vxi[0] * prt.qni_wni * dmprts.fnqs();

    // out-of-plane currents at intermediate time
    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.xi, lf, of, float(0.));
    lf[1] -= ci0[1];
    lf[2] -= ci0[2];

    scurr.add(0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnqx, ci0);
    scurr.add(0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnqx, ci0);
    scurr.add(0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnqx, ci0);
    scurr.add(0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnqx, ci0);

    // x^(n+1.0), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    advance.push_x(prt.xi, vxi, .5f);
    STORE_PARTICLE_POS(prt, d_xi4, n);
  } else if (DEPOSIT == DEPOSIT_VB_3D) {
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    advance.push_x(prt.xi, vxi);
    STORE_PARTICLE_POS(prt, d_xi4, n);
  }

  // has moved into which block? (given as relative shift)
  d_bidx[n] = dmprts.blockShift(prt.xi, p_nr, bid);

  // position xm at x^(n+.5)
  dmprts.find_idx_off_pos_1st(prt.xi, k, h1, xp, float(0.));

  // deposit xm -> xp
  int idiff[3] = { 0, k[1] - j[1], k[2] - j[2] };
  int i[3] = { 0, j[1] - ci0[1], j[2] - ci0[2] };
  float x[3] = { 0.f, xm[1] - j[1] - float(.5), xm[2] - j[2] - float(.5) };
  //float dx[3] = { 0.f, xp[1] - xm[1], xp[2] - xm[2] };
  float dx[3] = { vxi[0] * dmprts.dt() * dmprts.dxi(0), xp[1] - xm[1], xp[2] - xm[2] };
  
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
  curr_vb_cell<DEPOSIT, CURR>(dmprts, i, x, dx1, prt.qni_wni, scurr, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
  
  off[1] = idiff[1] - off[1];
  off[2] = idiff[2] - off[2];
  calc_dx1(dx1, x, dx, off);
  curr_vb_cell<DEPOSIT, CURR>(dmprts, i, x, dx1, prt.qni_wni, scurr, ci0);
  curr_vb_cell_upd(i, x, dx1, dx, off);
    
  curr_vb_cell<DEPOSIT, CURR>(dmprts, i, x, dx, prt.qni_wni, scurr, ci0);
}

// ----------------------------------------------------------------------
// push_mprts_ab

template<typename Config, bool REORDER, typename OPT_IP, int DEPOSIT>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(int block_start, DMparticlesCuda dmprts, DMFields d_mflds)
{
  using BS = typename Config::Bs;
  using Currmem = typename Config::Currmem;
  using Curr = typename Currmem::Curr<BS>;
  using FldCache_t = FldCache<BS::x::value, BS::y::value, BS::z::value>;

  int block_pos[3], ci0[3];
  int p = Currmem::template find_block_pos_patch<BS>(dmprts, block_pos, ci0, block_start);
  if (p < 0)
    return;

  int bid = Currmem::find_bid(dmprts, p, block_pos);
  int block_begin = dmprts.off_[bid];
  int block_end = dmprts.off_[bid + 1];
  __shared__ FldCache_t fld_cache;
  fld_cache.load(d_mflds[p], ci0);
  __shared__ float _scurr[Curr::shared_size];
  Curr scurr(_scurr, d_mflds[p]);
  
  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BS::x::value, BS::y::value, BS::z::value, REORDER, OPT_IP>
      (dmprts, prt, n, fld_cache);

    if (REORDER) {
      yz_calc_j<DEPOSIT_VB_2D, Curr>
	(dmprts, prt, n, dmprts.alt_xi4_, dmprts.alt_pxi4_, scurr, dmprts.n_blocks_, p, dmprts.bidx_, bid, ci0);
    } else {
      yz_calc_j<DEPOSIT_VB_2D, Curr>
	(dmprts, prt, n, dmprts.xi4_, dmprts.pxi4_, scurr, dmprts.n_blocks_, p, dmprts.bidx_, bid, ci0);
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

template<typename Config>
template<bool REORDER>
void CudaPushParticles_<Config>::push_mprts_ab(CudaMparticles* cmprts, struct cuda_mfields *cmflds)
{
  using Currmem = typename Config::Currmem;

  zero_currents(cmflds);

  dim3 dimGrid = Currmem::dimGrid(*cmprts);

  if (REORDER) {
    cmprts->d_alt_xi4.resize(cmprts->n_prts);
    cmprts->d_alt_pxi4.resize(cmprts->n_prts);
  }

  for (auto block_start : Currmem::block_starts()) {
    ::push_mprts_ab<Config, REORDER, typename Config::Ip, Config::Deposit::value>
      <<<dimGrid, THREADS_PER_BLOCK>>>(block_start, *cmprts, *cmflds);
    cuda_sync_if_enabled();
  }

  if (REORDER) {
    cmprts->swap_alt();
    cmprts->need_reorder = false;
  }
}

// ----------------------------------------------------------------------
// push_mprts_yz

template<typename Config>
void CudaPushParticles_<Config>::push_mprts_yz(CudaMparticles* cmprts, struct cuda_mfields *cmflds)
{
  if (!cmprts->need_reorder) {
    //    printf("INFO: yz_cuda_push_mprts: need_reorder == false\n");
    push_mprts_ab<false>(cmprts, cmflds);
  } else {
    push_mprts_ab<true>(cmprts, cmflds);
  }
}

template struct CudaPushParticles_<Config1vb>;
template struct CudaPushParticles_<Config1vbec3d>;
