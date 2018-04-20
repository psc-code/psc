
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

#include "cuda_fld_cache.cuh"
#include "cuda_currmem.cuh"

using real_t = float;

// FIXME
#define CUDA_BND_S_OOB (10)

// ----------------------------------------------------------------------

// OPT: use more shmem?

// OPT: passing shared memory cache etc around is probably sub-optimal
// OPT: fld cache is much bigger than needed
// OPT: precalculating IP coeffs could be a gain, too

template<typename Config, typename FC>
struct CudaPushParticles_yz
{
  using BS = typename Config::Bs;
  using DMparticles = DMparticlesCuda<BS>;
  using FldCache = FC;

  // ----------------------------------------------------------------------
  // push_part_one

  template<bool REORDER>
  __device__ static void
  push_part_one(DMparticles& dmprts, struct d_particle& prt, int n, const FldCache& fld_cache)
    
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
    InterpolateEM<FldCache, typename Config::Ip, dim_yz> ip;
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
  
  template<typename CURR>
  __device__ static void
  curr_vb_cell(DMparticlesCuda<BS144>& dmprts, int i[3], float x[3], float dx[3], float qni_wni,
	       CURR &scurr, int *ci0)
  {
    float xa[3] = { 0.,
		    x[1] + .5f * dx[1],
		    x[2] + .5f * dx[2], };
    if (Config::Deposit::value == DEPOSIT_VB_3D) {
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
  
  template<typename CURR>
  __device__ static void
  yz_calc_j(DMparticlesCuda<BS144>& dmprts, struct d_particle& prt, int n, float4 *d_xi4, float4 *d_pxi4,
	    CURR &scurr, int p_nr, int bid, int *ci0)
  {
    AdvanceParticle<real_t, dim> advance{dmprts.dt()};

    float vxi[3];
    advance.calc_v(vxi, prt.pxi);

    // position xm at x^(n+.5)
    float h0[3], h1[3];
    float xm[3], xp[3];
    int j[3], k[3];
  
    dmprts.find_idx_off_pos_1st(prt.xi, j, h0, xm, float(0.));

    if (Config::Deposit::value == DEPOSIT_VB_2D) {
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
    } else if (Config::Deposit::value == DEPOSIT_VB_3D) {
      // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
      advance.push_x(prt.xi, vxi);
      STORE_PARTICLE_POS(prt, d_xi4, n);
    }

    // has moved into which block? (given as relative shift)
    dmprts.bidx_[n] = dmprts.blockShift(prt.xi, p_nr, bid);

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
    curr_vb_cell(dmprts, i, x, dx1, prt.qni_wni, scurr, ci0);
    curr_vb_cell_upd(i, x, dx1, dx, off);
  
    off[1] = idiff[1] - off[1];
    off[2] = idiff[2] - off[2];
    calc_dx1(dx1, x, dx, off);
    curr_vb_cell(dmprts, i, x, dx1, prt.qni_wni, scurr, ci0);
    curr_vb_cell_upd(i, x, dx1, dx, off);
    
    curr_vb_cell(dmprts, i, x, dx, prt.qni_wni, scurr, ci0);
  }

  // ----------------------------------------------------------------------
  // push_mprts

  template<bool REORDER>
  __device__
  static void push_mprts(DMparticles& dmprts, DMFields& d_mflds, int block_start)
  {
    using BS = typename Config::Bs;
    using Currmem = typename Config::Currmem;
    using Curr = typename Currmem::Curr<BS>;
    
    int block_pos[3], ci0[3];
    int p = Currmem::template find_block_pos_patch<BS>(dmprts, block_pos, ci0, block_start);
    if (p < 0)
      return;
    
    int bid = Currmem::find_bid(dmprts, p, block_pos);
    int block_begin = dmprts.off_[bid];
    int block_end = dmprts.off_[bid + 1];
    
    using CURR = SCurr<BS>;

    __shared__ FldCache fld_cache;
    fld_cache.load(d_mflds[p], ci0);

    __shared__ float _scurr[CURR::shared_size];
    CURR scurr(_scurr, d_mflds[p]);

    __syncthreads();
    for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
      if (n < block_begin) {
	continue;
      }
      struct d_particle prt;
      push_part_one<REORDER>(dmprts, prt, n, fld_cache);
      
      if (REORDER) {
	yz_calc_j(dmprts, prt, n, dmprts.alt_xi4_, dmprts.alt_pxi4_, scurr, p, bid, ci0);
      } else {
	yz_calc_j(dmprts, prt, n, dmprts.xi4_, dmprts.pxi4_, scurr, p, bid, ci0);
      }
    }
    
    scurr.add_to_fld(ci0);
  }
};

// ----------------------------------------------------------------------
// push_mprts_ab

template<typename Config, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(int block_start, DMparticlesCuda<BS144> dmprts, DMFields d_mflds)
{
  using BS = typename Config::Bs;
  using FldCache_t = FldCache<BS::x::value, BS::y::value, BS::z::value>;
  
  CudaPushParticles_yz<Config, FldCache_t>::push_mprts<REORDER>(dmprts, d_mflds, block_start);
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
    ::push_mprts_ab<Config, REORDER>
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
