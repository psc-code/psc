
#include "cuda_iface.h"
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_push_particles.cuh"
#include "push_particles_cuda_impl.hxx"
#include "range.hxx"

#include "../psc_push_particles/inc_defs.h"

#include "psc.h" // FIXME

#include "dim.hxx"

#include "interpolate.hxx"
#include "pushp.hxx"

#define BND (2) // FIXME

#define THREADS_PER_BLOCK (512)

#include "cuda_fld_cache.cuh"
#include "cuda_currmem.cuh"

// FIXME
#define CUDA_BND_S_OOB (10)

// ----------------------------------------------------------------------

// OPT: use more shmem?

// OPT: passing shared memory cache etc around is probably sub-optimal
// OPT: fld cache is much bigger than needed
// OPT: precalculating IP coeffs could be a gain, too

template<typename Config, bool REORDER>
struct CudaPushParticles
{
  using dim = typename Config::dim;
  using BS = typename Config::Bs;
  using Currmem = typename Config::Currmem;
  using Block = typename Currmem::Block<BS, dim>;
  using Curr = typename Currmem::Curr<BS, dim>;
  using DMparticles = DMparticlesCuda<BS>;
  using real_t = typename DMparticles::real_t;
  using FldCache = FldCache<BS, dim>;

  // ----------------------------------------------------------------------
  // push_part_one

  __device__ static void
  push_part_one(DMparticles& dmprts, struct d_particle& prt, int n, const FldCache& fld_cache,
		const Block& current_block)
    
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
    dmprts.template scalePos<dim>(xm, prt.xi);
#if 1
    const int *ci0 = current_block.ci0;
    if (!dim::InvarX::value && ((xm[0] < ci0[0] || xm[0] > ci0[0] + BS::x::value)) ||
	xm[1] < ci0[1] || xm[1] > ci0[1] + BS::y::value ||
	xm[2] < ci0[2] || xm[2] > ci0[2] + BS::z::value) {
      printf("xm %g %g (xi %g %g n %d)\n", xm[1], xm[2], prt.xi[0], prt.xi[1], n);
    }
#endif
    InterpolateEM<FldCache, typename Config::Ip, dim> ip;
    AdvanceParticle<real_t, dim> advance{dmprts.dt()};
    
    ip.set_coeffs(xm);
    
    real_t E[3] = { ip.ex(fld_cache), ip.ey(fld_cache), ip.ez(fld_cache) };
    real_t H[3] = { ip.hx(fld_cache), ip.hy(fld_cache), ip.hz(fld_cache) };
#if 1
    if (!isfinite(E[0]) || !isfinite(E[1]) || !isfinite(E[2]) ||
	!isfinite(H[0]) || !isfinite(H[1]) || !isfinite(H[2])) {
      printf("CUDA_ERROR push_part_one: n = %d E %g %g %g H %g %g %g\n", n,
	     E[0], E[1], E[2], H[0], H[1], H[2]);
      printf("CUDA_ERROR xm %g %g\n", xm[1], xm[2]);
    }
#endif
    
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
#if 0
    if (!isfinite(prt.pxi[0]) || !isfinite(prt.pxi[1]) || !isfinite(prt.pxi[2])) {
      printf("CUDA_ERROR push_part_one: n = %d pxi %g %g %g\n", n,
	     prt.pxi[0], prt.pxi[1], prt.pxi[2]);
    }
#endif
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
  // curr_vb_cell -- dim_yz
  
  __device__ static void
  curr_vb_cell(DMparticles& dmprts, int i[3], float x[3], float dx[3], float qni_wni,
	       Curr &scurr, const Block& current_block, dim_yz tag)
  {
#if 0
    if (i[1] < -1 || i[1] >= int(BS::y::value) + 1 ||
	i[2] < -1 || i[2] >= int(BS::z::value) + 1) {
      printf("CUDA_ERROR curr_vb_cell jyz %d:%d:%d\n", i[1], i[2]);
    }
#endif
    float xa[3] = { 0.,
		    x[1] + .5f * dx[1],
		    x[2] + .5f * dx[2], };
    if (Config::Deposit::value == DEPOSIT_VB_3D) {
      if (dx[0] != 0.f) {
	float fnqx = qni_wni * dmprts.fnqxs();
	float h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
	scurr.add(0, i[1]  , i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), current_block.ci0);
	scurr.add(0, i[1]+1, i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), current_block.ci0);
	scurr.add(0, i[1]  , i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), current_block.ci0);
	scurr.add(0, i[1]+1, i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), current_block.ci0);
      }
    }
    if (dx[1] != 0.f) {
      float fnqy = qni_wni * dmprts.fnqys();
      scurr.add(1, i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]), current_block.ci0);
      scurr.add(1, i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]), current_block.ci0);
    }
    if (dx[2] != 0.f) {
      float fnqz = qni_wni * dmprts.fnqzs();
      scurr.add(2, i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]), current_block.ci0);
      scurr.add(2, i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]), current_block.ci0);
    }
  }

  // ----------------------------------------------------------------------
  // curr_vb_cell -- dim_xyz
  
  __device__ static void
  curr_vb_cell(DMparticles& dmprts, int i[3], float x[3], float dx[3], float qni_wni,
	       Curr &scurr, const Block& current_block, dim_xyz tag)
  {
#if 1
    if (i[0] < -1 || i[0] >= int(BS::x::value) + 1 ||
	i[1] < -1 || i[1] >= int(BS::y::value) + 1 ||
	i[2] < -1 || i[2] >= int(BS::z::value) + 1) {
      printf("CUDA_ERROR curr_vb_cell jxyz %d:%d:%d\n", i[0], i[1], i[2]);
    }
#endif
    float xa[3] = { x[0] + .5f * dx[0],
		    x[1] + .5f * dx[1],
		    x[2] + .5f * dx[2], };
    float h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
    if (dx[0] != 0.f) {
      float fnqx = qni_wni * dmprts.fnqxs();
      scurr.add(0, i[0], i[1]  , i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), current_block.ci0);
      scurr.add(0, i[0], i[1]+1, i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), current_block.ci0);
      scurr.add(0, i[0], i[1]  , i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), current_block.ci0);
      scurr.add(0, i[0], i[1]+1, i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), current_block.ci0);
    }
    if (dx[1] != 0.f) {
      float fnqy = qni_wni * dmprts.fnqys();
      scurr.add(1, i[0]  , i[1], i[2]  , fnqy * (dx[1] * (.5f - xa[0]) * (.5f - xa[2]) + h), current_block.ci0);
      scurr.add(1, i[0]+1, i[1], i[2]  , fnqy * (dx[1] * (.5f + xa[0]) * (.5f - xa[2]) - h), current_block.ci0);
      scurr.add(1, i[0]  , i[1], i[2]+1, fnqy * (dx[1] * (.5f - xa[0]) * (.5f + xa[2]) + h), current_block.ci0);
      scurr.add(1, i[0]+1, i[1], i[2]+1, fnqy * (dx[1] * (.5f + xa[0]) * (.5f + xa[2]) - h), current_block.ci0);
    }
    if (dx[2] != 0.f) {
      float fnqz = qni_wni * dmprts.fnqzs();
      scurr.add(2, i[0]  , i[1]  , i[2], fnqz * (dx[2] * (.5f - xa[0]) * (.5f - xa[1]) + h), current_block.ci0);
      scurr.add(2, i[0]+1, i[1]  , i[2], fnqz * (dx[2] * (.5f + xa[0]) * (.5f - xa[1]) - h), current_block.ci0);
      scurr.add(2, i[0]  , i[1]+1, i[2], fnqz * (dx[2] * (.5f - xa[0]) * (.5f + xa[1]) + h), current_block.ci0);
      scurr.add(2, i[0]+1, i[1]+1, i[2], fnqz * (dx[2] * (.5f + xa[0]) * (.5f + xa[1]) - h), current_block.ci0);
    }
  }

  // ----------------------------------------------------------------------
  // curr_vb_cell_upd
  
  __device__ static void
  curr_vb_cell_upd(int i[3], float x[3], float dx1[3], float dx[3], int off[3], dim_yz tag)
  {
    dx[0] -= dx1[0];
    dx[1] -= dx1[1];
    dx[2] -= dx1[2];
    x[1] += dx1[1] - off[1];
    x[2] += dx1[2] - off[2];
    i[1] += off[1];
    i[2] += off[2];
#if 0
    if (i[1] < -1 || i[1] >= int(BS::y::value) + 1 ||
	i[2] < -1 || i[2] >= int(BS::z::value) + 1) {
      printf("CUDA_ERROR cell_upd B jyz %d:%d\n", i[1], i[2]);
    }
#endif
  }
  
  // ----------------------------------------------------------------------
  // calc_j -- dispatched for dim_yz
  
  __device__ static void
  calc_j(DMparticles& dmprts, struct d_particle& prt, int n, float4 *d_xi4, float4 *d_pxi4,
	 Curr &scurr, const Block& current_block, dim_yz tag)
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
      lf[1] -= current_block.ci0[1];
      lf[2] -= current_block.ci0[2];

      scurr.add(0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnqx, current_block.ci0);
      scurr.add(0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnqx, current_block.ci0);
      scurr.add(0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnqx, current_block.ci0);
      scurr.add(0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnqx, current_block.ci0);

      // x^(n+1.0), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
      advance.push_x(prt.xi, vxi, .5f);
      STORE_PARTICLE_POS(prt, d_xi4, n);
    } else if (Config::Deposit::value == DEPOSIT_VB_3D) {
      // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
      advance.push_x(prt.xi, vxi);
      STORE_PARTICLE_POS(prt, d_xi4, n);
    }

    // has moved into which block? (given as relative shift)
    dmprts.bidx_[n] = dmprts.blockShift(prt.xi, current_block.p, current_block.bid);

    // position xm at x^(n+.5)
    dmprts.find_idx_off_pos_1st(prt.xi, k, h1, xp, float(0.));

    // deposit xm -> xp
    int idiff[3] = { 0, k[1] - j[1], k[2] - j[2] };
#if 0
    if (idiff[1] < -1 || idiff[1] > 1 ||
	idiff[2] < -1 || idiff[2] > 1) {
      printf("A idiff %d %d j %d %d k %d %d\n", idiff[1], idiff[2],
	     j[1], j[2], k[1], k[2]);
      printf("A prt.xi %g %g scaled %g %g k %d %d\n", prt.xi[1], prt.xi[2],
	     dmprts.scalePos(prt.xi[1], 1), dmprts.scalePos(prt.xi[2], 2), k[1], k[2]);
    }
#endif
    int i[3] = { 0, j[1] - current_block.ci0[1], j[2] - current_block.ci0[2] };
#if 0
    if (i[1] < -1 || i[1] >= int(BS::y::value) + 1 ||
	i[2] < -1 || i[2] >= int(BS::z::value) + 1) {
      printf("CUDA_ERROR deposit jyz %d:%d\n", i[1], i[2]);
    }
#endif
    float x[3] = { 0.f, xm[1] - j[1] - float(.5), xm[2] - j[2] - float(.5) };
    //float dx[3] = { 0.f, xp[1] - xm[1], xp[2] - xm[2] };
    float dx[3] = { dmprts.scalePos(vxi[0] * dmprts.dt(), 0), xp[1] - xm[1], xp[2] - xm[2] };
  
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
    curr_vb_cell(dmprts, i, x, dx1, prt.qni_wni, scurr, current_block, dim{});
    curr_vb_cell_upd(i, x, dx1, dx, off, dim{});
  
    off[1] = idiff[1] - off[1];
    off[2] = idiff[2] - off[2];
    calc_dx1(dx1, x, dx, off);
    curr_vb_cell(dmprts, i, x, dx1, prt.qni_wni, scurr, current_block, dim{});
    curr_vb_cell_upd(i, x, dx1, dx, off, dim{});
    
    curr_vb_cell(dmprts, i, x, dx, prt.qni_wni, scurr, current_block, dim{});
  }

  // ----------------------------------------------------------------------
  // calc_j -- dispatched for dim_xyz
  
  __device__ static void
  calc_j(DMparticles& dmprts, struct d_particle& prt, int n, float4 *d_xi4, float4 *d_pxi4,
	 Curr &scurr, const Block& current_block, dim_xyz tag)
  {
    AdvanceParticle<real_t, dim> advance{dmprts.dt()};

    float vxi[3];
    advance.calc_v(vxi, prt.pxi);

    // position xm at x^(n+.5)
    float h0[3];
    float xm[3], xp[3];
    int j[3];

    dmprts.find_idx_off_pos_1st(prt.xi, j, h0, xm, float(0.));

    static_assert(Config::Deposit::value == DEPOSIT_VB_3D, "calc_j dim_xyz needs 3d deposit");
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    advance.push_x(prt.xi, vxi);
    STORE_PARTICLE_POS(prt, d_xi4, n);

    // position xp at x^(n+1.5)
    for (int d = 0; d < 3; d++) {
      xp[d] = dmprts.scalePos(prt.xi[d], d);
    }

    // deposit xm -> xp
    int i[3] = { j[0] - current_block.ci0[0], j[1] - current_block.ci0[1], j[2] - current_block.ci0[2] };
#if 1
    if (i[0] < 0 || i[0] >= int(BS::x::value) ||
	i[1] < 0 || i[1] >= int(BS::y::value) ||
	i[2] < 0 || i[2] >= int(BS::z::value)) {
      printf("CUDA_ERROR deposit jxyz %d:%d:%d\n", i[0], i[1], i[2]);
      const int* ci0 = current_block.ci0;
      printf("j %d:%d:%d ci0 %d:%d:%d xm %g:%g:%g\n", j[0], j[1], j[2], ci0[0], ci0[1], ci0[2],
	     xm[0], xm[1], xm[2]);
    }
#endif
    // x will be orig position inside cell, ie [-.5,.5]
    float x[3] = { xm[0] - (j[0] + float(.5)),
		   xm[1] - (j[1] + float(.5)),
		   xm[2] - (j[2] + float(.5)) };
    //FIXME x and h0 are almost the same (shifted by .5)
    //printf("x %g %g %g h0 %g %g %g\n", x[0], x[1], x[2], h0[0], h0[1], h0[2]);
    float dx[3] = { xp[0] - xm[0],
		    xp[1] - xm[1],
		    xp[2] - xm[2] };
    for (;;) {
      float bnd[3] = { dx[0] > 0 ? 1.f : -1.f,
		       dx[1] > 0 ? 1.f : -1.f,
		       dx[2] > 0 ? 1.f : -1.f };
      float frac[3] = { dx[0] == 0.f ? 1.f : (.5f * bnd[0] - x[0]) / dx[0],
			dx[1] == 0.f ? 1.f : (.5f * bnd[1] - x[1]) / dx[1],
			dx[2] == 0.f ? 1.f : (.5f * bnd[2] - x[2]) / dx[2] };
      float step = 1.f;
      int dir = -1;
      if (frac[0] < step) { step = frac[0]; dir = 0; }
      if (frac[1] < step) { step = frac[1]; dir = 1; }
      if (frac[2] < step) { step = frac[2]; dir = 2; }
      
      float dx1[3] = { dx[0] * step, dx[1] * step, dx[2] * step };
      // printf("frac %g %g %g step %g dir %d\n", frac[0], frac[1], frac[2], step, dir);
      // printf("i %d:%d:%d dx1 %g %g %g x %g %g %g\n", i[0], i[1], i[2], dx1[0], dx1[1], dx1[2],
      // 	     x[0], x[1], x[2]);
      curr_vb_cell(dmprts, i, x, dx1, prt.qni_wni, scurr, current_block, dim{});
      if (dir < 0) {
	break;
      }
      
      i[dir] += dx[dir] > 0 ? 1 : -1;
      for (int d = 0; d < 3; d++) {
	x[d] += dx1[d];
	dx[d] -= dx1[d];
	x[dir] = - .5f * bnd[dir];
      }
      // printf("2 i %d:%d:%d dx %g %g %g x %g %g %g\n", i[0], i[1], i[2],
      // 	     dx[0], dx[1], dx[2], x[0], x[1], x[2]); 
    }

    // has moved into which block?
#if 0
    int cpos[3] = { (i[0] + current_block.ci0[0]),
		    (i[1] + current_block.ci0[1]),
		    (i[2] + current_block.ci0[2]) };
#else
    int cpos[3] = { __float2int_rd(xp[0]),
		    __float2int_rd(xp[1]),
		    __float2int_rd(xp[2]) };
#endif
    //dmprts.find_idx_off_pos_1st(prt.xi, i, h0, xp, float(0.));
    // printf("cpos %d %d %d -- i %d %d %d\n", cpos[0], cpos[1], cpos[2],
    // 	   __float2int_rd(xp[0]), __float2int_rd(xp[1]), __float2int_rd(xp[2]));

    dmprts.bidx_[n] = dmprts.blockIndexFromCellPosition(cpos, current_block.p);
    //    printf("cpos %d %d %d p %d bidx %d\n", cpos[0], cpos[1], cpos[2], current_block.p, dmprts.bidx_[n]);
  }

  // ----------------------------------------------------------------------
  // push_mprts

  __device__
  static void push_mprts(DMparticles& dmprts, DMFields& d_mflds, int block_start)
  {
    Block current_block;
    if (!current_block.init(dmprts, block_start)) {
      return;
    }
    
    __shared__ FldCache fld_cache;
    fld_cache.load(d_mflds[current_block.p], current_block.ci0);

    __shared__ float _scurr[Curr::shared_size];
    Curr scurr(_scurr, d_mflds[current_block.p]);
    __syncthreads();

    int block_begin = dmprts.off_[current_block.bid];
    int block_end = dmprts.off_[current_block.bid + 1];
    for (int n : in_block_loop(block_begin, block_end)) {
      if (n < block_begin) {
	continue;
      }
      struct d_particle prt;
      push_part_one(dmprts, prt, n, fld_cache, current_block);
      
      if (REORDER) {
	calc_j(dmprts, prt, n, dmprts.alt_xi4_, dmprts.alt_pxi4_, scurr, current_block, dim{});
      } else {
	calc_j(dmprts, prt, n, dmprts.xi4_, dmprts.pxi4_, scurr, current_block, dim{});
      }
    }
    
    scurr.add_to_fld(current_block.ci0);
  }
};

// ----------------------------------------------------------------------
// push_mprts_ab

template<typename Config, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(int block_start, DMparticlesCuda<typename Config::Bs> dmprts, DMFields d_mflds)
{
  CudaPushParticles<Config, REORDER>::push_mprts(dmprts, d_mflds, block_start);
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
  using Block = typename Currmem::Block<typename Config::Bs, typename Config::dim>;

  zero_currents(cmflds);

  dim3 dimGrid = Block::dimGrid(*cmprts);

  if (REORDER) {
    cmprts->d_alt_xi4.resize(cmprts->n_prts);
    cmprts->d_alt_pxi4.resize(cmprts->n_prts);
  }

  for (auto block_start : Block::block_starts()) {
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
// push_mprts

template<typename Config>
void CudaPushParticles_<Config>::push_mprts(CudaMparticles* cmprts, struct cuda_mfields *cmflds)
{
  if (!cmprts->need_reorder) {
    //    printf("INFO: push_mprts: need_reorder == false\n");
    push_mprts_ab<false>(cmprts, cmflds);
  } else {
    push_mprts_ab<true>(cmprts, cmflds);
  }
}

//template struct CudaPushParticles_<CudaConfig1vb<dim_yz>>;
template struct CudaPushParticles_<CudaConfig1vbec3d<dim_yz, BS144>>;
template struct CudaPushParticles_<CudaConfig1vbec3dGmem<dim_xyz, BS144>>;
template struct CudaPushParticles_<CudaConfig1vbec3dGmem<dim_xyz, BS444>>;
