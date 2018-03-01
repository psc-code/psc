
namespace {

#if DIM == DIM_XYZ

#include "inc_curr_1vb_split.c" // FIXME hack...

#else

// ======================================================================
// 1vb current deposition "var1"
//
// my original implementation, 2d (yz) only
// in theory less divergent, in particular for CUDA,
// but also definitely more complex

// ----------------------------------------------------------------------
// calc_3d_dx1

template<typename curr_cache_t, typename dim_t>
struct Current1vb
{
  using real_t = typename curr_cache_t::real_t;
  using Real3 = Vec3<real_t>;
  
  Current1vb(const Grid_t& grid)
    : dt_(grid.dt),
      dxi_{ Real3{1., 1. , 1.} / Real3(grid.dx) }
  {
    fnqxs_ = grid.dx[0] * grid.fnqs / grid.dt;
    fnqys_ = grid.dx[1] * grid.fnqs / grid.dt;
    fnqzs_ = grid.dx[2] * grid.fnqs / grid.dt;
  }
  
  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, particle_t *prt, real_t *vxi, dim_1 tag)
  {
    // FIXME
    //assert(0);
  }

  CUDA_DEVICE static inline void
  calc_3d_dx1(real_t dx1[3], real_t x[3], real_t dx[3], int off[3])
  {
#ifdef __CUDACC__
    if (off[1] == 0) {
      if (off[2] == 0 || dx[2] == 0.f) {
	dx1[0] = 0.f;
	dx1[1] = 0.f;
	dx1[2] = 0.f;
      } else {
	dx1[2] = .5f * off[2] - x[2];
	dx1[1] = dx[1] / dx[2] * dx1[2];
	dx1[0] = dx[0] / dx[2] * dx1[2];
      }
    } else { // off[1] != 0
      if (dx[1] == 0.f) {
	dx1[0] = 0.f;
	dx1[1] = 0.f;
	dx1[2] = 0.f;
      } else {
	dx1[1] = .5f * off[1] - x[1];
	dx1[2] = dx[2] / dx[1] * dx1[1];
	dx1[0] = dx[0] / dx[1] * dx1[1];
      }
    }

#else

    if (off[2] == 0) {
      dx1[1] = .5f * off[1] - x[1];
      if (dx[1] == 0.f) {
	dx1[0] = 0.f;
	dx1[2] = 0.f;
      } else {
	dx1[0] = dx[0] / dx[1] * dx1[1];
	dx1[2] = dx[2] / dx[1] * dx1[1];
      }
    } else {
      dx1[2] = .5f * off[2] - x[2];
      if (dx[2] == 0.f) {
	dx1[0] = 0.f;
	dx1[1] = 0.f;
      } else {
	dx1[0] = dx[0] / dx[2] * dx1[2];
	dx1[1] = dx[1] / dx[2] * dx1[2];
      }
    }
#endif
  }

  // ----------------------------------------------------------------------
  // curr_3d_vb_cell

  void curr_3d_vb_cell(curr_cache_t curr_cache, int i[3], real_t x[3], real_t dx[3],
		       real_t qni_wni)
  {
    real_t xa[3] = { 0.,
		     x[1] + .5f * dx[1],
		     x[2] + .5f * dx[2], };
#ifdef __CUDACC__
    if (dx[0] != 0.f)
#endif
      {
	real_t fnqx = qni_wni * fnqxs_;
	real_t h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
	curr_cache.add(0, 0,i[1]  ,i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h));
	curr_cache.add(0, 0,i[1]+1,i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h));
	curr_cache.add(0, 0,i[1]  ,i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) - h));
	curr_cache.add(0, 0,i[1]+1,i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) + h));
      }
#ifdef __CUDACC__
    if (dx[1] != 0.f)
#endif
      {
	real_t fnqy = qni_wni * fnqys_;
	curr_cache.add(1, 0,i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]));
	curr_cache.add(1, 0,i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]));
      }
#ifdef __CUDACC__
    if (dx[2] != 0.f)
#endif
      {
	real_t fnqz = qni_wni * fnqzs_;
	curr_cache.add(2, 0,i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]));
	curr_cache.add(2, 0,i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]));
      }
  }

  // ----------------------------------------------------------------------
  // curr_3d_vb_cell_upd

  CUDA_DEVICE static void
  curr_3d_vb_cell_upd(int i[3], real_t x[3], real_t dx1[3],
		      real_t dx[3], int off[3])
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

  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, particle_t *prt, real_t *vxi, dim_yz tag_dim)
  {
    // deposit xm -> xp
    int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			
    int i[3] = { 0, lg[1], lg[2] };					
    real_t dx[3] = { vxi[0] * dt_ * dxi_[0], xp[1] - xm[1], xp[2] - xm[2] };
    real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; 
    
    real_t dx1[3];
    int off[3];
#ifdef __CUDACC__
    
    real x1 = x[1] * idiff[1];
    real x2 = x[2] * idiff[2];
    int d_first = (std::abs(dx[2]) * (.5f - x1) >= std::abs(dx[1]) * (.5f - x2));
    
    if (d_first == 0) {
      off[1] = idiff[1];
      off[2] = 0;
    } else {
      off[1] = 0;
      off[2] = idiff[2];
    }
    
    calc_3d_dx1(dx1, x, dx, off);
    curr_3d_vb_cell(curr_cache, i, x, dx1, prt->qni_wni);
    curr_3d_vb_cell_upd(i, x, dx1, dx, off);
    
    off[1] = idiff[1] - off[1];
    off[2] = idiff[2] - off[2];
    calc_3d_dx1(dx1, x, dx, off);
    curr_3d_vb_cell(curr_cache, i, x, dx1, prt->qni_wni);
    curr_3d_vb_cell_upd(i, x, dx1, dx, off);
    
    curr_3d_vb_cell(curr_cache, i, x, dx, prt->qni_wni);
    
#else
    int first_dir, second_dir = -1;
    /* FIXME, make sure we never div-by-zero? */
    if (idiff[1] == 0 && idiff[2] == 0) {
      first_dir = -1;
    } else if (idiff[1] == 0) {
      first_dir = 2;
    } else if (idiff[2] == 0) {
      first_dir = 1;
    } else {
      dx1[1] = .5f * idiff[1] - x[1];
      if (dx[1] == 0.f) {
	dx1[2] = 0.f;
      } else {
	dx1[2] = dx[2] / dx[1] * dx1[1];
      }
      if (std::abs(x[2] + dx1[2]) > .5f) {
	first_dir = 2;
      } else {
	first_dir = 1;
      }
      second_dir = 3 - first_dir;
    }
    real_t qni_wni = particle_qni_wni(prt);
    
    if (first_dir >= 0) {
      off[3 - first_dir] = 0;
      off[first_dir] = idiff[first_dir];
      calc_3d_dx1(dx1, x, dx, off);
      curr_3d_vb_cell(curr_cache, i, x, dx1, qni_wni);
      curr_3d_vb_cell_upd(i, x, dx1, dx, off);
    }
    
    if (second_dir >= 0) {
      off[first_dir] = 0;
      off[second_dir] = idiff[second_dir];
      calc_3d_dx1(dx1, x, dx, off);
      curr_3d_vb_cell(curr_cache, i, x, dx1, qni_wni);
      curr_3d_vb_cell_upd(i, x, dx1, dx, off);
    }
    
    curr_3d_vb_cell(curr_cache, i, x, dx, qni_wni);
#endif
  }

  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, particle_t *prt, real_t *vxi)
  {
#if DIM == DIM_1
    calc_j(curr_cache, xm, xp, lf, lg, prt, vxi, dim_1{});
#elif DIM == DIM_YZ
    calc_j(curr_cache, xm, xp, lf, lg, prt, vxi, dim_yz{});
#endif
  }
  
private:
  real_t dt_;
  Real3 dxi_;
  real_t fnqxs_, fnqys_, fnqzs_;
};

  
#endif

// ======================================================================
// TBD: piece to save block_idx as we go for following sort

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

}
