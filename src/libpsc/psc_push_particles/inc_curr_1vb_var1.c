
// ======================================================================
// 1vb current deposition "var1"
//
// my original implementation, 2d (yz) only
// in theory less divergent, in particular for CUDA,
// but also definitely more complex

// ----------------------------------------------------------------------
// calc_3d_dx1

template<typename CURR_CACHE, typename dim_t>
struct Current1vbVar1
{
  using curr_cache_t = CURR_CACHE;
  using real_t = typename curr_cache_t::real_t;
  using Real3 = Vec3<real_t>;
  
  Current1vbVar1(const Grid_t& grid)
    : dt_(grid.dt),
      dxi_{ Real3{1., 1. , 1.} / Real3(grid.domain.dx) }
  {
    fnqxs_ = grid.domain.dx[0] * grid.norm.fnqs / grid.dt;
    fnqys_ = grid.domain.dx[1] * grid.norm.fnqs / grid.dt;
    fnqzs_ = grid.domain.dx[2] * grid.norm.fnqs / grid.dt;
  }
  
  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, real_t qni_wni, real_t *vxi, dim_1 tag)
  {
    // FIXME
    //assert(0);
  }

  static void calc_3d_dx1(real_t dx1[3], real_t x[3], real_t dx[3], int off[3])
  {
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
  }

  // ----------------------------------------------------------------------
  // curr_3d_vb_cell

  void curr_3d_vb_cell(curr_cache_t curr_cache, int i[3], real_t x[3], real_t dx[3],
		       real_t qni_wni)
  {
    real_t xa[3] = { 0.,
		     x[1] + .5f * dx[1],
		     x[2] + .5f * dx[2], };

    real_t fnqx = qni_wni * fnqxs_;
    real_t h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
    curr_cache.add(0, 0,i[1]  ,i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h));
    curr_cache.add(0, 0,i[1]+1,i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h));
    curr_cache.add(0, 0,i[1]  ,i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) - h));
    curr_cache.add(0, 0,i[1]+1,i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) + h));

    real_t fnqy = qni_wni * fnqys_;
    curr_cache.add(1, 0,i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]));
    curr_cache.add(1, 0,i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]));

    real_t fnqz = qni_wni * fnqzs_;
    curr_cache.add(2, 0,i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]));
    curr_cache.add(2, 0,i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]));
  }

  // ----------------------------------------------------------------------
  // curr_3d_vb_cell_upd

  static void curr_3d_vb_cell_upd(int i[3], real_t x[3], real_t dx1[3],
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
	      int *lf, int *lg, real_t qni_wni, real_t *vxi, dim_yz tag_dim)
  {
    // deposit xm -> xp
    int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			
    int i[3] = { 0, lg[1], lg[2] };					
    real_t dx[3] = { vxi[0] * dt_ * dxi_[0], xp[1] - xm[1], xp[2] - xm[2] };
    real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; 
    
    real_t dx1[3];
    int off[3];
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
  }

  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, real_t qni_wni, real_t *vxi, dim_xyz tag_dim)
  {
    assert(0); // FIXME
  }
  
  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, real_t qni_wni, real_t *vxi)
  {
    calc_j(curr_cache, xm, xp, lf, lg, qni_wni, vxi, dim_t{});
  }
  
private:
  real_t dt_;
  Real3 dxi_;
  real_t fnqxs_, fnqys_, fnqzs_;
};

  
