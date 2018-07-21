
// ======================================================================
// 1vb current deposition "2d"
//
// This works as 2d (yz) only
// It is different in the sense the the out-of-plane current (jx)
// deposited is not the same that the 3d V-B algorithm would produce,
// but charge conservation remains exactly satisfied, anyway.

template<typename CURR_CACHE, typename dim_t>
struct Current1vb2d
{
  using curr_cache_t = CURR_CACHE;
  using real_t = typename curr_cache_t::real_t;
  
  Current1vb2d(const Grid_t& grid)
    : dt_(grid.dt),
      fnqs_(grid.norm.fnqs)
  {
    fnqys_ = grid.domain.dx[1] * grid.norm.fnqs / grid.dt;
    fnqzs_ = grid.domain.dx[2] * grid.norm.fnqs / grid.dt;
  }
  
  // ----------------------------------------------------------------------
  // calc_dx1

  static inline void
  calc_dx1(real_t dx1[2], real_t x[2], real_t dx[2], int off[2])
  {
    if (off[1] == 0) {
      dx1[0] = .5f * off[0] - x[0];
      if (dx[0] == 0.f) {
	dx1[1] = 0.f;
      } else {
	dx1[1] = dx[1] / dx[0] * dx1[0];
      }
    } else {
      dx1[1] = .5f * off[1] - x[1];
      if (dx[1] == 0.f) {
	dx1[0] = 0.f;
      } else {
	dx1[0] = dx[0] / dx[1] * dx1[1];
      }
    }
  }

  static inline void
  curr_2d_vb_cell(curr_cache_t curr_cache, int i[2], real_t x[2], real_t dx[2],
		  real_t fnq[2], real_t dxt[2], int off[2])
  {
    curr_cache.add(1, 0,i[0]  ,i[1]  , fnq[0] * dx[0] * (.5f - x[1] - .5f * dx[1]));
    curr_cache.add(1, 0,i[0]  ,i[1]+1, fnq[0] * dx[0] * (.5f + x[1] + .5f * dx[1]));
    curr_cache.add(2, 0,i[0]  ,i[1]  , fnq[1] * dx[1] * (.5f - x[0] - .5f * dx[0]));
    curr_cache.add(2, 0,i[0]+1,i[1]  , fnq[1] * dx[1] * (.5f + x[0] + .5f * dx[0]));

    if (dxt) {
      dxt[0] -= dx[0];
      dxt[1] -= dx[1];
      x[0] += dx[0] - off[0];
      x[1] += dx[1] - off[1];
      i[0] += off[0];
      i[1] += off[1];
    }
  }

  // ----------------------------------------------------------------------
  // calc_j_oop

  void calc_j_oop(curr_cache_t curr_cache, real_t qni_wni, real_t *vxi, int lf[3], real_t of[3])
  {
    real_t fnqx = vxi[0] * qni_wni * fnqs_;
    curr_cache.add(JXI, 0,lf[1]  ,lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnqx);
    curr_cache.add(JXI, 0,lf[1]+1,lf[2]  , (      of[1]) * (1.f - of[2]) * fnqx);
    curr_cache.add(JXI, 0,lf[1]  ,lf[2]+1, (1.f - of[1]) * (      of[2]) * fnqx);
    curr_cache.add(JXI, 0,lf[1]+1,lf[2]+1, (      of[1]) * (      of[2]) * fnqx);
  }

  // ----------------------------------------------------------------------
  // calc_j

  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, real_t qni_wni, real_t *vxi)
  {
    static_assert(std::is_same<dim_t, dim_yz>::value, "only dim_yz supported");
    int i[2] = { lg[1], lg[2] };
    int idiff[2] = { lf[1] - lg[1], lf[2] - lg[2] };
    real_t dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };
    real_t x[2] = { xm[1] - (i[0] + .5f), xm[2] - (i[1] + .5f) }; 

    real_t dx1[2];
    int off[2];
    int first_dir, second_dir = -1;
    /* FIXME, make sure we never div-by-zero? */
    if (idiff[0] == 0 && idiff[1] == 0) {
      first_dir = -1;
    } else if (idiff[0] == 0) {
      first_dir = 1;							
    } else if (idiff[1] == 0) {
      first_dir = 0;
    } else {
      dx1[0] = .5f * idiff[0] - x[0];
      if (dx[0] == 0.f) {
	dx1[1] = 0.f;
      } else {
	dx1[1] = dx[1] / dx[0] * dx1[0];
      }
      if (std::abs(x[1] + dx1[1]) > .5f) {
	first_dir = 1;
      } else {
	first_dir = 0;
      }
      second_dir = 1 - first_dir;
    }

    real_t fnq[2] = { qni_wni * fnqys_, qni_wni * fnqzs_ };

    if (first_dir >= 0) {
      off[1-first_dir] = 0;
      off[first_dir] = idiff[first_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(curr_cache, i, x, dx1, fnq, dx, off);
    }

    if (second_dir >= 0) {
      off[first_dir] = 0;
      off[second_dir] = idiff[second_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(curr_cache, i, x, dx1, fnq, dx, off);
    }

    curr_2d_vb_cell(curr_cache, i, x, dx, fnq, NULL, NULL);
  }

private:
  real_t dt_;
  real_t fnqs_;
  real_t fnqys_, fnqzs_;
};

