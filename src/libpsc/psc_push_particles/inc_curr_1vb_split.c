
// ======================================================================

template<typename curr_cache_t>
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

  void calc_j2_one_cell(curr_cache_t curr_cache, real_t qni_wni,
			real_t xm[3], real_t xp[3], dim_yz tag_dim)
  {

    real_t dx[3] = { xp[0] - xm[0], xp[1] - xm[1], xp[2] - xm[2] };
    real_t xa[3]= { .5f * (xm[0] + xp[0]),
		    .5f * (xm[1] + xp[1]),
		    .5f * (xm[2] + xp[2]) };
    real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];

    int i[3];
    for (int d = 0; d < 3; d++) {
      i[d] = fint(xa[d]);
    }

    for (int d = 0; d < 3; d++) {
      xa[d] -= i[d];
    }

    real_t fnqx = qni_wni * fnqxs_;
    curr_cache.add(0, i[0]  ,i[1]  ,i[2]  , fnqx * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
    curr_cache.add(0, i[0]  ,i[1]+1,i[2]  , fnqx * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
    curr_cache.add(0, i[0]  ,i[1]  ,i[2]+1, fnqx * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
    curr_cache.add(0, i[0]  ,i[1]+1,i[2]+1, fnqx * (dx[0] * (      xa[1]) * (      xa[2]) + h));

    real_t fnqy = qni_wni * fnqys_;
    curr_cache.add(1, i[0]  ,i[1]  ,i[2]  , fnqy * (dx[1] * (1.f - xa[2])));
    curr_cache.add(1, i[0]  ,i[1]  ,i[2]+1, fnqy * (dx[1] * (      xa[2])));

    real_t fnqz = qni_wni * fnqzs_;
    curr_cache.add(2, i[0]  ,i[1]  ,i[2]  , fnqz * (dx[2] * (1.f - xa[1])));
    curr_cache.add(2, i[0]  ,i[1]+1,i[2]  , fnqz * (dx[2] * (      xa[1])));
  }

  void calc_j2_one_cell(curr_cache_t curr_cache, real_t qni_wni,
			real_t xm[3], real_t xp[3], dim_xyz tag_dim)
  {
    real_t dx[3] = { xp[0] - xm[0], xp[1] - xm[1], xp[2] - xm[2] };
    real_t xa[3]= { .5f * (xm[0] + xp[0]),
		    .5f * (xm[1] + xp[1]),
		    .5f * (xm[2] + xp[2]) };
    real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];

    int i[3] = {};
    for (int d = 0; d < 3; d++) {
      i[d] = fint(xa[d]);
      xa[d] -= i[d];
    }

    real_t fnqx = qni_wni * fnqxs_;
    curr_cache.add(0, i[0]  ,i[1]  ,i[2]  , fnqx * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
    curr_cache.add(0, i[0]  ,i[1]+1,i[2]  , fnqx * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
    curr_cache.add(0, i[0]  ,i[1]  ,i[2]+1, fnqx * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
    curr_cache.add(0, i[0]  ,i[1]+1,i[2]+1, fnqx * (dx[0] * (      xa[1]) * (      xa[2]) + h));

    real_t fnqy = qni_wni * fnqys_;
    curr_cache.add(1, i[0]  ,i[1]  ,i[2]  , fnqy * (dx[1] * (1.f - xa[0]) * (1.f - xa[2]) + h));
    curr_cache.add(1, i[0]+1,i[1]  ,i[2]  , fnqy * (dx[1] * (      xa[0]) * (1.f - xa[2]) - h));
    curr_cache.add(1, i[0]  ,i[1]  ,i[2]+1, fnqy * (dx[1] * (1.f - xa[0]) * (      xa[2]) - h));
    curr_cache.add(1, i[0]+1,i[1]  ,i[2]+1, fnqy * (dx[1] * (      xa[0]) * (      xa[2]) + h));

    real_t fnqz = qni_wni * fnqzs_;
    curr_cache.add(2, i[0]  ,i[1]  ,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (1.f - xa[1]) + h));
    curr_cache.add(2, i[0]+1,i[1]  ,i[2]  , fnqz * (dx[2] * (      xa[0]) * (1.f - xa[1]) - h));
    curr_cache.add(2, i[0]  ,i[1]+1,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (      xa[1]) - h));
    curr_cache.add(2, i[0]+1,i[1]+1,i[2]  , fnqz * (dx[2] * (      xa[0]) * (      xa[1]) + h));
  }

  void calc_j2_one_cell(curr_cache_t curr_cache, real_t qni_wni,
			real_t xm[3], real_t xp[3])
  {
#if DIM == DIM_YZ
    calc_j2_one_cell(curr_cache, qni_wni, xm, xp, dim_yz{});
#elif DIM == DIM_XYZ
    calc_j2_one_cell(curr_cache, qni_wni, xm, xp, dim_xyz{});
#endif
  }
  
  CUDA_DEVICE __forceinline__ static void
  calc_j2_split_along_dim(int dim, int im, real_t x1[3],
			  real_t xm[3], real_t xp[3])
  {
    real_t bnd = 0.f; // quell warning
    if (xp[dim] > im + 1) { // crossed boundary to right
      bnd = im + 1;
    } else if (xp[dim] < im) { // crosses boundary to left
      bnd = im;
    }
    real_t frac = (bnd - xm[dim]) / (xp[dim] - xm[dim]);
    // FIXME, set d == dim value to exact boundary?
    for (int d = 0; d < 3; d++) {
      if (d == dim) {
	x1[d] = bnd;
      } else {
	x1[d] = xm[d] + frac * (xp[d] - xm[d]);
      }
    }
  }

  void calc_j2_split_dim_y(curr_cache_t curr_cache, real_t qni_wni,
			   real_t *xm, real_t *xp, dim_yz tag_dim)
  {
    const int dim = 1;
    int im = fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_one_cell(curr_cache, qni_wni, xm, x1);
      calc_j2_one_cell(curr_cache, qni_wni, x1, xp);
    } else {
      calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
    }
  }

  void calc_j2_split_dim_z(curr_cache_t curr_cache, real_t qni_wni,
			   real_t *xm, real_t *xp, dim_yz tag_dim)
  {
    const int dim = 2;
    int im = fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_split_dim_y(curr_cache, qni_wni, xm, x1, tag_dim);
      calc_j2_split_dim_y(curr_cache, qni_wni, x1, xp, tag_dim);
    } else {
      calc_j2_split_dim_y(curr_cache, qni_wni, xm, xp, tag_dim);
    }
  }

  void calc_j2_split_dim_x(curr_cache_t curr_cache, real_t qni_wni,
			   real_t *xm, real_t *xp)
  {
    const int dim = 0;
    int im = fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_one_cell(curr_cache, qni_wni, xm, x1);
      calc_j2_one_cell(curr_cache, qni_wni, x1, xp);
    } else {
      calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
    }
  }

  void calc_j2_split_dim_y(curr_cache_t curr_cache, real_t qni_wni,
			   real_t *xm, real_t *xp, dim_xyz tag_dim)
  {
    const int dim = 1;
    int im = fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_split_dim_x(curr_cache, qni_wni, xm, x1);
      calc_j2_split_dim_x(curr_cache, qni_wni, x1, xp);
    } else {
      calc_j2_split_dim_x(curr_cache, qni_wni, xm, xp);
    }
  }

  void calc_j2_split_dim_z(curr_cache_t curr_cache, real_t qni_wni,
			   real_t *xm, real_t *xp, dim_xyz tag_dim)
  {
    const int dim = 2;
    int im = fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_split_dim_y(curr_cache, qni_wni, xm, x1, tag_dim);
      calc_j2_split_dim_y(curr_cache, qni_wni, x1, xp, tag_dim);
    } else {
      calc_j2_split_dim_y(curr_cache, qni_wni, xm, xp, tag_dim);
    }
  }

  void calc_j2_split_dim_z(curr_cache_t curr_cache, real_t qni_wni,
			   real_t *xm, real_t *xp)
  {
#if DIM == DIM_YZ
    calc_j2_split_dim_z(curr_cache, qni_wni, xm, xp, dim_yz{});
#elif DIM == DIM_XYZ
    calc_j2_split_dim_z(curr_cache, qni_wni, xm, xp, dim_xyz{});
#endif
  }
  
  // ----------------------------------------------------------------------
  // calc_j

  void calc_j(curr_cache_t curr_cache, real_t *xm, real_t *xp,
	      int *lf, int *lg, particle_t *prt, real_t *vxi)
  {
    real_t qni_wni = particle_qni_wni(prt);
    
#if DIM == DIM_YZ
    xm[0] = .5f; // this way, we guarantee that the average position will remain in the 0th cell
    xp[0] = xm[0] + vxi[0] * dt_ * dxi_[0];
    calc_j2_split_dim_z(curr_cache, qni_wni, xm, xp);
#else
    calc_j2_split_dim_z(curr_cache, qni_wni, xm, xp);
#endif
  }
  
private:
  real_t dt_;
  real_t fnqxs_, fnqys_, fnqzs_;
  Real3 dxi_;
};

