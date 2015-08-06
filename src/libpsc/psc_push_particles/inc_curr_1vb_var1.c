
// ======================================================================
// 1vb current deposition "var1"
//
// my original implementation, 2d (yz) only
// in theory less divergent, in particular for CUDA,
// but also definitely more complex

// ----------------------------------------------------------------------
// calc_3d_dx1

#if DIM == DIM_YZ

CUDA_DEVICE static inline void
calc_3d_dx1(particle_real_t dx1[3], particle_real_t x[3], particle_real_t dx[3], int off[3])
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

#ifdef __CUDACC__

// ----------------------------------------------------------------------
// curr_3d_vb_cell

CUDA_DEVICE static void
curr_3d_vb_cell(real *d_flds, int i[3], real x[3], real dx[3], real qni_wni,
		int *ci0)
{
  real xa[3] = { 0.,
		 x[1] + .5f * dx[1],
		 x[2] + .5f * dx[2], };
  if (dx[0] != 0.f) {
    real fnqx = qni_wni * prm.fnqxs;
    real h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
    curr_add(d_flds, 0, 0,i[1]  ,i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h), ci0);
    curr_add(d_flds, 0, 0,i[1]+1,i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h), ci0);
    curr_add(d_flds, 0, 0,i[1]  ,i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h), ci0);
    curr_add(d_flds, 0, 0,i[1]+1,i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h), ci0);
  }
  if (dx[1] != 0.f) {
    real fnqy = qni_wni * prm.fnqys;
    curr_add(d_flds, 1, 0,i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]), ci0);
    curr_add(d_flds, 1, 0,i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]), ci0);
  }
  if (dx[2] != 0.f) {
    real fnqz = qni_wni * prm.fnqzs;
    curr_add(d_flds, 2, 0,i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]), ci0);
    curr_add(d_flds, 2, 0,i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]), ci0);
  }
}

#else

CUDA_DEVICE static inline void
curr_3d_vb_cell(struct psc_fields *flds, int i[3], particle_real_t x[3], particle_real_t dx[3],
		particle_real_t fnq[3])
{
  particle_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];
  particle_real_t xa[3] = { 0.,
			    x[1] + .5f * dx[1],
			    x[2] + .5f * dx[2], };
  F3_CURR(flds, JXI, 0,i[1]  ,i[2]  ) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h);
  F3_CURR(flds, JXI, 0,i[1]+1,i[2]  ) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h);
  F3_CURR(flds, JXI, 0,i[1]  ,i[2]+1) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h);
  F3_CURR(flds, JXI, 0,i[1]+1,i[2]+1) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h);

  F3_CURR(flds, JYI, 0,i[1]  ,i[2]  ) += fnq[1] * dx[1] * (.5f - xa[2]);
  F3_CURR(flds, JYI, 0,i[1]  ,i[2]+1) += fnq[1] * dx[1] * (.5f + xa[2]);

  F3_CURR(flds, JZI, 0,i[1]  ,i[2]  ) += fnq[2] * dx[2] * (.5f - xa[1]);
  F3_CURR(flds, JZI, 0,i[1]+1,i[2]  ) += fnq[2] * dx[2] * (.5f + xa[1]);
}

#endif

// ----------------------------------------------------------------------
// curr_3d_vb_cell_upd

CUDA_DEVICE static void
curr_3d_vb_cell_upd(int i[3], particle_real_t x[3], particle_real_t dx1[3],
		    particle_real_t dx[3], int off[3])
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

#ifdef __CUDACC__

CUDA_DEVICE static void
calc_j(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
       real *d_flds, int p_nr, int bid, int *ci0)
{
  real vxi[3];
  calc_vxi(vxi, prt);

  // position xm at x^(n+.5)
  real h0[3], h1[3];
  real xm[3], xp[3];
  int j[3], k[3];
  
  find_idx_off_pos_1st_rel(prt->xi, j, h0, xm, real(0.));

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, prm.dt);
  STORE_PARTICLE_POS_(*prt, d_xi4, n);

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

  // position xp at x^(n+.5)
  find_idx_off_pos_1st_rel(prt->xi, k, h1, xp, real(0.));

  // deposit xm -> xp
  int idiff[3] = { 0, k[1] - j[1], k[2] - j[2] };
  int i[3] = { 0, j[1], j[2] };
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
  calc_3d_dx1(dx1, x, dx, off);
  curr_3d_vb_cell(d_flds, i, x, dx1, prt->qni_wni, ci0);
  curr_3d_vb_cell_upd(i, x, dx1, dx, off);
  
  off[1] = idiff[1] - off[1];
  off[2] = idiff[2] - off[2];
  calc_3d_dx1(dx1, x, dx, off);
  curr_3d_vb_cell(d_flds, i, x, dx1, prt->qni_wni, ci0);
  curr_3d_vb_cell_upd(i, x, dx1, dx, off);
    
  curr_3d_vb_cell(d_flds, i, x, dx, prt->qni_wni, ci0);
}

#else

CUDA_DEVICE static inline void
calc_j(struct psc_fields *flds, particle_real_t *xm, particle_real_t *xp,
       int *lf, int *lg, particle_t *prt, particle_real_t *vxi)
{
  int i[3] = { 0, lg[1], lg[2] };					
  int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			
  particle_real_t dx[3] = { 0., xp[1] - xm[1], xp[2] - xm[2] };		
  particle_real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; 
  									
  particle_real_t dx1[3];
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
    if (particle_real_abs(x[2] + dx1[2]) > .5f) {
      first_dir = 2;
    } else {
      first_dir = 1;
    }
    second_dir = 3 - first_dir;
  }
  particle_real_t fnq[3] = { particle_qni_wni(prt) * prm.fnqxs,
			     particle_qni_wni(prt) * prm.fnqys,
			     particle_qni_wni(prt) * prm.fnqzs };
  dx[0] = vxi[0] * prm.dt * prm.dxi[0];
  if (first_dir >= 0) {
    off[3 - first_dir] = 0;
    off[first_dir] = idiff[first_dir];
    calc_3d_dx1(dx1, x, dx, off);
    curr_3d_vb_cell(flds, i, x, dx1, fnq);
    curr_3d_vb_cell_upd(i, x, dx1, dx, off);
  }

  if (second_dir >= 0) {
    off[first_dir] = 0;
    off[second_dir] = idiff[second_dir];
    calc_3d_dx1(dx1, x, dx, off);
    curr_3d_vb_cell(flds, i, x, dx1, fnq);
    curr_3d_vb_cell_upd(i, x, dx1, dx, off);
  }

  curr_3d_vb_cell(flds, i, x, dx, fnq);
}

#endif

#elif DIM == DIM_XYZ

#ifdef __CUDACC__

CUDA_DEVICE static void
calc_j(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
       real *d_flds, int p_nr, int bid, int *ci0)
{
  real vxi[3];
  calc_vxi(vxi, prt);

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, prm.dt);
  STORE_PARTICLE_POS_(*prt, d_xi4, n);

  assert(0);
}

#endif

#endif // DIM

