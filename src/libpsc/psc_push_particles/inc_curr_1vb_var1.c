
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

CUDA_DEVICE static void
curr_3d_vb_cell(curr_cache_t curr_cache, int i[3], particle_real_t x[3], particle_real_t dx[3],
		particle_real_t qni_wni)
{
  real xa[3] = { 0.,
		 x[1] + .5f * dx[1],
		 x[2] + .5f * dx[2], };
#ifdef __CUDACC__
  if (dx[0] != 0.f)
#endif
    {
      real fnqx = qni_wni * prm.fnqxs;
      real h = (1.f / 12.f) * dx[0] * dx[1] * dx[2];
      curr_add(curr_cache, 0, 0,i[1]  ,i[2]  , fnqx * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h));
      curr_add(curr_cache, 0, 0,i[1]+1,i[2]  , fnqx * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h));
      curr_add(curr_cache, 0, 0,i[1]  ,i[2]+1, fnqx * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) - h));
      curr_add(curr_cache, 0, 0,i[1]+1,i[2]+1, fnqx * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) + h));
    }
#ifdef __CUDACC__
  if (dx[1] != 0.f)
#endif
    {
      real fnqy = qni_wni * prm.fnqys;
      curr_add(curr_cache, 1, 0,i[1],i[2]  , fnqy * dx[1] * (.5f - xa[2]));
      curr_add(curr_cache, 1, 0,i[1],i[2]+1, fnqy * dx[1] * (.5f + xa[2]));
    }
#ifdef __CUDACC__
  if (dx[2] != 0.f)
#endif
    {
      real fnqz = qni_wni * prm.fnqzs;
      curr_add(curr_cache, 2, 0,i[1]  ,i[2], fnqz * dx[2] * (.5f - xa[1]));
      curr_add(curr_cache, 2, 0,i[1]+1,i[2], fnqz * dx[2] * (.5f + xa[1]));
    }
}

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

CUDA_DEVICE static void
calc_j(curr_cache_t curr_cache, particle_real_t *xm, particle_real_t *xp,
       int *lf, int *lg, particle_t *prt, particle_real_t *vxi)

{
  // deposit xm -> xp
  int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			
  int i[3] = { 0, lg[1], lg[2] };					
  particle_real_t dx[3] = { vxi[0] * prm.dt * prm.dxi[0], xp[1] - xm[1], xp[2] - xm[2] };
  particle_real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; 

  particle_real_t dx1[3];
  int off[3];
#ifdef __CUDACC__

  real x1 = x[1] * idiff[1];
  real x2 = x[2] * idiff[2];
  int d_first = (fabsf(dx[2]) * (.5f - x1) >= fabsf(dx[1]) * (.5f - x2));

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
    if (particle_real_abs(x[2] + dx1[2]) > .5f) {
      first_dir = 2;
    } else {
      first_dir = 1;
    }
    second_dir = 3 - first_dir;
  }
  particle_real_t qni_wni = particle_qni_wni(prt);

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


#elif DIM == DIM_XYZ

#ifdef __CUDACC__

CUDA_DEVICE static void
calc_j(curr_cache_t curr_cache, particle_real_t *xm, particle_real_t *xp,
       int *lf, int *lg, particle_t *prt, particle_real_t *vxi)
{
  assert(0);
}

#endif

#endif // DIM

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

