
// ======================================================================

#if DIM == DIM_YZ

CUDA_DEVICE __forceinline__ static void
calc_j2_one_cell(flds_curr_t flds_curr, particle_real_t qni_wni,
		 particle_real_t xm[3], particle_real_t xp[3])
{

  particle_real_t dx[3] = { xp[0] - xm[0], xp[1] - xm[1], xp[2] - xm[2] };
  particle_real_t xa[3]= { .5f * (xm[0] + xp[0]),
			   .5f * (xm[1] + xp[1]),
			   .5f * (xm[2] + xp[2]) };
  particle_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];

  int i[3] = {};
  for (int d = 0; d < 3; d++) {
    i[d] = particle_real_fint(xa[d]);
    xa[d] -= i[d];
  }

  particle_real_t fnqx = qni_wni * prm.fnqxs;
  curr_add(flds_curr, 0, i[0]  ,i[1]  ,i[2]  , fnqx * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
  curr_add(flds_curr, 0, i[0]  ,i[1]+1,i[2]  , fnqx * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
  curr_add(flds_curr, 0, i[0]  ,i[1]  ,i[2]+1, fnqx * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
  curr_add(flds_curr, 0, i[0]  ,i[1]+1,i[2]+1, fnqx * (dx[0] * (      xa[1]) * (      xa[2]) + h));

  particle_real_t fnqy = qni_wni * prm.fnqys;
  curr_add(flds_curr, 1, i[0]  ,i[1]  ,i[2]  , fnqy * (dx[1] * (1.f - xa[0]) * (1.f - xa[2]) + h));
  curr_add(flds_curr, 1, i[0]+0,i[1]  ,i[2]  , fnqy * (dx[1] * (      xa[0]) * (1.f - xa[2]) - h));
  curr_add(flds_curr, 1, i[0]  ,i[1]  ,i[2]+1, fnqy * (dx[1] * (1.f - xa[0]) * (      xa[2]) - h));
  curr_add(flds_curr, 1, i[0]+0,i[1]  ,i[2]+1, fnqy * (dx[1] * (      xa[0]) * (      xa[2]) + h));

  particle_real_t fnqz = qni_wni * prm.fnqzs;
  curr_add(flds_curr, 2, i[0]  ,i[1]  ,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (1.f - xa[1]) + h));
  curr_add(flds_curr, 2, i[0]+0,i[1]  ,i[2]  , fnqz * (dx[2] * (      xa[0]) * (1.f - xa[1]) - h));
  curr_add(flds_curr, 2, i[0]  ,i[1]+1,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (      xa[1]) - h));
  curr_add(flds_curr, 2, i[0]+0,i[1]+1,i[2]  , fnqz * (dx[2] * (      xa[0]) * (      xa[1]) + h));
}

#elif DIM == DIM_XYZ

CUDA_DEVICE static void
calc_j2_one_cell(flds_curr_t flds_curr, particle_real_t qni_wni,
		 particle_real_t xm[3], particle_real_t xp[3])
{
  particle_real_t dx[3] = { xp[0] - xm[0], xp[1] - xm[1], xp[2] - xm[2] };
  particle_real_t xa[3]= { .5f * (xm[0] + xp[0]),
			   .5f * (xm[1] + xp[1]),
			   .5f * (xm[2] + xp[2]) };
  particle_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];

  int i[3] = {};
  for (int d = 0; d < 3; d++) {
    i[d] = particle_real_fint(xa[d]);
    xa[d] -= i[d];
  }

  particle_real_t fnqx = qni_wni * prm.fnqxs;
  curr_add(flds_curr, 0, i[0]  ,i[1]  ,i[2]  , fnqx * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
  curr_add(flds_curr, 0, i[0]  ,i[1]+1,i[2]  , fnqx * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
  curr_add(flds_curr, 0, i[0]  ,i[1]  ,i[2]+1, fnqx * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
  curr_add(flds_curr, 0, i[0]  ,i[1]+1,i[2]+1, fnqx * (dx[0] * (      xa[1]) * (      xa[2]) + h));

  particle_real_t fnqy = qni_wni * prm.fnqys;
  curr_add(flds_curr, 1, i[0]  ,i[1]  ,i[2]  , fnqy * (dx[1] * (1.f - xa[0]) * (1.f - xa[2]) + h));
  curr_add(flds_curr, 1, i[0]+1,i[1]  ,i[2]  , fnqy * (dx[1] * (      xa[0]) * (1.f - xa[2]) - h));
  curr_add(flds_curr, 1, i[0]  ,i[1]  ,i[2]+1, fnqy * (dx[1] * (1.f - xa[0]) * (      xa[2]) - h));
  curr_add(flds_curr, 1, i[0]+1,i[1]  ,i[2]+1, fnqy * (dx[1] * (      xa[0]) * (      xa[2]) + h));

  particle_real_t fnqz = qni_wni * prm.fnqzs;
  curr_add(flds_curr, 2, i[0]  ,i[1]  ,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (1.f - xa[1]) + h));
  curr_add(flds_curr, 2, i[0]+1,i[1]  ,i[2]  , fnqz * (dx[2] * (      xa[0]) * (1.f - xa[1]) - h));
  curr_add(flds_curr, 2, i[0]  ,i[1]+1,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (      xa[1]) - h));
  curr_add(flds_curr, 2, i[0]+1,i[1]+1,i[2]  , fnqz * (dx[2] * (      xa[0]) * (      xa[1]) + h));
}

#endif

CUDA_DEVICE __forceinline__ static void
calc_j2_split_along_dim(int dim, int im, particle_real_t x1[3],
			particle_real_t xm[3], particle_real_t xp[3])
{
  particle_real_t bnd;
  if (xp[dim] > im + 1) { // crossed boundary to right
    bnd = im + 1;
  } else if (xp[dim] < im) { // crosses boundary to left
    bnd = im;
  }
  particle_real_t frac = (bnd - xm[dim]) / (xp[dim] - xm[dim]);
  // FIXME, set d == dim value to exact boundary?
  for (int d = 0; d < 3; d++) {
    if (d == dim) {
      x1[d] = bnd;
    } else {
      x1[d] = xm[d] + frac * (xp[d] - xm[d]);
    }
  }
}

#if DIM == DIM_YZ

CUDA_DEVICE __forceinline__ static inline void 
calc_j2_split_dim(flds_curr_t flds_curr, particle_real_t qni_wni,
		  particle_real_t *xm, particle_real_t *xp, int dim)
{
  if (dim == 0) {
    calc_j2_one_cell(flds_curr, qni_wni, xm, xp);
  } else {
    int im = particle_real_fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      particle_real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_split_dim(flds_curr, qni_wni, xm, x1, dim - 1);
      calc_j2_split_dim(flds_curr, qni_wni, x1, xp, dim - 1);
    } else {
      calc_j2_split_dim(flds_curr, qni_wni, xm, xp, dim - 1);
    }
  }
}

#elif DIM == DIM_XYZ

CUDA_DEVICE static inline void
calc_j2_split_dim(flds_curr_t flds_curr, particle_real_t qni_wni,
		  particle_real_t *xm, particle_real_t *xp, int dim)
{
  if (dim == -1) {
    calc_j2_one_cell(flds_curr, qni_wni, xm, xp);
  } else {
    int im = particle_real_fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      particle_real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_split_dim(flds_curr, qni_wni, xm, x1, dim - 1);
      calc_j2_split_dim(flds_curr, qni_wni, x1, xp, dim - 1);
    } else {
      calc_j2_split_dim(flds_curr, qni_wni, xm, xp, dim - 1);
    }
  }
}

#endif

// ----------------------------------------------------------------------
// calc_j

CUDA_DEVICE static inline void
calc_j(flds_curr_t flds_curr, particle_real_t *xm, particle_real_t *xp,
       int *lf, int *lg, particle_t *prt, particle_real_t *vxi)
{
  particle_real_t qni_wni = particle_qni_wni(prt);

#if DIM == DIM_YZ
  xm[0] = .5f; // this way, we guarantee that the average position will remain in the 0th cell
  xp[0] = xm[0] + vxi[0] * prm.dt * prm.dxi[0];
#endif

  calc_j2_split_dim(flds_curr, qni_wni, xm, xp, 2);
}

