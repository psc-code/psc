
// ======================================================================

#if DIM == DIM_YZ

CUDA_DEVICE __forceinline__ static void
calc_j2_one_cell(curr_cache_t curr_cache, particle_real_t qni_wni,
		 particle_real_t xm[3], particle_real_t xp[3])
{

  particle_real_t dx[3] = { xp[0] - xm[0], xp[1] - xm[1], xp[2] - xm[2] };
  particle_real_t xa[3]= { .5f * (xm[0] + xp[0]),
			   .5f * (xm[1] + xp[1]),
			   .5f * (xm[2] + xp[2]) };
  particle_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];

  int i[3];
  for (int d = 0; d < 3; d++) {
    i[d] = particle_real_fint(xa[d]);
  }

#ifdef CURR_CACHE_HAVE_SHIFT
  curr_cache = curr_cache_shift(curr_cache, 0, i[0], i[1], i[2]);
#endif

  for (int d = 0; d < 3; d++) {
    xa[d] -= i[d];
  }

  particle_real_t fnqx = qni_wni * prm.fnqxs;
#ifdef CURR_CACHE_HAVE_SHIFT
  curr_cache_add(curr_cache, 0, 0,0,0, fnqx * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
  curr_cache_add(curr_cache, 0, 0,1,0, fnqx * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
  curr_cache_add(curr_cache, 0, 0,0,1, fnqx * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
  curr_cache_add(curr_cache, 0, 0,1,1, fnqx * (dx[0] * (      xa[1]) * (      xa[2]) + h));
#else
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]  ,i[2]  , fnqx * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]+1,i[2]  , fnqx * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]  ,i[2]+1, fnqx * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]+1,i[2]+1, fnqx * (dx[0] * (      xa[1]) * (      xa[2]) + h));
#endif

  particle_real_t fnqy = qni_wni * prm.fnqys;
#ifdef CURR_CACHE_HAVE_SHIFT
  curr_cache_add(curr_cache, 1, 0,0,0, fnqy * (dx[1] * (1.f - xa[2])));
  curr_cache_add(curr_cache, 1, 0,0,1, fnqy * (dx[1] * (      xa[2])));
#else
  curr_cache_add(curr_cache, 1, i[0]  ,i[1]  ,i[2]  , fnqy * (dx[1] * (1.f - xa[2])));
  curr_cache_add(curr_cache, 1, i[0]  ,i[1]  ,i[2]+1, fnqy * (dx[1] * (      xa[2])));
#endif

  particle_real_t fnqz = qni_wni * prm.fnqzs;
#ifdef CURR_CACHE_HAVE_SHIFT
  curr_cache_add(curr_cache, 2, 0,0,0, fnqz * (dx[2] * (1.f - xa[1])));
  curr_cache_add(curr_cache, 2, 0,1,0, fnqz * (dx[2] * (      xa[1])));
#else
  curr_cache_add(curr_cache, 2, i[0]  ,i[1]  ,i[2]  , fnqz * (dx[2] * (1.f - xa[1])));
  curr_cache_add(curr_cache, 2, i[0]  ,i[1]+1,i[2]  , fnqz * (dx[2] * (      xa[1])));
#endif
}

#elif DIM == DIM_XYZ

CUDA_DEVICE __forceinline__ static void
calc_j2_one_cell(curr_cache_t curr_cache, particle_real_t qni_wni,
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
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]  ,i[2]  , fnqx * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]+1,i[2]  , fnqx * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]  ,i[2]+1, fnqx * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
  curr_cache_add(curr_cache, 0, i[0]  ,i[1]+1,i[2]+1, fnqx * (dx[0] * (      xa[1]) * (      xa[2]) + h));

  particle_real_t fnqy = qni_wni * prm.fnqys;
  curr_cache_add(curr_cache, 1, i[0]  ,i[1]  ,i[2]  , fnqy * (dx[1] * (1.f - xa[0]) * (1.f - xa[2]) + h));
  curr_cache_add(curr_cache, 1, i[0]+1,i[1]  ,i[2]  , fnqy * (dx[1] * (      xa[0]) * (1.f - xa[2]) - h));
  curr_cache_add(curr_cache, 1, i[0]  ,i[1]  ,i[2]+1, fnqy * (dx[1] * (1.f - xa[0]) * (      xa[2]) - h));
  curr_cache_add(curr_cache, 1, i[0]+1,i[1]  ,i[2]+1, fnqy * (dx[1] * (      xa[0]) * (      xa[2]) + h));

  particle_real_t fnqz = qni_wni * prm.fnqzs;
  curr_cache_add(curr_cache, 2, i[0]  ,i[1]  ,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (1.f - xa[1]) + h));
  curr_cache_add(curr_cache, 2, i[0]+1,i[1]  ,i[2]  , fnqz * (dx[2] * (      xa[0]) * (1.f - xa[1]) - h));
  curr_cache_add(curr_cache, 2, i[0]  ,i[1]+1,i[2]  , fnqz * (dx[2] * (1.f - xa[0]) * (      xa[1]) - h));
  curr_cache_add(curr_cache, 2, i[0]+1,i[1]+1,i[2]  , fnqz * (dx[2] * (      xa[0]) * (      xa[1]) + h));
}

#endif

CUDA_DEVICE __forceinline__ static void
calc_j2_split_along_dim(int dim, int im, particle_real_t x1[3],
			particle_real_t xm[3], particle_real_t xp[3])
{
  particle_real_t bnd = 0.f; // quell warning
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

#if 0
CUDA_DEVICE __forceinline__ static void 
calc_j2_split_dim(curr_cache_t curr_cache, particle_real_t qni_wni,
		  particle_real_t *xm, particle_real_t *xp, int dim)
{
  if (dim == 0) {
    calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
  } else {
    int im = particle_real_fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      particle_real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_split_dim(curr_cache, qni_wni, xm, x1, dim - 1);
      calc_j2_split_dim(curr_cache, qni_wni, x1, xp, dim - 1);
    } else {
      calc_j2_split_dim(curr_cache, qni_wni, xm, xp, dim - 1);
    }
  }
}
#endif

CUDA_DEVICE __forceinline__ static void
calc_j2_split_dim_y(curr_cache_t curr_cache, particle_real_t qni_wni,
		    particle_real_t *xm, particle_real_t *xp)
{
  const int dim = 1;
  int im = particle_real_fint(xm[dim]);
  if (xp[dim] > im + 1 || xp[dim] < im) {
    particle_real_t x1[3];
    calc_j2_split_along_dim(dim, im, x1, xm, xp);
    calc_j2_one_cell(curr_cache, qni_wni, xm, x1);
    calc_j2_one_cell(curr_cache, qni_wni, x1, xp);
  } else {
    calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
  }
}

CUDA_DEVICE __forceinline__ static void
calc_j2_split_dim_z(curr_cache_t curr_cache, particle_real_t qni_wni,
		    particle_real_t *xm, particle_real_t *xp)
{
  const int dim = 2;
  int im = particle_real_fint(xm[dim]);
  if (xp[dim] > im + 1 || xp[dim] < im) {
    particle_real_t x1[3];
    calc_j2_split_along_dim(dim, im, x1, xm, xp);
    calc_j2_split_dim_y(curr_cache, qni_wni, xm, x1);
    calc_j2_split_dim_y(curr_cache, qni_wni, x1, xp);
  } else {
    calc_j2_split_dim_y(curr_cache, qni_wni, xm, xp);
  }
}

#elif DIM == DIM_XYZ

#if 0
CUDA_DEVICE __forceinline__ static void
calc_j2_split_dim(curr_cache_t curr_cache, particle_real_t qni_wni,
		  particle_real_t *xm, particle_real_t *xp, int dim)
{
  if (dim == -1) {
    calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
  } else {
    int im = particle_real_fint(xm[dim]);
    if (xp[dim] > im + 1 || xp[dim] < im) {
      particle_real_t x1[3];
      calc_j2_split_along_dim(dim, im, x1, xm, xp);
      calc_j2_split_dim(curr_cache, qni_wni, xm, x1, dim - 1);
      calc_j2_split_dim(curr_cache, qni_wni, x1, xp, dim - 1);
    } else {
      calc_j2_split_dim(curr_cache, qni_wni, xm, xp, dim - 1);
    }
  }
}
#endif

CUDA_DEVICE __forceinline__ static void
calc_j2_split_dim_x(curr_cache_t curr_cache, particle_real_t qni_wni,
		    particle_real_t *xm, particle_real_t *xp)
{
  const int dim = 0;
  int im = particle_real_fint(xm[dim]);
  if (xp[dim] > im + 1 || xp[dim] < im) {
    particle_real_t x1[3];
    calc_j2_split_along_dim(dim, im, x1, xm, xp);
    calc_j2_one_cell(curr_cache, qni_wni, xm, x1);
    calc_j2_one_cell(curr_cache, qni_wni, x1, xp);
  } else {
    calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
  }
}

CUDA_DEVICE __forceinline__ static void
calc_j2_split_dim_y(curr_cache_t curr_cache, particle_real_t qni_wni,
		    particle_real_t *xm, particle_real_t *xp)
{
  const int dim = 1;
  int im = particle_real_fint(xm[dim]);
  if (xp[dim] > im + 1 || xp[dim] < im) {
    particle_real_t x1[3];
    calc_j2_split_along_dim(dim, im, x1, xm, xp);
    calc_j2_split_dim_x(curr_cache, qni_wni, xm, x1);
    calc_j2_split_dim_x(curr_cache, qni_wni, x1, xp);
  } else {
    calc_j2_split_dim_x(curr_cache, qni_wni, xm, xp);
  }
}

CUDA_DEVICE __forceinline__ static void
calc_j2_split_dim_z(curr_cache_t curr_cache, particle_real_t qni_wni,
		    particle_real_t *xm, particle_real_t *xp)
{
  const int dim = 2;
  int im = particle_real_fint(xm[dim]);
  if (xp[dim] > im + 1 || xp[dim] < im) {
    particle_real_t x1[3];
    calc_j2_split_along_dim(dim, im, x1, xm, xp);
    calc_j2_split_dim_y(curr_cache, qni_wni, xm, x1);
    calc_j2_split_dim_y(curr_cache, qni_wni, x1, xp);
  } else {
    calc_j2_split_dim_y(curr_cache, qni_wni, xm, xp);
  }
}

#endif

// ----------------------------------------------------------------------
// calc_j

CUDA_DEVICE __forceinline__ static void
calc_j(curr_cache_t curr_cache, particle_real_t *xm, particle_real_t *xp,
       int *lf, int *lg, particle_t *prt, particle_real_t *vxi)
{
  particle_real_t qni_wni = particle_qni_wni(prt);

#if DIM == DIM_YZ
  xm[0] = .5f; // this way, we guarantee that the average position will remain in the 0th cell
  xp[0] = xm[0] + vxi[0] * prm.dt * prm.dxi[0];
  calc_j2_split_dim_z(curr_cache, qni_wni, xm, xp);
#else
  calc_j2_split_dim_z(curr_cache, qni_wni, xm, xp);
#endif
}

