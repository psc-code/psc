
// ----------------------------------------------------------------------
// calc_j

__device__ static void
calc_j(const int *ci, int i, particles_cuda_dev_t d_particles,
       real *vxi, SHAPE_INFO_ARGS, real *qni_wni, int cell_end)
{
#if DIM == DIM_Z  
  short int shift0z;
  short int shift1z;
#elif DIM == DIM_YZ
  short int shift0y;
  short int shift0z;
  short int shift1y;
  short int shift1z;
#endif
  real h0[3], h1[3];
  if (i < cell_end) {
    struct d_particle p;
    LOAD_PARTICLE(p, d_particles, i);
    *qni_wni = p.qni_wni;
    int j[3], k[3];

    calc_vxi(vxi, p);

    // x^(n+1.0), p^(n+1.0) -> x^(n+0.5), p^(n+1.0) 
    push_xi(&p, vxi, -.5f * d_dt);
    find_idx_off(p.xi, j, h0, real(0.));
    
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(&p, vxi, d_dt);
    find_idx_off(p.xi, k, h1, real(0.));

#if DIM == DIM_Z  
    shift0z = j[2] - ci[2];
    shift1z = k[2] - ci[2];
#elif DIM == DIM_YZ
    shift0y = j[1] - ci[1];
    shift0z = j[2] - ci[2];
    shift1y = k[1] - ci[1];
    shift1z = k[2] - ci[2];
#endif
  } else {
    *qni_wni = real(0.);
    vxi[0] = real(0.);
    vxi[1] = real(0.);
    vxi[2] = real(0.);
#if DIM == DIM_Z  
    shift0z = 0;
    shift1z = 0;
    h0[2] = real(0.);
    h1[2] = real(0.);
#elif DIM == DIM_YZ
    shift0y = 0;
    shift0z = 0;
    shift1y = 0;
    shift1z = 0;
    h0[1] = real(0.);
    h1[1] = real(0.);
    h0[2] = real(0.);
    h1[2] = real(0.);
#endif
  }
#if DIM == DIM_Z
  cache_shape_arrays(SHAPE_INFO_PARAMS, h0, h1, shift0z, shift1z);
#elif DIM == DIM_YZ
  cache_shape_arrays(SHAPE_INFO_PARAMS, h0, h1, shift0y, shift0z,
		     shift1y, shift1z);
#endif
}

// ======================================================================

#if DIM == DIM_YZ

// ----------------------------------------------------------------------

__shared__ real scurr[(BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW) * 3];

#define scurr(m,jy,jz) (scurr[(m * (BLOCKSIZE_Z + 2*SW) + (jz)+SW)	\
			      * (BLOCKSIZE_Y + 2*SW) + (jy)+SW])

// ----------------------------------------------------------------------
// yz_calc_jx

__device__ static void
yz_calc_jx(int bid, real *d_scratch, real vxi, real qni_wni, SHAPE_INFO_ARGS)
{
  int tid = threadIdx.x;

  for (int jz = -SW; jz <= SW; jz++) {
    real fnqx = vxi * qni_wni * d_fnqs;
    
    // FIXME, can be simplified
    real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
    real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z) - s0z;
    
    for (int jy = -SW; jy <= SW; jy++) {
      real s0y = pick_shape_coeff(0, y, jy, SI_SHIFT0Y);
      real s1y = pick_shape_coeff(1, y, jy, SI_SHIFT1Y) - s0y;
      real wx = s0y * s0z
	+ real(.5) * s1y * s0z
	+ real(.5) * s0y * s1z
	+ real(.3333333333) * s1y * s1z;
      
      reduce_sum(fnqx * wx);
#ifdef CALC_CURRENT
      if (tid == 0) {
	scurr(0,jy,jz) += sdata1[0];
      }
#endif
    }
  }
}

// ----------------------------------------------------------------------
// yz_calc_jy

__device__ static void
yz_calc_jy(int bid, real *d_scratch, real qni_wni, SHAPE_INFO_ARGS)
{
  int tid = threadIdx.x;
  
  for (int jz = -SW; jz <= SW; jz++) {
    real fnqy = qni_wni * d_fnqys;
    
    real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
    real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z);
    real tmp1 = real(.5) * (s0z + s1z);
    
    real last;
    { int jy = -2;
      if (SI_SHIFT0Y >= 0 && SI_SHIFT1Y >= 0) {
	last = 0.f;
      } else {
	real s0y = pick_shape_coeff(0, y, jy, SI_SHIFT0Y);
	real s1y = pick_shape_coeff(1, y, jy, SI_SHIFT1Y) - s0y;
	real wy = s1y * tmp1;
	last = -fnqy*wy;
      }
      reduce_sum(last);
#ifdef CALC_CURRENT
      if (tid == 0) {
	scurr(1,jy,jz) += sdata1[0];
      }
#endif
    }
    for (int jy = -1; jy <= 0; jy++) {
      real s0y = pick_shape_coeff(0, y, jy, SI_SHIFT0Y);
      real s1y = pick_shape_coeff(1, y, jy, SI_SHIFT1Y) - s0y;
      real wy = s1y * tmp1;
      last -= fnqy*wy;
      reduce_sum(last);
#ifdef CALC_CURRENT
      if (tid == 0) {
	scurr(1,jy,jz) += sdata1[0];
      }
#endif
    }
    { int jy = 1;
      if (SI_SHIFT0Y <= 0 && SI_SHIFT1Y <= 0) {
	last = 0.f;
      } else {
	real s0y = pick_shape_coeff(0, y, jy, SI_SHIFT0Y);
	real s1y = pick_shape_coeff(1, y, jy, SI_SHIFT1Y) - s0y;
	real wy = s1y * tmp1;
	last -= fnqy*wy;
      }
      reduce_sum(last);
#ifdef CALC_CURRENT
      if (tid == 0) {
	scurr(1,jy,jz) += sdata1[0];
      }
#endif
    }
  }
}

// ----------------------------------------------------------------------
// yz_calc_jz

__device__ static void
yz_calc_jz(int bid, real *d_scratch, real qni_wni, SHAPE_INFO_ARGS)
{
  int tid = threadIdx.x;
  
  for (int jy = -SW; jy <= SW; jy++) {
    real fnqz = qni_wni * d_fnqzs;
    
    real s0y = pick_shape_coeff(0, y, jy, SI_SHIFT0Y);
    real s1y = pick_shape_coeff(1, y, jy, SI_SHIFT1Y);
    real tmp1 = real(.5) * (s0y + s1y);
    
    real last;
    { int jz = -2;
      if (SI_SHIFT0Z >= 0 && SI_SHIFT1Z >= 0) {
	last = 0.f;
      } else {
	real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
	real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z) - s0z;
	real wz = s1z * tmp1;
	last = -fnqz*wz;
      }
      reduce_sum(last);
#ifdef CALC_CURRENT
      if (tid == 0) {
	scurr(2,jy,jz) += sdata1[0];
      }
#endif
    }
    for (int jz = -1; jz <= 0; jz++) {
      real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
      real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z) - s0z;
      real wz = s1z * tmp1;
      last -= fnqz*wz;
      reduce_sum(last);
#ifdef CALC_CURRENT
      if (tid == 0) {
	scurr(2,jy,jz) += sdata1[0];
      }
#endif
    }
    { int jz = 1;
      if (SI_SHIFT0Z <= 0 && SI_SHIFT1Z <= 0) {
	last = 0.f;
      } else {
	real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
	real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z) - s0z;
	real wz = s1z * tmp1;
	last -= fnqz*wz;
      }
      reduce_sum(last);
#ifdef CALC_CURRENT
      if (tid == 0) {
	scurr(2,jy,jz) += sdata1[0];
      }
#endif
    }
  }
}

__device__ static void
zero_scurr()
{
  int i = threadIdx.x;
  int N = (BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW) * 3;
  while (i < N) {
    scurr[i] = real(0.);
    i += THREADS_PER_BLOCK;
  }
}

__device__ static void
add_scurr_to_scratch(real *d_scratch, int bid)
{
  int i = threadIdx.x;
  int stride = (BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW) * 3;
  while (i < stride) {
    int rem = i;
    int m = rem / ((BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Y + 2*SW));
    rem -= m * ((BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Y + 2*SW));
    int jz = rem / (BLOCKSIZE_Y + 2*SW);
    rem -= jz * (BLOCKSIZE_Y + 2*SW);
    int jy = rem;
    jz -= SW;
    jy -= SW;
    real *scratch = d_scratch + (bid * 3 + m) * BLOCKSTRIDE;
    scratch(0,jy,jz) += scurr(m, jy, jz);
    i += THREADS_PER_BLOCK;
  }
}

__global__ static void
push_part_p2(int n_particles, particles_cuda_dev_t d_particles, real *d_flds,
	     real *d_scratch, int block_stride, int block_start)
{
  int tid = threadIdx.x, bid = blockIdx.x * block_stride + block_start;

  zero_scurr();
  __syncthreads();

  int cell_begin = d_particles.offsets[bid];
  int cell_end   = d_particles.offsets[bid+1];
  int ci[3];

  blockIdx_to_cellPos(&d_particles, bid, ci);

  int nr_loops = (cell_end - cell_begin + THREADS_PER_BLOCK-1) / THREADS_PER_BLOCK;
  for (int l = 0; l < nr_loops; l++) {
    int i = cell_begin + tid + l * THREADS_PER_BLOCK;
    DECLARE_SHAPE_INFO;
    real vxi[3], qni_wni;
    calc_j(ci, i, d_particles, vxi, SHAPE_INFO_PARAMS, &qni_wni, cell_end);
    yz_calc_jx(bid, d_scratch, vxi[0], qni_wni, SHAPE_INFO_PARAMS);
    yz_calc_jy(bid, d_scratch, qni_wni, SHAPE_INFO_PARAMS);
    yz_calc_jz(bid, d_scratch, qni_wni, SHAPE_INFO_PARAMS);
  }

  __syncthreads();
  add_scurr_to_scratch(d_scratch, bid);
}

#endif

