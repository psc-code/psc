
// ----------------------------------------------------------------------
// calc_j

__device__ static void
calc_j(const int *ci, struct d_particle *p, particles_cuda_dev_t d_particles,
       int i, real *vxi, SHAPE_INFO_ARGS, real *qni_wni, int cell_end)
{
  int j[3], k[3];
  real h0[3], h1[3];
  if (i < cell_end) {
    calc_vxi(vxi, *p);

    // x^(n+1.0), p^(n+1.0) -> x^(n+0.5), p^(n+1.0) 
    push_xi(p, vxi, -.5f * d_dt);
    find_idx_off(p->xi, j, h0, real(0.));
    
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(p, vxi, d_dt);
    find_idx_off(p->xi, k, h1, real(0.));

    push_xi(p, vxi, -.5f * d_dt); // !!! FIXME, back to before
  } else {
#if DIM == DIM_Z
    j[2] = ci[2];
    k[2] = ci[2];
#elif DIM == DIM_YZ
    j[1] = ci[1];
    k[1] = ci[1];
    j[2] = ci[2];
    k[2] = ci[2];
    h0[1] = real(0.);
    h1[1] = real(0.);
    h0[2] = real(0.);
    h1[2] = real(0.);
#endif
    vxi[0] = real(0.);
    vxi[1] = real(0.);
    vxi[2] = real(0.);
  }
#if DIM == DIM_Z  
  short int shift0z = j[2] - ci[2];
  short int shift1z = k[2] - ci[2];
  cache_shape_arrays(SHAPE_INFO_PARAMS, h0, h1, shift0z, shift1z);
#elif DIM == DIM_YZ
  short int shift0y = j[1] - ci[1];
  short int shift0z = j[2] - ci[2];
  short int shift1y = k[1] - ci[1];
  short int shift1z = k[2] - ci[2];
  cache_shape_arrays(SHAPE_INFO_PARAMS, h0, h1, shift0y, shift0z,
		     shift1y, shift1z);
#endif
  if (i < cell_end) {
    *qni_wni = p->qni_wni;
  } else {
    *qni_wni = real(0.);
  }
}

// ----------------------------------------------------------------------
// add_current_to_scratch

__device__ static void
add_current_to_scratch_z(real *scratch, int jy, int len)
{
  int tid = threadIdx.x;

  //__syncthreads(); // not necessary if len < WARPSIZE (?)
  if (tid < len) {
    int jz = tid - SW;
    scratch(0,jy,jz) += SDATA(0,jz);
  }
  //__syncthreads(); // CHECKME, necessary?
}

__device__ static void
add_current_to_scratch_y(real *scratch, int jz, int len)
{
  int tid = threadIdx.x;

  //  __syncthreads(); // not necessary if BLOCKSIZE < WARPSIZE (?)
  if (tid < len) {
    int jy = tid - SW;
    scratch(0,jy,jz) += SDATA(0,jy);
  }
  //__syncthreads(); // CHECKME, necessary?
}

// ======================================================================

#if DIM == DIM_Z

// ----------------------------------------------------------------------

__device__ static void
z_calc_jxh(real qni_wni, real vxi, SHAPE_INFO_ARGS, int jz)
{
  int tid = threadIdx.x;
  real fnqx = vxi * qni_wni * d_fnqs;

  real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
  real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z);
  real wx = .5f * (s0z + s1z);
    
  SDATA(tid,jz) = fnqx * wx;
}

__device__ static void
z_calc_jyh(real qni_wni, real vyi, SHAPE_INFO_ARGS, int jz)
{
  int tid = threadIdx.x;
  real fnqy = vyi * qni_wni * d_fnqs;

  real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
  real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z);
  real wy = .5f * (s0z + s1z);

  SDATA(tid,jz) = fnqy * wy;
}

__device__ static void
z_calc_jzh(real qni_wni, SHAPE_INFO_ARGS)
{
  int tid = threadIdx.x;
  real fnqz = qni_wni * d_fnqzs;

  real last;
  { int jz = -2;
    if (SI_SHIFT0Z >= 0 && SI_SHIFT1Z >= 0) {
      last = 0;
    } else {
      real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
      real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z);
      real wz = s1z - s0z;
      last = -fnqz*wz;
    }
    SDATA(tid,jz) = last;
  }
  for (int jz = -1; jz <= 0; jz++) {
    real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
    real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z);
    real wz = s1z - s0z;
    last -= fnqz*wz;
    SDATA(tid,jz) = last;
  }
  { int jz = 1;
    if (SI_SHIFT0Z <= 0 && SI_SHIFT1Z <= 0) {
      last = 0;
    } else {
      real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
      real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z);
      real wz = s1z - s0z;
      last -= fnqz*wz;
    }
    SDATA(tid,jz) = last;
  }
}

// ----------------------------------------------------------------------

EXTERN_C __global__ static void
push_part_p2(int n_particles, particles_cuda_dev_t d_particles, real *d_flds,
	     real *d_scratch, int block_stride, int block_start)
{
  int tid = threadIdx.x, bid = blockIdx.x * block_stride + block_start;
  int cell_begin = d_particles.offsets[bid];
  int cell_end   = d_particles.offsets[bid+1];
  int ci[3];

  blockIdx_to_blockCrd(blockIdx.x, ci);
  ci[0] *= BLOCKSIZE_X;
  ci[1] *= BLOCKSIZE_Y;
  ci[2] *= BLOCKSIZE_Z;

  int nr_loops = (cell_end - cell_begin + THREADS_PER_BLOCK-1) / THREADS_PER_BLOCK;
  for (int l = 0; l < nr_loops; l++) {
    int i = cell_begin + tid + l * THREADS_PER_BLOCK;
    DECLARE_SHAPE_INFO;
    real vxi[3], qni_wni = 0.;
    if (i < cell_end) {
      struct d_particle p;
      LOAD_PARTICLE(p, d_particles, i);
      calc_j(ci, &p, d_particles, i, vxi, SHAPE_INFO_PARAMS, &qni_wni);
    } else {
      SI_SHIFT0Z = 0;
    }

    // ----------------------------------------------------------------------
    // JX

    for (int jz = -SW; jz <= SW; jz++) {
      if (i < cell_end) {
	z_calc_jxh(qni_wni, vxi[0], SHAPE_INFO_PARAMS, jz);
      } else {
	SDATA(tid,jz) = real(0.);
      }
    }
#ifdef CALC_CURRENT
    reduce_sum_sdata(sdata);
    real *scratch = d_scratch + bid * 3 * BLOCKSTRIDE;
    add_current_to_scratch_z(scratch, 0, 5);
#endif

    // ----------------------------------------------------------------------
    // JY

    for (int jz = -SW; jz <= SW; jz++) {
      if (i < cell_end) {
	z_calc_jyh(qni_wni, vxi[1], SHAPE_INFO_PARAMS, jz);
      } else {
	SDATA(tid,jz) = real(0.);
      }
    }
#ifdef CALC_CURRENT
    reduce_sum_sdata(sdata);
    scratch = d_scratch + (bid * 3 + 1) * BLOCKSTRIDE;
    add_current_to_scratch_z(scratch, 0, 5);
#endif

    // ----------------------------------------------------------------------
    // JZ

    if (i < cell_end) {
      z_calc_jzh(qni_wni, SHAPE_INFO_PARAMS);
    } else {
      for (int jz = -SW; jz < SW; jz++) SDATA(tid,jz) = real(0.);
    }
#ifdef CALC_CURRENT
    reduce_sum_sdata4(sdata);
    scratch = d_scratch + (bid * 3 + 2) * BLOCKSTRIDE;
    add_current_to_scratch_z(scratch, 0, 4);
#endif
  }
}

#elif DIM == DIM_YZ

// ----------------------------------------------------------------------

__device__ static void
yz_calc_jxh_z(real qni_wni, real vxi, SHAPE_INFO_ARGS, int jz)
{
  int tid = threadIdx.x;
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
    
    SDATA(tid,jy) = fnqx * wx;
  }
}

__device__ static void
yz_calc_jyh_z(real qni_wni, SHAPE_INFO_ARGS, int jz)
{
  int tid = threadIdx.x;
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
    SDATA(tid,jy) = last;
  }
  for (int jy = -1; jy <= 0; jy++) {
    real s0y = pick_shape_coeff(0, y, jy, SI_SHIFT0Y);
    real s1y = pick_shape_coeff(1, y, jy, SI_SHIFT1Y) - s0y;
    real wy = s1y * tmp1;
    last -= fnqy*wy;
    SDATA(tid,jy) = last;
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
    SDATA(tid,jy) = last;
  }
}

__device__ static void
yz_calc_jzh_y(real qni_wni, SHAPE_INFO_ARGS, int jy)
{
  int tid = threadIdx.x;
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
    SDATA(tid,jz) = last;
  }
  for (int jz = -1; jz <= 0; jz++) {
    real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
    real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z) - s0z;
    real wz = s1z * tmp1;
    last -= fnqz*wz;
    SDATA(tid,jz) = last;
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
    SDATA(tid,jz) = last;
  }
}

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

__global__ static void
push_part_p2(int n_particles, particles_cuda_dev_t d_particles, real *d_flds,
	     real *d_scratch, int block_stride, int block_start)
{
  int tid = threadIdx.x, bid = blockIdx.x * block_stride + block_start;

  int stride = (BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW);
  // stride must be <= THREADS_PER_BLOCK!
  if (tid < stride) {
    for (int m = 0; m < 3; m++) {
      scurr[m * stride + tid] = real(0.);
    }
  }
  __syncthreads();

  int cell_begin = d_particles.offsets[bid];
  int cell_end   = d_particles.offsets[bid+1];
  int ci[3];

  blockIdx_to_blockCrd(bid, ci);
  ci[0] *= BLOCKSIZE_X;
  ci[1] *= BLOCKSIZE_Y;
  ci[2] *= BLOCKSIZE_Z;

  int nr_loops = (cell_end - cell_begin + THREADS_PER_BLOCK-1) / THREADS_PER_BLOCK;
  for (int l = 0; l < nr_loops; l++) {
    int i = cell_begin + tid + l * THREADS_PER_BLOCK;
    DECLARE_SHAPE_INFO;
    real vxi[3], qni_wni;
    struct d_particle p;
    if (i < cell_end) {
      LOAD_PARTICLE(p, d_particles, i);
    }
    calc_j(ci, &p, d_particles, i, vxi, SHAPE_INFO_PARAMS, &qni_wni, cell_end);
    yz_calc_jx(bid, d_scratch, vxi[0], qni_wni, SHAPE_INFO_PARAMS);
    yz_calc_jy(bid, d_scratch, qni_wni, SHAPE_INFO_PARAMS);
    yz_calc_jz(bid, d_scratch, qni_wni, SHAPE_INFO_PARAMS);
  }

  __syncthreads();
  if (tid < stride) {
    for (int m = 0; m < 3; m++) {
      real *scratch = d_scratch + (bid * 3 + m) * BLOCKSTRIDE;
      int jz = tid / (BLOCKSIZE_Y + 2*SW) - SW;
      int jy = tid % (BLOCKSIZE_Y + 2*SW) - SW;
      scratch(0,jy,jz) += scurr[m*stride + tid];
    }
  }
}

#endif

