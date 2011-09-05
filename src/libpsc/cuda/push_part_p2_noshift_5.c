
struct shapeinfo {
  real hy[2], hz[2];
  short int shift0[2], shift1[2];
};

#define DECLARE_SHAPE_INFO			\
  struct shapeinfo _si, *si = &_si

#define SHAPE_INFO_PARAMS si
#define SHAPE_INFO_ARGS struct shapeinfo *si

#define SI_SHIFT0Y si->shift0[0]
#define SI_SHIFT1Y si->shift1[0]
#define SI_SHIFT10Y (si->shift1[0] - si->shift0[0])
#define SI_SHIFT0Z si->shift0[1]
#define SI_SHIFT1Z si->shift1[1]
#define SI_SHIFT10Z (si->shift1[1] - si->shift0[1])

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *h0, real *h1,
		   short int shift0y, short int shift0z,
		   short int shift1y, short int shift1z)
{
  si->hy[0] = h0[1];
  si->hy[1] = h1[1];
  si->hz[0] = h0[2];
  si->hz[1] = h1[2];
  si->shift0[0] = shift0y;
  si->shift0[1] = shift0z;
  si->shift1[0] = shift1y;
  si->shift1[1] = shift1z;
}

#define pick_shape_coeff(t, comp, j, shift) ({				\
      const int __y __attribute__((unused)) = 1;			\
      const int __z __attribute__((unused)) = 2;			\
      __pick_shape_coeff(j, shift, __ ## comp, si->h ## comp[t]);		\
  })


__device__ static real
__pick_shape_coeff(int j, int shift, int d, real h)
{
  real s;
  if (j == shift - 1) {
    s = find_shape_coeff_d_shift(-1, h, 0);
  } else if (j == shift + 0) {
    s = find_shape_coeff_d_shift( 0, h, 0);
  } else if (j == shift + 1) {
    s = find_shape_coeff_d_shift(+1, h, 0);
  } else {
    s = real(0.);
  }
  return s;
}

// based on p2_noshift_4.c

#if DIM != DIM_YZ
#error TBD
#endif

__shared__ int ci0[3]; // cell index of lower-left cell in block

// ======================================================================

// max blocksize_[xyz] = 256!

__device__ static int
encode_ci1(int ci1[3])
{
  return ci1[0] | (ci1[1] << 8) | (ci1[2] << 16);
}

__device__ static void
decode_ci1(int val, int ci1[3])
{
  ci1[0] = val & 0xff;
  ci1[1] = (val >> 8) & 0xff;
  ci1[2] = (val >> 16) & 0xff;
}

// ----------------------------------------------------------------------
// calc_shape_info

__device__ static void
calc_shape_info(int *ci1, int i, particles_cuda_dev_t d_particles,
		real *vxi, real *qni_wni, SHAPE_INFO_ARGS, int cell_end)
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
    find_idx(p.xi, ci1, real(0.));
    ci1[1] -= ci0[1];
    ci1[2] -= ci0[2];

    int j[3], k[3];
    calc_vxi(vxi, p);

    // x^(n+1.0), p^(n+1.0) -> x^(n+0.5), p^(n+1.0) 
    push_xi(&p, vxi, -.5f * d_dt);
    find_idx_off(p.xi, j, h0, real(0.));
    
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(&p, vxi, d_dt);
    find_idx_off(p.xi, k, h1, real(0.));

    shift0y = j[1] - (ci0[1] + ci1[1]);
    shift0z = j[2] - (ci0[2] + ci1[2]);
    shift1y = k[1] - (ci0[1] + ci1[1]);
    shift1z = k[2] - (ci0[2] + ci1[2]);
  } else {
    *qni_wni = real(0.);
    vxi[0] = real(0.);
    vxi[1] = real(0.);
    vxi[2] = real(0.);
    shift0y = 0;
    shift0z = 0;
    shift1y = 0;
    shift1z = 0;
    h0[1] = real(0.);
    h1[1] = real(0.);
    h0[2] = real(0.);
    h1[2] = real(0.);
  }
  cache_shape_arrays(SHAPE_INFO_PARAMS, h0, h1, shift0y, shift0z,
		     shift1y, shift1z);
}

__global__ static void
push_part_p1_5(int n_particles, particles_cuda_dev_t d_particles,
	       struct shapeinfo *d_shapeinfo, real *d_vxi, real *d_qni,
	       int *d_ci1,
	       int block_stride, int block_start)
{
  int tid = threadIdx.x, bid = blockIdx.x * block_stride + block_start;
  const int cells_per_block = BLOCKSIZE_Y * BLOCKSIZE_Z;

  if (tid == 0) {
    blockIdx_to_cellPos(&d_particles, bid, ci0);
  }
  __syncthreads();

  int cell_begin = d_particles.c_offsets[bid * cells_per_block];
  int cell_end   = d_particles.c_offsets[(bid+1) * cells_per_block];
  
  int nr_loops = (cell_end - cell_begin + THREADS_PER_BLOCK-1) / THREADS_PER_BLOCK;
  int imax = cell_begin + nr_loops * THREADS_PER_BLOCK;
  
  for (int i = cell_begin + tid; i < imax; i += THREADS_PER_BLOCK) {
    DECLARE_SHAPE_INFO;
    int ci1[3];
    real vxi[3], qni_wni;
    calc_shape_info(ci1, i, d_particles, vxi, &qni_wni, SHAPE_INFO_PARAMS, cell_end);
    if (i < cell_end) {
      d_shapeinfo[i] = *si;
      d_vxi[i] = vxi[0];
      d_qni[i] = qni_wni;
      d_ci1[i] = encode_ci1(ci1);
    }
  }
}

// ----------------------------------------------------------------------
// calc_j

__device__ static void
calc_j(int *ci1, int i, particles_cuda_dev_t d_particles,
       real *vxi, real *qni_wni, SHAPE_INFO_ARGS, int cell_end)
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
    find_idx(p.xi, ci1, real(0.));
    ci1[1] -= ci0[1];
    ci1[2] -= ci0[2];

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
    shift0y = j[1] - (ci0[1] + ci1[1]);
    shift0z = j[2] - (ci0[2] + ci1[2]);
    shift1y = k[1] - (ci0[1] + ci1[1]);
    shift1z = k[2] - (ci0[2] + ci1[2]);
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
// current_add

__device__ static void
current_add(int m, int jy, int jz, real val, int ci1[3])
{
  float *addr = &scurr(m, jy + ci1[1], jz + ci1[2]);
#if __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
  atomicAdd(addr, val);
#else
  while ((val = atomicExch(addr, atomicExch(addr, 0.0f)+val))!=0.0f);
#endif
}

// ----------------------------------------------------------------------
// yz_calc_jx

__device__ static void
yz_calc_jx(real vxi, real qni_wni, int ci1[3], SHAPE_INFO_ARGS)
{
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
      
      current_add(0, jy, jz, fnqx * wx, ci1);
    }
  }
}

// ----------------------------------------------------------------------
// yz_calc_jy

__device__ static void
yz_calc_jy(real qni_wni, int ci1[3], SHAPE_INFO_ARGS)
{
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
      current_add(1, jy, jz, last, ci1);
    }
    for (int jy = -1; jy <= 0; jy++) {
      real s0y = pick_shape_coeff(0, y, jy, SI_SHIFT0Y);
      real s1y = pick_shape_coeff(1, y, jy, SI_SHIFT1Y) - s0y;
      real wy = s1y * tmp1;
      last -= fnqy*wy;
      current_add(1, jy, jz, last, ci1);
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
      current_add(1, jy, jz, last, ci1);
    }
  }
}

// ----------------------------------------------------------------------
// yz_calc_jz

__device__ static void
yz_calc_jz(real qni_wni, int ci1[3], SHAPE_INFO_ARGS)
{
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
      current_add(2, jy, jz, last, ci1);
    }
    for (int jz = -1; jz <= 0; jz++) {
      real s0z = pick_shape_coeff(0, z, jz, SI_SHIFT0Z);
      real s1z = pick_shape_coeff(1, z, jz, SI_SHIFT1Z) - s0z;
      real wz = s1z * tmp1;
      last -= fnqz*wz;
      current_add(2, jy, jz, last, ci1);
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
      current_add(2, jy, jz, last, ci1);
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
    real val = scurr(m, jy, jz);
    scratch(0,jy,jz) += val;
    i += THREADS_PER_BLOCK;
  }
}

__device__ static void
add_scurr_to_flds(real *d_flds)
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
    real val = scurr(m, jy, jz);
    F3_DEV(JXI+m, 0,jy+ci0[1],jz+ci0[2]) += val;
    i += THREADS_PER_BLOCK;
  }
}

__global__ static void
push_part_p2(int n_particles, particles_cuda_dev_t d_particles,
	     struct shapeinfo *d_shapeinfo,
	     real *d_vxi, real *d_qni, int *d_ci1,
	     real *d_flds,
#ifdef USE_SCRATCH
	     real *d_scratch,
#endif
	     int block_stride, int block_start)
{
  int tid = threadIdx.x;
  const int cells_per_block = BLOCKSIZE_Y * BLOCKSIZE_Z;

  zero_scurr();

  __shared__ int bid;
  if (tid == 0) {
    bid = blockIdx.x * block_stride + block_start;
    blockIdx_to_cellPos(&d_particles, bid, ci0);
  }
  __syncthreads();

  {
    int cell_begin = d_particles.c_offsets[bid * cells_per_block];
    int cell_end   = d_particles.c_offsets[(bid+1) * cells_per_block];
    
    int nr_loops = (cell_end - cell_begin + THREADS_PER_BLOCK-1) / THREADS_PER_BLOCK;
    int imax = cell_begin + nr_loops * THREADS_PER_BLOCK;
    
    for (int i = cell_begin + tid; i < imax; i += THREADS_PER_BLOCK) {
      DECLARE_SHAPE_INFO;
      int ci1[3];
      real vxi[3];
      real qni_wni;
      calc_j(ci1, i, d_particles, vxi, &qni_wni, SHAPE_INFO_PARAMS, cell_end);
      *si = d_shapeinfo[i];
      vxi[0] = d_vxi[i];
      qni_wni = d_qni[i];
      decode_ci1(d_ci1[i], ci1);
      yz_calc_jx(vxi[0], qni_wni, ci1, SHAPE_INFO_PARAMS);
      yz_calc_jy(qni_wni, ci1, SHAPE_INFO_PARAMS);
      yz_calc_jz(qni_wni, ci1, SHAPE_INFO_PARAMS);
    }
  }

  __syncthreads();
#ifdef USE_SCRATCH
  add_scurr_to_scratch(d_scratch, bid);
#else
  add_scurr_to_flds(d_flds);
#endif
}

#endif

