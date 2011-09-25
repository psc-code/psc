
// FIXME -> common.c

__device__ static void
find_idx_off_pos_1st(const real xi[3], int j[3], real h[3], real pos[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    pos[d] = xi[d] * d_dxi[d] + shift;
    j[d] = cuda_fint(pos[d]);
    h[d] = pos[d] - j[d];
  }
}


__shared__ volatile bool do_read;
__shared__ volatile bool do_write;
__shared__ volatile bool do_reduce;
__shared__ volatile bool do_calc_jx;
__shared__ volatile bool do_calc_jy;
__shared__ volatile bool do_calc_jz;

// OPT: take i < cell_end condition out of load
// OPT: reduce two at a time
// OPT: try splitting current calc / measuring by itself
// OPT: get rid of block_stride

#define WARPS_PER_BLOCK (THREADS_PER_BLOCK / 32)

__shared__ int _cell_end[WARPS_PER_BLOCK]; // last particle in current cell valid in p2x

#define xBLOCKSTRIDE ((((BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW) + 31) / 32) * 32)
__shared__ real scurr1[WARPS_PER_BLOCK * xBLOCKSTRIDE];
__shared__ real scurr2[WARPS_PER_BLOCK * xBLOCKSTRIDE];

#define w_scurr_(scurr, wid, jy, jz) (scurr[((jz)+SW) * (BLOCKSIZE_Y + 2*SW) + (jy)+SW \
					     + (wid) * xBLOCKSTRIDE])

#define w_scurr(scurr, jy, jz) w_scurr_(scurr, threadIdx.x >> 5, jy, jz)

#define w_cell_end (_cell_end[threadIdx.x >> 5])

#if CACHE_SHAPE_ARRAYS == 7

struct shapeinfo_h {
  real xm[2], xp[2];
};

struct shapeinfo_i {
  char4 shiftyz; // x: y[0] y: y[1] z: z[0] w: z[1]
};

#define DECLARE_SHAPE_INFO			\
  struct shapeinfo_h _si_h, *si_h = &_si_h;	\
  struct shapeinfo_i _si_i, *si_i = &_si_i

#define SHAPE_INFO_PARAMS si_h, si_i
#define SHAPE_INFO_ARGS struct shapeinfo_h *si_h, struct shapeinfo_i *si_i

#define D_SHAPEINFO_PARAMS d_si_h, d_si_i
#define D_SHAPEINFO_ARGS struct shapeinfo_h *d_si_h, struct shapeinfo_i *d_si_i

#define SI_SHIFT0Y si_i->shiftyz.x
#define SI_SHIFT1Y si_i->shiftyz.y
#define SI_SHIFT10Y (si_i->shiftyz.y - si_i->shiftyz.x)
#define SI_SHIFT0Z si_i->shiftyz.z
#define SI_SHIFT1Z si_i->shiftyz.w
#define SI_SHIFT10Z (si_i->shiftyz.w - si_i->shiftyz.z)

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *xm, real *xp,
		   short int shift0y, short int shift0z,
		   short int shift1y, short int shift1z)
{
  si_h->xm[0] = xm[1];
  si_h->xm[1] = xm[2];
  si_h->xp[0] = xp[1];
  si_h->xp[1] = xp[2];
  SI_SHIFT0Y = shift0y;
  SI_SHIFT1Y = shift1y;
  SI_SHIFT0Z = shift0z;
  SI_SHIFT1Z = shift1z;
}

// ======================================================================
#else
#error
#endif

__shared__ int ci0[3]; // cell index of lower-left cell in block

// ======================================================================

// ----------------------------------------------------------------------
// calc_shape_info

__device__ static void
calc_shape_info(int i, particles_cuda_dev_t d_particles,
		real *vxi, real *qni_wni, SHAPE_INFO_ARGS)
{
  short int shift0y;
  short int shift0z;
  short int shift1y;
  short int shift1z;
  real h0[3], h1[3];
  real xm[3], xp[3];

  struct d_particle p;
  LOAD_PARTICLE(p, d_particles, i);
  *qni_wni = p.qni_wni;
  
  int j[3], k[3];
  calc_vxi(vxi, p);
  
  // x^(n+1.0), p^(n+1.0) -> x^(n+0.5), p^(n+1.0) 
  push_xi(&p, vxi, -.5f * d_dt);
  find_idx_off_pos_1st(p.xi, j, h0, xm, real(0.));
  
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(&p, vxi, d_dt);
  find_idx_off_pos_1st(p.xi, k, h1, xp, real(0.));
  
  shift0y = j[1] - ci0[1];
  shift0z = j[2] - ci0[2];
  shift1y = k[1] - ci0[1];
  shift1z = k[2] - ci0[2];

  cache_shape_arrays(SHAPE_INFO_PARAMS, xm, xp, shift0y, shift0z,
		     shift1y, shift1z);
}

// ======================================================================

// ----------------------------------------------------------------------
// current_add

__device__ static void
current_add(real *scurr, int jy, int jz, real val)
{
  int lid = threadIdx.x & 31;
  float *addr = &w_scurr(scurr, jy, jz);
  if (!do_write) {
    if (do_write)
      *addr += val;
  } else if (do_reduce) {
#if __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
    atomicAdd(addr, val);
#else
#if 0
    while ((val = atomicExch(addr, atomicExch(addr, 0.0f)+val))!=0.0f);
#else
    for (int i = 0; i < 32; i++) {
      if (lid == i) {
	*addr += val;
      }
    }
#endif
#endif
  } else {
    *addr += val;
  }
}

// ----------------------------------------------------------------------
// yz_calc_jx

__device__ static void
yz_calc_jx(int i, particles_cuda_dev_t d_particles)
{
  struct d_particle p;
  if (do_read) {
    LOAD_PARTICLE(p, d_particles, i);
  }

  if (do_calc_jx) {
    real vxi[3];
    calc_vxi(vxi, p);
    real fnqx = vxi[0] * p.qni_wni * d_fnqs;
    
    int lf[3];
    real of[3];
    find_idx_off_1st(p.xi, lf, of, real(0.));
    lf[1] -= ci0[1];
    lf[2] -= ci0[2];
    current_add(scurr1, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnqx);
    current_add(scurr1, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnqx);
    current_add(scurr1, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnqx);
    current_add(scurr1, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnqx);
  }
}

// ----------------------------------------------------------------------
// yz_calc_jy

__device__ static void
calc_dx1(real dx1[2], real x[2], real dx[2], int off[2])
{
  if (off[1] == 0) {
    dx1[0] = .5 * off[0] - x[0];
    dx1[1] = dx[1] / dx[0] * dx1[0];
  } else {
    dx1[1] = .5 * off[1] - x[1];
    dx1[0] = dx[0] / dx[1] * dx1[1];
  }
}

__device__ static void
curr_2d_vb_cell(int i[2], real x[2], real dx[2], real qni_wni,
		real dxt[2], int off[2])
{
  real fnqy = qni_wni * d_fnqys;
  real fnqz = qni_wni * d_fnqzs;
  current_add(scurr1, i[0],i[1]  , fnqy * dx[0] * (.5f - x[1] - .5f * dx[1]));
  current_add(scurr1, i[0],i[1]+1, fnqy * dx[0] * (.5f + x[1] + .5f * dx[1]));
  current_add(scurr2, i[0],i[1]  , fnqz * dx[1] * (.5f - x[0] - .5f * dx[0]));
  current_add(scurr2, i[0]+1,i[1], fnqz * dx[1] * (.5f + x[0] + .5f * dx[0]));
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    x[0] += dx[0] - off[0];
    x[1] += dx[1] - off[1];
    i[0] += off[0];
    i[1] += off[1];
  }
}

__device__ static void
yz_calc_jyjz(real qni_wni, SHAPE_INFO_ARGS)
{
  int i[2] = { SI_SHIFT0Y + ci0[1], SI_SHIFT0Z + ci0[2] };
  int idiff[2] = { SI_SHIFT1Y - SI_SHIFT0Y, SI_SHIFT1Z - SI_SHIFT0Z };
  real *xp = si_h->xp, *xm = si_h->xm;
  real dx[2] = { xp[0] - xm[0], xp[1] - xm[1] };
  real x[2] = { xm[0] - (i[0] + real(.5)), xm[1] - (i[1] + real(.5)) };

  i[0] -= ci0[1]; i[1] -= ci0[2];

  real dx1[2];
  int off[2];
  int first_dir, second_dir = -1;
  // FIXME, make sure we never div-by-zero?
  if (idiff[0] == 0 && idiff[1] == 0) {
    first_dir = -1;
  } else if (idiff[0] == 0) {
    first_dir = 1;
  } else if (idiff[1] == 0) {
    first_dir = 0;
  } else {
    dx1[0] = .5 * idiff[0] - x[0];
    dx1[1] = dx[1] / dx[0] * dx1[0];
    if (fabsf(x[1] + dx1[1]) > .5f) {
      first_dir = 1;
    } else {
      first_dir = 0;
    }
    second_dir = 1 - first_dir;
  }

  if (first_dir >= 0) {
    off[1-first_dir] = 0;
    off[first_dir] = idiff[first_dir];
    calc_dx1(dx1, x, dx, off);
    curr_2d_vb_cell(i, x, dx1, qni_wni, dx, off);
  }

  if (second_dir >= 0) {
    off[first_dir] = 0;
    off[second_dir] = idiff[second_dir];
    calc_dx1(dx1, x, dx, off);
    curr_2d_vb_cell(i, x, dx1, qni_wni, dx, off);
  }
    
  curr_2d_vb_cell(i, x, dx, qni_wni, NULL, NULL);
}

__device__ static void
zero_scurr(real *scurr)
{
  int i = threadIdx.x;
  int N = xBLOCKSTRIDE * WARPS_PER_BLOCK;
  while (i < N) {
    scurr[i] = real(0.);
    i += THREADS_PER_BLOCK;
  }
}

__device__ static void
add_scurr_to_flds1(real *scurr, real *d_flds, int m)
{
  int i = threadIdx.x;
  int stride = (BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW);
  while (i < stride) {
    int rem = i;
    int jz = rem / (BLOCKSIZE_Y + 2*SW);
    rem -= jz * (BLOCKSIZE_Y + 2*SW);
    int jy = rem;
    jz -= SW;
    jy -= SW;
    real val = real(0.);
    for (int wid = 0; wid < WARPS_PER_BLOCK; wid++) {
      val += w_scurr_(scurr, wid, jy, jz);
    }
    F3_DEV(JXI+m, 0,jy+ci0[1],jz+ci0[2]) += val;
    i += THREADS_PER_BLOCK;
  }
}

#define push_part_p2_0							\
  do_read = true;							\
  do_reduce = true;							\
  do_write = true;							\
  do_calc_jx = true;							\
  do_calc_jy = true;							\
  do_calc_jz = true

#define push_part_p2_1							\
  int tid = threadIdx.x;						\
  const int cells_per_block = BLOCKSIZE_Y * BLOCKSIZE_Z;		\
									\
  __shared__ int bid;							\
  if (tid == 0) {							\
    bid = blockIdx.x * block_stride + block_start;			\
    blockIdx_to_cellPos(&d_particles, bid, ci0);			\
  }									\
  __syncthreads();							\
									\
  for (int cid = bid * cells_per_block + (tid >> 5);			\
       cid < (bid + 1) * cells_per_block; cid += WARPS_PER_BLOCK)

#define push_part_p2_2							\
  int cell_begin = d_particles.c_offsets[cid];				\
  w_cell_end = d_particles.c_offsets[cid+1];				\
  									\
  for (int i = cell_begin + (tid & 31); i < w_cell_end; i += 32)

// ======================================================================

__global__ static void
push_part_p2x(int n_particles, particles_cuda_dev_t d_particles, real *d_flds,
	      int block_stride, int block_start)
{
  push_part_p2_0;

  if (do_write) {
    zero_scurr(scurr1);
  }

  push_part_p2_1 {
    push_part_p2_2 {
      yz_calc_jx(i, d_particles);
    }
  }

  if (do_write) {
    __syncthreads();
    add_scurr_to_flds1(scurr1, d_flds, 0);
  }
}

// ======================================================================

__global__ static void
push_part_p2y(int n_particles, particles_cuda_dev_t d_particles, real *d_flds,
	      int block_stride, int block_start)
{
  push_part_p2_0;

  if (do_write) {
    zero_scurr(scurr1);
    zero_scurr(scurr2);
  }

  push_part_p2_1 {
    push_part_p2_2 {
      DECLARE_SHAPE_INFO;
      real vxi[3];
      real qni_wni;
      if (do_read) {
	calc_shape_info(i, d_particles, vxi, &qni_wni, SHAPE_INFO_PARAMS);
      }
      if (do_calc_jy) {
	yz_calc_jyjz(qni_wni, SHAPE_INFO_PARAMS);
      }
    }
  }
  if (do_write) {
    __syncthreads();
    add_scurr_to_flds1(scurr1, d_flds, 1);
    add_scurr_to_flds1(scurr2, d_flds, 2);
  }
}

// ======================================================================

__global__ static void
push_part_p2z(int n_particles, particles_cuda_dev_t d_particles, real *d_flds,
	      int block_stride, int block_start)
{
}



