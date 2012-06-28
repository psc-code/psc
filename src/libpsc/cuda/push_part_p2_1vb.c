
// FIXME -> common.c

__device__ static void
find_idx_off_pos_1st(const real xi[3], int j[3], real h[3], real pos[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    pos[d] = xi[d] * d_consts.dxi[d] + shift;
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

#define xBLOCKSTRIDE ((((BLOCKSIZE_Y + 2*SW) * (BLOCKSIZE_Z + 2*SW) + 31) / 32) * 32)
__shared__ real scurr1[WARPS_PER_BLOCK * xBLOCKSTRIDE];
__shared__ real scurr2[WARPS_PER_BLOCK * xBLOCKSTRIDE];
__shared__ real scurr3[WARPS_PER_BLOCK * xBLOCKSTRIDE];

#define w_scurr_(scurr, wid, jy, jz) (scurr[((jz)+SW) * (BLOCKSIZE_Y + 2*SW) + (jy)+SW \
					     + (wid) * xBLOCKSTRIDE])

#define w_scurr(scurr, jy, jz) w_scurr_(scurr, threadIdx.x >> 5, jy, jz)

__shared__ int ci0[3]; // cell index of lower-left cell in block

// ======================================================================

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
    real fnqx = vxi[0] * p.qni_wni * d_consts.fnqs;
    
    int lf[3];
    real of[3];
    find_idx_off_1st(p.xi, lf, of, real(0.), d_consts.dxi);
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
    dx1[0] = .5f * off[0] - x[0];
    if (dx[0] != 0.f)
      dx1[1] = dx[1] / dx[0] * dx1[0];
    else
      dx1[1] = 0.f;
  } else {
    dx1[1] = .5f * off[1] - x[1];
    if (dx[1] != 0.f)
      dx1[0] = dx[0] / dx[1] * dx1[1];
    else
      dx1[0] = 0.f;
  }
}

__device__ static void
curr_2d_vb_cell(int i[2], real x[2], real dx[2], real qni_wni,
		real dxt[2], int off[2])
{
  real fnqy = qni_wni * d_consts.fnqys;
  real fnqz = qni_wni * d_consts.fnqzs;
  current_add(scurr2, i[0],i[1]  , fnqy * dx[0] * (.5f - x[1] - .5f * dx[1]));
  current_add(scurr2, i[0],i[1]+1, fnqy * dx[0] * (.5f + x[1] + .5f * dx[1]));
  current_add(scurr3, i[0],i[1]  , fnqz * dx[1] * (.5f - x[0] - .5f * dx[0]));
  current_add(scurr3, i[0]+1,i[1], fnqz * dx[1] * (.5f + x[0] + .5f * dx[0]));
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
yz_calc_jyjz(int i, particles_cuda_dev_t d_particles)
{
  struct d_particle p;

  if (do_read) {
    LOAD_PARTICLE(p, d_particles, i);
  }

  if (do_calc_jy) {
    real vxi[3];
    real h0[3], h1[3];
    real xm[3], xp[3];
    
    int j[3], k[3];
    calc_vxi(vxi, p);
    
    // x^(n+1.0), p^(n+1.0) -> x^(n+0.5), p^(n+1.0) 
    push_xi(&p, vxi, -.5f * d_consts.dt);
    find_idx_off_pos_1st(p.xi, j, h0, xm, real(0.));
    
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    push_xi(&p, vxi, d_consts.dt);
    find_idx_off_pos_1st(p.xi, k, h1, xp, real(0.));
    
    int idiff[2] = { k[1] - j[1], k[2] - j[2] };
    real dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };
    real x[2] = { xm[1] - (j[1] + real(.5)), xm[2] - (j[2] + real(.5)) };
    int i[2] = { j[1] - ci0[1], j[2] - ci0[2] };
  
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
      real dx1[2];
      dx1[0] = .5f * idiff[0] - x[0];
      dx1[1] = dx[1] / dx[0] * dx1[0];
      if (fabsf(x[1] + dx1[1]) > .5f) {
	first_dir = 1;
      } else {
	first_dir = 0;
      }
      second_dir = 1 - first_dir;
    }
    
    if (first_dir >= 0) {
      real dx1[2];
      off[1-first_dir] = 0;
      off[first_dir] = idiff[first_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(i, x, dx1, p.qni_wni, dx, off);
    }
    
    if (second_dir >= 0) {
      real dx1[2];
      off[first_dir] = 0;
      off[second_dir] = idiff[second_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(i, x, dx1, p.qni_wni, dx, off);
    }
    
    curr_2d_vb_cell(i, x, dx, p.qni_wni, NULL, NULL);
  }
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

#define push_part_p2_2							\

// ======================================================================

__global__ static void
push_part_p2x(int n_particles, particles_cuda_dev_t d_particles, real *d_flds,
	      int block_stride, int block_start)
{
  do_read = true;
  do_reduce = true;
  do_write = true;
  do_calc_jx = true;
  do_calc_jy = true;
  do_calc_jz = true;

  if (do_write) {
    zero_scurr(scurr1);
    zero_scurr(scurr2);
    zero_scurr(scurr3);
  }

  int tid = threadIdx.x;
  int block_pos[3];
  block_pos[1] = blockIdx.x * 2;
  block_pos[2] = blockIdx.y * 2;
  block_pos[1] += block_start & 1;
  block_pos[2] += block_start >> 1;
  if (block_pos[1] >= d_consts.b_mx[1] ||
      block_pos[2] >= d_consts.b_mx[2])
    return;

  int bid = block_pos_to_block_idx(block_pos, d_consts.b_mx);
  __shared__ int s_block_end;
  if (tid == 0) {
    ci0[0] = 0;
    ci0[1] = block_pos[1] * BLOCKSIZE_Y;
    ci0[2] = block_pos[2] * BLOCKSIZE_Z;
    s_block_end = d_particles.offsets[bid + 1];
  }
  __syncthreads();

  int block_begin = d_particles.offsets[bid];

  for (int i = block_begin + tid; i < s_block_end; i += THREADS_PER_BLOCK) {
    yz_calc_jx(i, d_particles);
    yz_calc_jyjz(i, d_particles);
  }
  
  if (do_write) {
    __syncthreads();
    add_scurr_to_flds1(scurr1, d_flds, 0);
    add_scurr_to_flds1(scurr2, d_flds, 1);
    add_scurr_to_flds1(scurr3, d_flds, 2);
  }
}
