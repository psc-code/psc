
#define F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz, wid)				\
  (((((fldnr)								\
      *BLOCKGSIZE_Z + (jz))						\
     *BLOCKGSIZE_Y + (jy))						\
    *BLOCKGSIZE_X + (jx))						\
   *CURR_CACHE_N_REDUNDANT + (wid))

#define F3_DEV_SHIFT(d_flds, fldnr, jx,jy,jz, wid)	\
  ((d_flds)[F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz, wid)])

#define CURR_CACHE_SIZE (3 * BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z * CURR_CACHE_N_REDUNDANT)

typedef fields_real_t * flds_curr_t;

CUDA_DEVICE static inline flds_curr_t
flds_curr_shift(flds_curr_t flds_curr, int m, int dx, int dy, int dz)
{
  return flds_curr + F3_DEV_SHIFT_OFF(m, dx,dy,dz, 0);
}

CUDA_DEVICE static fields_real_t *
init_curr_cache(fields_real_t *flds_curr_block, int ci0[3])
{
#ifdef __CUDACC__
  for (int i = threadIdx.x; i < CURR_CACHE_SIZE; i += THREADS_PER_BLOCK) {
    flds_curr_block[i] = 0.f;
  }
#else
  if (threadIdx.x == 0) {
    for (int i = 0; i < CURR_CACHE_SIZE; i++) {
      flds_curr_block[i] = 0.f;
    }
  }
#endif
			 
  return flds_curr_shift(flds_curr_block, -JXI,
			 -ci0[0] + BLOCKBND_X,
			 -ci0[1] + BLOCKBND_Y,
			 -ci0[2] + BLOCKBND_Z);
}

#if CURR_CACHE_GMEM
#define NR_BLOCKS ((512/4) * (512/4))

__device__ static fields_real_t flds_curr_blocks[CURR_CACHE_SIZE * NR_BLOCKS];

#define DECLARE_CURR_CACHE(d_flds, ci0)					\
  ({									\
    assert(find_bid() < NR_BLOCKS);					\
    init_curr_cache(flds_curr_blocks + find_bid() * CURR_CACHE_SIZE, ci0); \
  })

#else

CUDA_SHARED fields_real_t flds_curr_block[CURR_CACHE_SIZE];

#define DECLARE_CURR_CACHE(d_flds, ci0)					\
  ({									\
    init_curr_cache(flds_curr_block, ci0);				\
  })

#endif

CUDA_DEVICE static inline void
curr_add(flds_curr_t flds_curr, int m, int jx, int jy, int jz, real val)
{
  real *addr = &F3_DEV_SHIFT(flds_curr, m, jx,jy,jz, threadIdx.x % CURR_CACHE_N_REDUNDANT);
#ifdef __CUDACC__
  atomicAdd(addr, val);
#else
  *addr += val;
#endif
}

CUDA_DEVICE static void
curr_cache_add(flds_curr_t flds_curr, fields_real_t *d_flds, int ci0[3])
{
  CUDA_SYNCTHREADS();

#ifdef __CUDACC__
  for (int i = threadIdx.x; i < BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z; i += THREADS_PER_BLOCK) {
    int rem = i;
    int ix = rem % BLOCKGSIZE_X;
    rem /= BLOCKGSIZE_X;
    int iy = rem % BLOCKGSIZE_Y;
    rem /= BLOCKGSIZE_Y;
    int iz = rem;
    ix -= BLOCKBND_X;
    iy -= BLOCKBND_Y;
    iz -= BLOCKBND_Z;
    for (int m = 0; m < 3; m++) {
      fields_real_t val = 0.f;
      for (int wid = 0; wid < CURR_CACHE_N_REDUNDANT; wid++) {
	val += F3_DEV_SHIFT(flds_curr, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2], wid);
      }
      fields_real_t *addr = &F3_DEV(d_flds, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]);
      atomicAdd(addr, val);
    }
  }
#else
  if (threadIdx.x != THREADS_PER_BLOCK - 1) {
    return;
  }
  for (int m = 0; m < 3; m++) {
    for (int iz = -BLOCKBND_Z; iz < BLOCKSIZE_Z + BLOCKBND_Z; iz++) {
      for (int iy = -BLOCKBND_Y; iy < BLOCKSIZE_Y + BLOCKBND_Y; iy++) {
	for (int ix = -BLOCKBND_X; ix < BLOCKSIZE_X + BLOCKBND_X; ix++) {
	  fields_real_t val = 0.f;
	  for (int wid = 0; wid < CURR_CACHE_N_REDUNDANT; wid++) {
	    val += F3_DEV_SHIFT(flds_curr, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2], wid);
	  }
	  F3_DEV(d_flds, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]) += val;
	}
      }
    }
  }
#endif
}

