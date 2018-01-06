
#define CURR_CACHE_HAVE_SHIFT

#define F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz, wid)				\
  (((((fldnr)								\
      *BLOCKGSIZE_Z + (jz))						\
     *BLOCKGSIZE_Y + (jy))						\
    *BLOCKGSIZE_X + (jx))						\
   *CURR_CACHE_N_REDUNDANT + (wid))

#define F3_DEV_SHIFT(d_flds, fldnr, jx,jy,jz, wid)	\
  ((d_flds)[F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz, wid)])

#define CURR_CACHE_SIZE (3 * BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z * CURR_CACHE_N_REDUNDANT)

typedef fields_t::real_t * curr_cache_t;

CUDA_DEVICE static inline curr_cache_t
curr_cache_shift(curr_cache_t curr_cache, int m, int dx, int dy, int dz)
{
  return curr_cache + F3_DEV_SHIFT_OFF(m, dx,dy,dz, 0);
}

CUDA_DEVICE static fields_t::real_t *
init_curr_cache(fields_t::real_t *curr_cache, int ci0[3])
{
#ifdef __CUDACC__
  for (int i = threadIdx.x; i < CURR_CACHE_SIZE; i += THREADS_PER_BLOCK) {
    curr_cache[i] = 0.f;
  }
#else
  if (threadIdx.x == 0) {
    for (int i = 0; i < CURR_CACHE_SIZE; i++) {
      curr_cache[i] = 0.f;
    }
  }
#endif
			 
  return curr_cache_shift(curr_cache, -JXI,
			  -ci0[0] + BLOCKBND_X,
			  -ci0[1] + BLOCKBND_Y,
			  -ci0[2] + BLOCKBND_Z);
}

#if CURR_CACHE_GMEM
#define NR_BLOCKS ((512/4) * (512/4))

__device__ static fields_t::real_t curr_cache_blocks[CURR_CACHE_SIZE * NR_BLOCKS];

CUDA_DEVICE static int find_bid();

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(flds_curr_t flds_curr, int ci0[3])
{
  assert(find_bid() < NR_BLOCKS);
  return init_curr_cache(curr_cache_blocks + find_bid() * CURR_CACHE_SIZE, ci0);
}

#else

CUDA_SHARED fields_t::real_t curr_cache_block[CURR_CACHE_SIZE];

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(flds_curr_t flds_curr, int ci0[3])
{
  return init_curr_cache(curr_cache_block, ci0);
}

#endif

CUDA_DEVICE static inline void
curr_cache_add(curr_cache_t curr_cache, int m, int jx, int jy, int jz, real val)
{
  int wid = threadIdx.x % CURR_CACHE_N_REDUNDANT;
  real *addr = &F3_DEV_SHIFT(curr_cache, m, jx,jy,jz, wid);
  atomicAdd(addr, val);
}

CUDA_DEVICE static void
curr_cache_destroy(curr_cache_t curr_cache, flds_curr_t flds_curr, int ci0[3])
{
  Fields3d<flds_curr_t> F(flds_curr);

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
      fields_t::real_t val = 0.f;
      for (int wid = 0; wid < CURR_CACHE_N_REDUNDANT; wid++) {
	val += F3_DEV_SHIFT(curr_cache, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2], wid);
      }
      fields_t::real_t *addr = &F(JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]);
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
	  fields_t::real_t val = 0.f;
	  for (int wid = 0; wid < CURR_CACHE_N_REDUNDANT; wid++) {
	    val += F3_DEV_SHIFT(curr_cache, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2], wid);
	  }
	  FJXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]) += val;
	}
      }
    }
  }
#endif
}

