
// ======================================================================
// field caching

#ifndef EM_CACHE
#define EM_CACHE EM_CACHE_NONE
#endif

#ifndef F3_EM
#define F3_EM _F3
typedef fields_t flds_em_t;
#endif

// FIXME, shared between em and curr cache currently

// OPT, shouldn't we be able to do with less ghosts?
#if DIM == DIM_YZ
#define BLOCKBND_X 0
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#elif DIM == DIM_XYZ
#define BLOCKBND_X 2
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#endif

#define BLOCKGSIZE_X (BLOCKSIZE_X + 2 * BLOCKBND_X)
#define BLOCKGSIZE_Y (BLOCKSIZE_Y + 2 * BLOCKBND_Y)
#define BLOCKGSIZE_Z (BLOCKSIZE_Z + 2 * BLOCKBND_Z)

// ----------------------------------------------------------------------
#if EM_CACHE == EM_CACHE_NONE

typedef flds_em_t em_cache_t;

#ifndef EM_CACHE_DIM
#define EM_CACHE_DIM DIM_XYZ
#endif

#if EM_CACHE_DIM == DIM_XYZ
#define F3_CACHE(em_cache, m, i,j,k)  (F3_EM(em_cache, m, i,j,k))
#elif EM_CACHE_DIM == DIM_XZ
#define F3_CACHE(em_cache, m, i,j,k)  (F3_EM(em_cache, m, i,0,k))
#elif EM_CACHE_DIM == DIM_1
#define F3_CACHE(em_cache, m, i,j,k)  (F3_EM(em_cache, m, 0,0,0))
#else
#error unhandled EM_CACHE_DIM
#endif

CUDA_DEVICE static inline em_cache_t
em_cache_create(flds_em_t flds_em, int ci0[3])
{
  return flds_em;
}

// ----------------------------------------------------------------------
#elif EM_CACHE == EM_CACHE_CUDA

typedef fields_real_t *em_cache_t;

#if DIM == DIM_YZ
#define F3_CACHE(em_cache, m, jx, jy, jz)				\
  ((em_cache)[(((m)							\
	       *BLOCKGSIZE_Z + (jz))					\
	      *BLOCKGSIZE_Y + (jy))])
#elif DIM == DIM_XYZ
#define F3_CACHE(em_cache, m, jx, jy, jz)				\
  ((em_cache)[((((m)							\
		*BLOCKGSIZE_Z + (jz))					\
	       *BLOCKGSIZE_Y + (jy))					\
	      *BLOCKGSIZE_X + (jx))])
#endif

__device__ static em_cache_t
cache_fields(em_cache_t flds_em_block, flds_em_t flds_em, int *ci0)
{
  em_cache_t em_cache = flds_em_block + ((((-EX) * 
				       BLOCKGSIZE_Z + -ci0[2] + BLOCKBND_Z) *
				      BLOCKGSIZE_Y + -ci0[1] + BLOCKBND_Y) *
				     BLOCKGSIZE_X + -ci0[0] + BLOCKBND_X);

  int n = BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z;
  // if we're not actually running on the GPU, we're not multi-threaded, so
  // the caching wouldn't all be initialized first (and then __syncthreads()),
  // so instead we have the first "thread" do all of the caching.
#ifdef __CUDACC__
  for (int ti = threadIdx.x; ti < n; ti += THREADS_PER_BLOCK) {
#else
  if (threadIdx.x == 0) for (int ti = 0; ti < n; ti ++) {
#endif
    int tmp = ti;
#if DIM == DIM_XYZ
    int jx = tmp % BLOCKGSIZE_X - BLOCKBND_X;
    tmp /= BLOCKGSIZE_X;
#endif
    int jy = tmp % BLOCKGSIZE_Y - BLOCKBND_Y;
    tmp /= BLOCKGSIZE_Y;
    int jz = tmp % BLOCKGSIZE_Z - BLOCKBND_Z;
    // OPT? currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3_CACHE(em_cache, m, jx+ci0[0],jy+ci0[1],jz+ci0[2]) = 
	F3_EM(flds_em, m, jx+ci0[0],jy+ci0[1],jz+ci0[2]);
    }
  }
  return em_cache;
}

CUDA_SHARED fields_real_t flds_em_block[6 * BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z];

CUDA_DEVICE static inline em_cache_t
em_cache_create(flds_em_t flds_em, int ci0[3])
{
  return cache_fields(flds_em_block, flds_em, ci0);
}

#endif

