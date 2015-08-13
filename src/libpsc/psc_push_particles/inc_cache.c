
// ======================================================================
// field caching

// OPT: use more shmem?

#if EM_CACHE == EM_CACHE_NONE

#define F3_CACHE(flds_em, m, jx, jy, jz)	\
  (F3_DEV(flds_em, m, jx,jy,jz))

#define DECLARE_EM_CACHE(flds_em, d_flds, size, ci0)	\
  fields_real_t *flds_em = d_flds

#elif EM_CACHE == EM_CACHE_CUDA

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

#if DIM == DIM_YZ
#define F3_CACHE(flds_em, m, jx, jy, jz)				\
  ((flds_em)[(((m)							\
	       *BLOCKGSIZE_Z + (jz))					\
	      *BLOCKGSIZE_Y + (jy))])
#elif DIM == DIM_XYZ
#define F3_CACHE(flds_em, m, jx, jy, jz)				\
  ((flds_em)[((((m)							\
		*BLOCKGSIZE_Z + (jz))					\
	       *BLOCKGSIZE_Y + (jy))					\
	      *BLOCKGSIZE_X + (jx))])
#endif

__device__ static float *
cache_fields(float *flds_em_shared, float *d_flds, int *ci0)
{
  float *flds_em = flds_em_shared + ((((-EX) * 
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
      F3_CACHE(flds_em, m, jx+ci0[0],jy+ci0[1],jz+ci0[2]) = 
	F3_DEV(d_flds, m, jx+ci0[0],jy+ci0[1],jz+ci0[2]);
    }
  }
  return flds_em;
}

#define DECLARE_EM_CACHE(flds_em, d_flds, ci0)	\
  __shared__ real flds_em_shared[6 * BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z]; \
  float *flds_em = cache_fields(flds_em_shared, d_flds, ci0)

#endif

