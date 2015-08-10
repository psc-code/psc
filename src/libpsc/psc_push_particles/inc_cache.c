
// ======================================================================
// field caching

// OPT: precalc offsets into fld_cache (including ci[])
// OPT: use more shmem?

#if EM_CACHE == EM_CACHE_NONE

#define F3_CACHE(flds_em, m, jx, jy, jz)	\
  (F3_DEV(flds_em, m, jx,jy,jz))

#define DECLARE_EM_CACHE(flds_em, d_flds, size, ci0)	\
  real *flds_em = d_flds

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
	       *BLOCKGSIZE_Z + ((jz)-ci0[2]))				\
	      *BLOCKGSIZE_Y + ((jy)-ci0[1]))])
#elif DIM == DIM_XYZ
#define F3_CACHE(flds_em, m, jx, jy, jz)				\
  ((flds_em)[((((m)							\
		*BLOCKGSIZE_Z + ((jz-ci0[2])))				\
	       *BLOCKGSIZE_Y + ((jy-ci0[1])))				\
	      *BLOCKGSIZE_X + ((jx-ci0[0])))])
#endif

__device__ static float *
cache_fields(float *flds_em_shared, float *d_flds, int size, int *ci0)
{
  float *flds_em = flds_em_shared + ((((-EX) * 
				       BLOCKGSIZE_Z + BLOCKBND_Z) *
				      BLOCKGSIZE_Y + BLOCKBND_Y) *
				     BLOCKGSIZE_X + BLOCKBND_X);
  int ti = threadIdx.x;
  int n = BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z;
  while (ti < n) {
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
    ti += THREADS_PER_BLOCK;
  }
  return flds_em;
}

#define DECLARE_EM_CACHE(flds_em, d_flds, size, ci0)	\
  __shared__ real flds_em_shared[6 * BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z]; \
  float *flds_em = cache_fields(flds_em_shared, d_flds, size, ci0)

#endif

