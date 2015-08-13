
// ----------------------------------------------------------------------
// curr_add

#if PSC_FIELDS_AS_CUDA2

// ----------------------------------------------------------------------
#if CURR_CACHE == CURR_CACHE_NONE

#if DIM == DIM_YZ

#define F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz)				\
  ((((fldnr)								\
     *prm.mx[2] + (jz))							\
    *prm.mx[1] + (jy)))

#else

#define F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz)				\
  ((((fldnr)								\
     *prm.mx[2] + (jz))							\
    *prm.mx[1] + (jy))							\
   *prm.mx[0] + (jx))

#endif

#define F3_DEV_SHIFT(d_flds, fldnr, jx,jy,jz)	\
  ((d_flds)[F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz)])

typedef struct { fields_real_t * arr_shift; } flds_curr_t;

CUDA_DEVICE static inline void
curr_add(flds_curr_t flds_curr, int m, int jx, int jy, int jz, real val)
{
  real *addr = &F3_DEV_SHIFT(flds_curr.arr_shift, JXI+m, jx,jy,jz);
#ifdef __CUDACC__
  atomicAdd(addr, val);
#else
  *addr += val;
#endif
}

#define DECLARE_CURR_CACHE(d_flds, ci0)					\
  ({ (flds_curr_t) {							\
    .arr_shift = d_flds + ((((0)					\
			     * prm.mx[2] - prm.ilg[2])			\
			    * prm.mx[1] - prm.ilg[1])			\
			   * prm.mx[0] - prm.ilg[0]) };			\
  })

static inline void
curr_cache_add(flds_curr_t flds_curr, fields_real_t *d_flds, int ci0[3])
{
}

// ----------------------------------------------------------------------
#elif CURR_CACHE == CURR_CACHE_CUDA

#if DIM == DIM_YZ

#define F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz)				\
  ((((fldnr)								\
     *BLOCKGSIZE_Z + (jz))						\
    *BLOCKGSIZE_Y + (jy)))

#else

#define F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz)				\
  ((((fldnr)								\
     *BLOCKGSIZE_Z + (jz))						\
    *BLOCKGSIZE_Y + (jy))						\
   *BLOCKGSIZE_X + (jx))

#endif

#define F3_DEV_SHIFT(d_flds, fldnr, jx,jy,jz)	\
  ((d_flds)[F3_DEV_SHIFT_OFF(fldnr, jx,jy,jz)])

typedef fields_real_t * flds_curr_t;

CUDA_DEVICE static inline void
curr_add(flds_curr_t flds_curr, int m, int jx, int jy, int jz, real val)
{
  real *addr = &F3_DEV_SHIFT(flds_curr, m, jx,jy,jz);
#ifdef __CUDACC__
  atomicAdd(addr, val);
#else
  *addr += val;
#endif
}

#define CURR_CACHE_SIZE (3 * BLOCKGSIZE_X * BLOCKGSIZE_Y * BLOCKGSIZE_Z)

CUDA_DEVICE static fields_real_t *
init_curr_cache(fields_real_t *flds_curr_shared, int ci0[3])
{
#ifdef __CUDACC__
#warning TBD
#else
  if (threadIdx.x == 0) {
    for (int i = 0; i < CURR_CACHE_SIZE; i++) {
      flds_curr_shared[i] = 0.f;
    }
  }
#endif
  return flds_curr_shared + ((((-JXI)
			       * BLOCKGSIZE_Z - ci0[2] + BLOCKBND_Z)
			      * BLOCKGSIZE_Y - ci0[1] + BLOCKBND_Y)
			     * BLOCKGSIZE_X - ci0[0] + BLOCKBND_X);
}

#define DECLARE_CURR_CACHE(d_flds, ci0)					\
  ({									\
    CUDA_SHARED fields_real_t flds_curr_shared[CURR_CACHE_SIZE];	\
    init_curr_cache(flds_curr_shared, ci0);				\
  })

static inline void
curr_cache_add(flds_curr_t flds_curr, fields_real_t *d_flds, int ci0[3])
{
  CUDA_SYNCTHREADS();

#ifdef __CUDACC__
#warning TBD
#else
  if (threadIdx.x != THREADS_PER_BLOCK - 1) {
    return;
  }
  for (int m = 0; m < 3; m++) {
    for (int iz = -BLOCKBND_Z; iz < BLOCKSIZE_Z + BLOCKBND_Z; iz++) {
      for (int iy = -BLOCKBND_Y; iy < BLOCKSIZE_Y + BLOCKBND_Y; iy++) {
	for (int ix = -BLOCKBND_X; ix < BLOCKSIZE_X + BLOCKBND_X; ix++) {
	  F3_DEV(d_flds, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]) +=
	    F3_DEV_SHIFT(flds_curr, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]);
	}
      }
    }
  }
#endif
}

#endif

// ----------------------------------------------------------------------
#else

typedef struct psc_fields * flds_curr_t;

CUDA_DEVICE static inline void
curr_add(flds_curr_t flds_curr, int m, int jx, int jy, int jz, real val)
{
  F3_CURR(flds_curr, JXI+m, jx,jy,jz) += val;
}

#endif

// ----------------------------------------------------------------------

#if CALC_J == CALC_J_1VB_SPLIT
#include "inc_curr_1vb_split.c"
#elif CALC_J == CALC_J_1VB_VAR1
#include "inc_curr_1vb_var1.c"
#elif CALC_J == CALC_J_1VB_2D
#include "inc_curr_1vb_2d.c"
#endif

