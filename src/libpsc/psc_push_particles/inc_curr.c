
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

#define DECLARE_CURR_CACHE(flds_curr, d_flds, size, ci0)		\
  flds_curr_t flds_curr = {						\
    .arr_shift = d_flds + ((((0)					\
			     * prm.mx[2] - prm.ilg[2])			\
			    * prm.mx[1] - prm.ilg[1])			\
			   * prm.mx[0] - prm.ilg[0]) }

static inline void
curr_cache_add(flds_curr_t flds_curr, fields_real_t *d_flds, int ci0[3])
{
}

// ----------------------------------------------------------------------
#elif CURR_CACHE == CURR_CACHE_CUDA

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

typedef fields_real_t * flds_curr_t;

CUDA_DEVICE static inline void
curr_add(flds_curr_t flds_curr, int m, int jx, int jy, int jz, real val)
{
  real *addr = &F3_DEV_SHIFT(flds_curr, JXI+m, jx,jy,jz);
#ifdef __CUDACC__
  atomicAdd(addr, val);
#else
  *addr += val;
#endif
}

#define DECLARE_CURR_CACHE(flds_curr, d_flds, size, ci0)		\
  unsigned int curr_size = 3 * prm.mx[0] * prm.mx[1] * prm.mx[2];	\
  CUDA_SHARED fields_real_t flds_curr_shared[curr_size];		\
  for (int i = 0; i < curr_size; i++) {					\
    flds_curr_shared[i] = 0.f;						\
  }									\
  flds_curr_t flds_curr = flds_curr_shared + ((((0)			\
						* prm.mx[2] - ci0[2] - prm.ilg[2]) \
					       * prm.mx[1] - ci0[1] - prm.ilg[1]) \
					      * prm.mx[0] - ci0[0] - prm.ilg[0]) 
//flds_curr_shared

static inline void
curr_cache_add(flds_curr_t flds_curr, fields_real_t *d_flds, int ci0[3])
{
  CUDA_SYNCTHREADS();

  for (int m = 0; m < 3; m++) {
    for (int iz = -BLOCKBND_Z; iz < BLOCKGSIZE_Z + BLOCKBND_Z; iz++) {
      for (int iy = -BLOCKBND_Y; iy < BLOCKGSIZE_Y + BLOCKBND_Y; iy++) {
	for (int ix = -BLOCKBND_X; ix < BLOCKGSIZE_X + BLOCKBND_X; ix++) {
	  F3_DEV(d_flds, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]) +=
	    F3_DEV_SHIFT(flds_curr, JXI + m, ix+ci0[0],iy+ci0[1],iz+ci0[2]);
	}
      }
    }
  }
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

