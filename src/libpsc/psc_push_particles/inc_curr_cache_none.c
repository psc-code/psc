

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

typedef struct { fields_real_t * arr_shift; } curr_cache_t;

CUDA_DEVICE static inline void
curr_cache_add(curr_cache_t curr_cache, int m, int jx, int jy, int jz, real val)
{
  real *addr = &F3_DEV_SHIFT(curr_cache.arr_shift, JXI+m, jx,jy,jz);
  atomicAdd(addr, val);
}

CUDA_DEVICE static inline curr_cache_t
curr_cache_shift(curr_cache_t curr_cache, int m, int dx, int dy, int dz)
{
  return (curr_cache_t) { .arr_shift = curr_cache.arr_shift
      + F3_DEV_SHIFT_OFF(m, dx, dy, dz) };
}

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(float *curr_cache, int ci0[3])
{
  return curr_cache_shift((curr_cache_t) { .arr_shift = curr_cache },
			  0, -prm.ilg[0], -prm.ilg[1], -prm.ilg[2]);
}

CUDA_DEVICE static void
curr_cache_destroy(curr_cache_t curr_cache, fields_real_t *d_flds, int ci0[3])
{
}

