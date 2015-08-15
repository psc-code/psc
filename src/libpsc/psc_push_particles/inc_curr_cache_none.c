

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
  atomicAdd(addr, val);
}

CUDA_DEVICE static inline flds_curr_t
flds_curr_shift(flds_curr_t flds_curr, int m, int dx, int dy, int dz)
{
  return (flds_curr_t) { .arr_shift = flds_curr.arr_shift
      + F3_DEV_SHIFT_OFF(m, dx, dy, dz) };
}

CUDA_DEVICE static inline flds_curr_t
curr_cache_create(float *flds_curr, int ci0[3])
{
  return flds_curr_shift((flds_curr_t) { .arr_shift = flds_curr },
			 0, -prm.ilg[0], -prm.ilg[1], -prm.ilg[2]);
}

CUDA_DEVICE static void
curr_cache_add(flds_curr_t flds_curr, fields_real_t *d_flds, int ci0[3])
{
}

