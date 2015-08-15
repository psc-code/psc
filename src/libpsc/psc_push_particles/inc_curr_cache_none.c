
typedef flds_curr_t curr_cache_t;

#ifdef PSC_FIELDS_AS_CUDA2
#define CURR_CACHE_HAVE_SHIFT
#endif

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

CUDA_DEVICE static inline curr_cache_t
curr_cache_shift(curr_cache_t curr_cache, int m, int dx, int dy, int dz)
{
  return curr_cache + F3_DEV_SHIFT_OFF(m, dx, dy, dz);
}

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(flds_curr_t flds_curr, int ci0[3])
{
  return flds_curr;
}

CUDA_DEVICE static inline void
curr_cache_add(curr_cache_t curr_cache, int m, int jx, int jy, int jz,
	       fields_real_t val)
{
#if PSC_FIELDS_AS_CUDA2
  fields_real_t *addr = &F3_DEV(curr_cache, JXI+m, jx,jy,jz);
#else
  fields_real_t *addr = &F3_CURR(curr_cache, JXI+m, jx,jy,jz);
#endif
  atomicAdd(addr, val);
}

CUDA_DEVICE static inline void
curr_cache_destroy(curr_cache_t curr_cache, flds_curr_t flds_curr, int ci0[3])
{
}

