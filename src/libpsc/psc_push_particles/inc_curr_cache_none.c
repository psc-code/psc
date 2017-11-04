
typedef flds_curr_t curr_cache_t;

#ifndef CURR_CACHE_DIM
#define CURR_CACHE_DIM DIM_XYZ
#endif

#if CURR_CACHE_DIM == DIM_XYZ
#define F3_CURR_CACHE(curr_cache, m, i,j,k)  (F3_CURR(curr_cache, m, i,j,k))
#elif CURR_CACHE_DIM == DIM_1
#define F3_CURR_CACHE(curr_cache, m, i,j,k)  (F3_CURR(curr_cache, m, 0,0,0))
#else
#error unhandled CURR_CACHE_DIM
#endif

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(flds_curr_t flds_curr, int ci0[3])
{
  return flds_curr;
}

CUDA_DEVICE static inline void
curr_cache_add(curr_cache_t curr_cache, int m, int i, int j, int k,
	       fields_real_t val)
{
  fields_real_t *addr = &F3_CURR_CACHE(curr_cache, JXI+m, i,j,k);
  atomicAdd(addr, val);
}

CUDA_DEVICE static inline void
curr_cache_destroy(curr_cache_t curr_cache, flds_curr_t flds_curr, int ci0[3])
{
}

