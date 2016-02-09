
typedef flds_curr_t curr_cache_t;

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(flds_curr_t flds_curr, int ci0[3])
{
  return flds_curr;
}

CUDA_DEVICE static inline void
curr_cache_add(curr_cache_t curr_cache, int m, int jx, int jy, int jz,
	       fields_real_t val)
{
  fields_real_t *addr = &F3_CURR(curr_cache, JXI+m, jx,jy,jz);
  atomicAdd(addr, val);
}

CUDA_DEVICE static inline void
curr_cache_destroy(curr_cache_t curr_cache, flds_curr_t flds_curr, int ci0[3])
{
}

