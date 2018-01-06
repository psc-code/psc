
typedef flds_curr_t curr_cache_t;

#ifndef CURR_CACHE_DIM
#define CURR_CACHE_DIM DIM_XYZ
#endif

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(flds_curr_t flds_curr, int ci0[3])
{
  return flds_curr;
}

CUDA_DEVICE static inline void
curr_cache_add(curr_cache_t curr_cache, int m, int i, int j, int k,
	       fields_t::real_t val)
{
#if CURR_CACHE_DIM == DIM_XYZ
  using FieldsJ = Fields3d<flds_curr_t, dim_xyz>;
#elif CURR_CACHE_DIM == DIM_XZ
  using FieldsJ = Fields3d<flds_curr_t, dim_xz>;
#elif CURR_CACHE_DIM == DIM_1
  using FieldsJ = Fields3d<flds_curr_t, dim_1>;
#else
#error unhandled CURR_CACHE_DIM
#endif

  FieldsJ J(curr_cache);
  fields_t::real_t *addr = &J(JXI+m, i,j,k);
  atomicAdd(addr, val);
}

CUDA_DEVICE static inline void
curr_cache_destroy(curr_cache_t curr_cache, flds_curr_t flds_curr, int ci0[3])
{
}

