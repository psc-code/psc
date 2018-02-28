
struct curr_cache_t
{
  using real_t = typename flds_curr_t::real_t;

  curr_cache_t(flds_curr_t& flds)
    : f_(flds)
  {}
  
  real_t* data()      { return f_.data(); }
  int n_comps() const { return f_.n_comps(); }
  Int3 ib() const     { return f_.ib(); }
  Int3 im() const     { return f_.im(); }

private:
  flds_curr_t f_;
};

#ifndef CURR_CACHE_DIM
#define CURR_CACHE_DIM DIM_XYZ
#endif

CUDA_DEVICE static inline curr_cache_t
curr_cache_create(flds_curr_t flds_curr, int ci0[3])
{
  return curr_cache_t{flds_curr};
}

CUDA_DEVICE static inline void
curr_cache_add(curr_cache_t curr_cache, int m, int i, int j, int k,
	       curr_cache_t::real_t val)
{
#if CURR_CACHE_DIM == DIM_XYZ
  using FieldsJ = Fields3d<curr_cache_t, dim_xyz>;
#elif CURR_CACHE_DIM == DIM_XZ
  using FieldsJ = Fields3d<curr_cache_t, dim_xz>;
#elif CURR_CACHE_DIM == DIM_1
  using FieldsJ = Fields3d<curr_cache_t, dim_1>;
#else
#error unhandled CURR_CACHE_DIM
#endif

  FieldsJ J(curr_cache);
  curr_cache_t::real_t *addr = &J(JXI+m, i,j,k);
  atomicAdd(addr, val);
}

CUDA_DEVICE static inline void
curr_cache_destroy(curr_cache_t curr_cache, flds_curr_t flds_curr, int ci0[3])
{
}

