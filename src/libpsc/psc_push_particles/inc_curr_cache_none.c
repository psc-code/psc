
#ifndef CURR_CACHE_DIM
#define CURR_CACHE_DIM DIM_XYZ
#endif

struct curr_cache_t : flds_curr_t
{
  using Self = curr_cache_t;
  
  curr_cache_t(const flds_curr_t& f)
    : flds_curr_t(f.ib(), f.im(), f.n_comps(), f.data_)
  {}
  
  void add(int m, int i, int j, int k, real_t val)
  {
#if CURR_CACHE_DIM == DIM_XYZ
    using FieldsJ = Fields3d<Self, dim_xyz>;
#elif CURR_CACHE_DIM == DIM_XZ
    using FieldsJ = Fields3d<Self, dim_xz>;
#elif CURR_CACHE_DIM == DIM_1
    using FieldsJ = Fields3d<Self, dim_1>;
#else
#error unhandled CURR_CACHE_DIM
#endif
    
    FieldsJ J(*this);
    real_t *addr = &J(JXI+m, i,j,k);
    atomicAdd(addr, val);
  }
};

