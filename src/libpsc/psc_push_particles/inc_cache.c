
// ======================================================================
// field caching

#ifndef EM_CACHE
#define EM_CACHE EM_CACHE_NONE
#endif

// OPT, shouldn't we be able to do with less ghosts?
#if DIM == DIM_YZ
#define BLOCKBND_X 0
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#elif DIM == DIM_XYZ
#define BLOCKBND_X 2
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#endif

#define BLOCKGSIZE_X (BLOCKSIZE_X + 2 * BLOCKBND_X)
#define BLOCKGSIZE_Y (BLOCKSIZE_Y + 2 * BLOCKBND_Y)
#define BLOCKGSIZE_Z (BLOCKSIZE_Z + 2 * BLOCKBND_Z)

// ----------------------------------------------------------------------

typedef fields_t em_cache_t;

#ifndef EM_CACHE_DIM
#define EM_CACHE_DIM DIM_XYZ
#endif

#if EM_CACHE_DIM == DIM_XYZ
using FieldsEM = Fields3d<em_cache_t, dim_xyz>;
#elif EM_CACHE_DIM == DIM_XZ
using FieldsEM = Fields3d<em_cache_t, dim_xz>;
#elif EM_CACHE_DIM == DIM_1
using FieldsEM = Fields3d<em_cache_t, dim_1>;
#else
#error unhandled EM_CACHE_DIM
#endif

CUDA_DEVICE static inline em_cache_t
em_cache_create(fields_t flds_em, int ci0[3])
{
  return flds_em;
}

