
#ifndef PUSH_DIM
#define PUSH_DIM DIM
#endif

#include <cmath>

#include "pushp.hxx"

// ----------------------------------------------------------------------
// calc_v

template<typename R>
CUDA_DEVICE static inline void
calc_v(R v[3], const R p[3])
{
#ifdef __CUDACC__
  R root = rsqrt(1.f + sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
#else
  R root = 1.f / std::sqrt(1.f + sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
#endif
  for (int d = 0; d < 3; d++) {
    v[d] = p[d] * root;
  }
}

// ----------------------------------------------------------------------
// push_x

template<typename R>
CUDA_DEVICE static inline void
push_x(R x[3], const R v[3], R dt)
{
  IF_DIM_X( x[0] += v[0] * dt; );
  IF_DIM_Y( x[1] += v[1] * dt; );
  IF_DIM_Z( x[2] += v[2] * dt; );
}

// ----------------------------------------------------------------------
// find_idx_off_1st_rel

template<typename R>
CUDA_DEVICE static inline void
find_idx_off_1st_rel(R xi[3], int lg[3], R og[3], R shift)
{
  for (int d = 0; d < 3; d++) {
    R pos = xi[d] * c_prm.dxi[d] + shift;
    lg[d] = fint(pos);
    og[d] = pos - lg[d];
  }
}

// ----------------------------------------------------------------------
// find_idx_off_pos_1st_rel

template<typename R>
CUDA_DEVICE static inline void
find_idx_off_pos_1st_rel(R xi[3], int lg[3], R og[3], R pos[3], R shift)
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * c_prm.dxi[d] + shift;
    lg[d] = fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

