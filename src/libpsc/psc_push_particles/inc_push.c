
#include <cmath>

#include "inc_params.c"
#include "pushp.hxx"

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

