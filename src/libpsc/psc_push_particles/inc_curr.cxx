
#pragma once

#include <psc/gtensor.h>
#include <psc/current_deposition.hxx>

// ----------------------------------------------------------------------

template <typename _Order, typename _Dim, typename _fields_t>
struct CurrentNone;

// ----------------------------------------------------------------------

template <typename Curr, typename real_t = typename Curr::real_t>
GT_INLINE void deposit(Curr& curr, const int _i[3], const real_t fnqs[3],
                       real_t qni_wni, const real_t dx[3], const real_t xa[3],
                       real_t h, const int off[3])
{
  psc::CurrentDeposition1vb<Curr> deposition{{fnqs[0], fnqs[1], fnqs[2]}};

  int i[3];
  for (int d = 0; d < 3; d++) {
    i[d] = _i[d] + off[d];
  }
  deposition(curr, i, qni_wni, dx, xa, dim_xyz{});
}

#include "inc_curr_1vb_split.cxx"
#include "inc_curr_1vb_var1.cxx"
#include "inc_curr_1vb_2d.cxx"
#include "inc_curr_zigzag.cxx"

template <typename _Order, typename _Dim, typename _fields_t>
struct Current1vb;
