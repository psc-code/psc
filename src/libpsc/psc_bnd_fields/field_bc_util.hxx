#pragma once

#include <limits>

#include "kg/VecRange.hxx"

namespace psc
{
namespace bnd
{
namespace field
{
namespace detail
{

/**
 * @brief Set E or B lower ghosts to the given constants (each component has
 * its own constant).
 * @param mflds mflds
 * @param p patch index
 * @param d which dimension to set the ghosts of
 * @param mb `EX` or `HX`; note that `mb+1` and `mb+2` are also set
 * @param val the constants
 * @param include_edge whether or not values located on exact domain edges
 * should be considered "ghosts"
 */
template <typename dim_t, typename MfieldsState>
void set_lower_ghosts(MfieldsState& mflds, int p, int d, int mb,
                      typename MfieldsState::Real3 val, bool include_edge)
{
  auto F = make_Fields3d<dim_t>(mflds[p]);
  Int3 start = mflds.ib();
  Int3 stop = mflds.ib() + mflds.im();
  stop[d] = 0;

  // TODO use gtensor views instead of VecRange

  for (int m = mb; m < mb + 3; m++) {
    for (Int3 i3 : VecRange(start, stop)) {
      F(m, i3) = val[m - mb];
    }
  }

  if (!include_edge) {
    return;
  }

  Int3 edge_start = mflds.ib();
  Int3 edge_stop = mflds.ib() + mflds.im();
  edge_start[d] = 0;
  edge_stop[d] = 1;

  for (int m = mb; m < mb + 3; m++) {
    bool edge_ec = mb == EX && m - mb != d;
    bool edge_fc = mb == HX && m - mb == d;

    if (edge_ec || edge_fc) {
      for (Int3 i3 : VecRange(edge_start, edge_stop)) {
        F(m, i3) = val[m - mb];
      }
    }
  }
}

/**
 * @brief Set E or B upper ghosts to the given constants (each component has
 * its own constant).
 * @param mflds mflds
 * @param p patch index
 * @param d which dimension to set the ghosts of
 * @param mb `EX` or `HX`; note that `mb+1` and `mb+2` are also set
 * @param val the constants
 * @param include_edge whether or not values located on exact domain edges
 * should be considered "ghosts"
 */
template <typename dim_t, typename MfieldsState>
void set_upper_ghosts(MfieldsState& mflds, int p, int d, int mb,
                      typename MfieldsState::Real3 val, bool include_edge)
{
  auto F = make_Fields3d<dim_t>(mflds[p]);
  Int3 start = mflds.ib();
  Int3 stop = mflds.ib() + mflds.im();
  start[d] = mflds.grid().ldims[d] + 1;

  // TODO use gtensor views instead of VecRange

  for (int m = mb; m < mb + 3; m++) {
    for (Int3 i3 : VecRange(start, stop)) {
      F(m, i3) = val[m - mb];
    }
  }

  Int3 edge_start = mflds.ib();
  Int3 edge_stop = mflds.ib() + mflds.im();
  edge_start[d] = mflds.grid().ldims[d];
  edge_stop[d] = mflds.grid().ldims[d] + 1;

  for (int m = mb; m < mb + 3; m++) {
    bool not_edge_ec = mb == EX && m - mb == d;
    bool not_edge_fc = mb == HX && m - mb != d;

    if (not_edge_ec || not_edge_fc || include_edge)
      for (Int3 i3 : VecRange(edge_start, edge_stop)) {
        F(m, i3) = val[m - mb];
      }
  }
}

template <typename dim_t, typename MfieldsState>
void set_lower_ghosts_to_nan(MfieldsState& mflds, int p, int d, int mb,
                             bool include_edge)
{
#ifndef DEBUG
  return;
#endif
  using real_t = typename MfieldsState::real_t;
  real_t nan = std::numeric_limits<real_t>::quiet_NaN();
  set_lower_ghosts<dim_t>(mflds, p, d, mb, {nan, nan, nan}, include_edge);
}

template <typename dim_t, typename MfieldsState>
void set_upper_ghosts_to_nan(MfieldsState& mflds, int p, int d, int mb,
                             bool include_edge)
{
#ifndef DEBUG
  return;
#endif
  using real_t = typename MfieldsState::real_t;
  real_t nan = std::numeric_limits<real_t>::quiet_NaN();
  set_upper_ghosts<dim_t>(mflds, p, d, mb, {nan, nan, nan}, include_edge);
}

} // namespace detail
} // namespace field
} // namespace bnd
} // namespace psc