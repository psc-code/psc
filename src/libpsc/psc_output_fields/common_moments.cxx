
#pragma once

#include "fields_item.hxx"

#include <string>

#include <psc/deposit.hxx>

#define DEPOSIT_TO_GRID_2ND_NC(prt, flds, m, val)                              \
  do {                                                                         \
    auto xi = prt.x(); /* don't shift back in time */                          \
    psc::deposit::code::Deposit2ndNc<real_t, dim_t> deposit(                   \
      {grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},               \
      grid.norm.fnqs);                                                         \
                                                                               \
    auto fld = flds.storage().view(_all, _all, _all, m);                       \
    deposit(fld, flds.ib(), {xi[0], xi[1], xi[2]}, prt.w());                   \
  } while (0)

// FIXME, this function exists about 100x all over the place, should
// be consolidated

template <typename Particle>
static inline void particle_calc_vxi(const Particle& prt,
                                     typename Particle::real_t vxi[3])
{
  typename Particle::real_t root =
    1.f / std::sqrt(1.f + sqr(prt.u()[0]) + sqr(prt.u()[1]) + sqr(prt.u()[2]));
  vxi[0] = prt.u()[0] * root;
  vxi[1] = prt.u()[1] * root;
  vxi[2] = prt.u()[2] * root;
}
