
#pragma once

#include "fields_item.hxx"

#include <string>

#include <psc/deposit.hxx>

#define DEPOSIT_TO_GRID_2ND_NC(prt, flds, m, val)                              \
  do {                                                                         \
    auto xi = prt.x(); /* don't shift back in time */                          \
    real_t u = xi[0] * dxi;                                                    \
    real_t v = xi[1] * dyi;                                                    \
    real_t w = xi[2] * dzi;                                                    \
    psc::deposit::norm::Deposit2ndNc<real_t, dim_t> deposit;                   \
    auto fld = flds.storage().view(_all, _all, _all, m);                       \
    real_t fnq = prt.w() * fnqs;                                               \
    deposit(fld, flds.ib(), {u, v, w}, fnq* val);                              \
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
