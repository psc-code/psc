
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
    int jx = nint(u);                                                          \
    int jy = nint(v);                                                          \
    int jz = nint(w);                                                          \
    real_t h1 = jx - u;                                                        \
    real_t h2 = jy - v;                                                        \
    real_t h3 = jz - w;                                                        \
    psc::deposit::norm::Deposit2nd<real_t, dim_t> deposit;                     \
    assert(jx >= 0 && jx <= grid.ldims[0]);                                    \
    assert(jy >= 0 && jy <= grid.ldims[1]);                                    \
    assert(jz >= 0 && jz <= grid.ldims[2]);                                    \
                                                                               \
    real_t fnq = prt.w() * fnqs;                                               \
                                                                               \
    auto fld = flds.storage().view(_all, _all, _all, m);                       \
    deposit(fld, {jx - flds.ib()[0], jy - flds.ib()[1], jz - flds.ib()[2]},    \
            {h1, h2, h3}, fnq* val);                                           \
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
