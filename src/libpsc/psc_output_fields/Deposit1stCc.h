
#pragma once

#include "psc_bits.h"
#include <psc/deposit.hxx>

#include <kg/Vec3.h>

// FIXME, this is still too intermingled, both doing the actual deposit as well
// as the particle / patch processing OTOH, there is still opportunity for
// optimization, in particular when operator() gets called back multiple times,
// we don't have to find the IP coefficients again Obviously, the rest of the IP
// macro should be converted, too

template <typename Mfields, typename D>
class Deposit1stCc
{
public:
  using dim_t = D;
  using FE = typename Mfields::fields_view_t;
  using R = typename Mfields::real_t;
  using real_t = typename Mfields::real_t;

  Deposit1stCc(Mfields& mflds)
    : mflds_{mflds},
      flds_{mflds[0]},
      deposit_({mflds.grid().domain.dx[0], mflds.grid().domain.dx[1],
                mflds.grid().domain.dx[2]},
               mflds.grid().norm.fnqs)

  {}

  template <typename PRT>
  void operator()(const PRT& prt, int m, R val)
  {
    auto ib = flds_.ib();
    deposit_(prt, flds_.storage().view(_all, _all, _all, m), ib, val);
  }

  template <typename Mparticles, typename F>
  void process(const Mparticles& mprts, F&& func)
  {
    auto accessor = mprts.accessor();

    for (int p = 0; p < mprts.n_patches(); p++) {
      flds_ = mflds_[p];
      for (auto prt : accessor[p]) {
        func(prt);
      }
    }
  }

private:
  Mfields& mflds_;
  FE flds_;
  psc::deposit::Deposit1stCc<real_t, dim_t> deposit_;
};
