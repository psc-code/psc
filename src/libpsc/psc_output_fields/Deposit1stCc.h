
#pragma once

#include "psc_bits.h"
#include <psc/deposit.hxx>

#include <kg/Vec3.h>

// FIXME, this is still too intermingled, both doing the actual deposit as well
// as the particle / patch processing OTOH, there is still opportunity for
// optimization, in particular when operator() gets called back multiple times,
// we don't have to find the IP coefficients again Obviously, the rest of the IP
// macro should be converted, too

struct DepositContext
{
  int p;
};

template <typename T, typename D>
class Deposit1stCc
{
public:
  using dim_t = D;
  using real_t = T;

  Deposit1stCc(const Grid_t& grid)
    : deposit_({grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
               grid.norm.fnqs)

  {}

  template <typename MF, typename PRT>
  void operator()(DepositContext& ctx, MF& mflds, const PRT& prt, int m,
                  real_t val)
  {
    auto ib = mflds.ib();
    deposit_(mflds.storage().view(_all, _all, _all, m, ctx.p), ib, prt.x(),
             val);
  }

  template <typename Mparticles, typename F>
  void process(const Mparticles& mprts, F&& func)
  {
    auto accessor = mprts.accessor();

    DepositContext ctx;
    for (ctx.p = 0; ctx.p < mprts.n_patches(); ctx.p++) {
      for (auto prt : accessor[ctx.p]) {
        func(ctx, prt);
      }
    }
  }

private:
  psc::deposit::code::Deposit1stCc<real_t, dim_t> deposit_;
};
