
#pragma once

#include "psc_bits.h"
#include <psc/deposit.hxx>

#include <kg/Vec3.h>

// FIXME, this is still too intermingled, both doing the actual deposit as well
// as the particle / patch processing OTOH, there is still opportunity for
// optimization, in particular when operator() gets called back multiple times,
// we don't have to find the IP coefficients again Obviously, the rest of the IP
// macro should be converted, too

template <typename T, typename D>
struct DepositContext
{
  using real_t = T;
  using dim_t = D;
  using real3_t = gt::sarray<real_t, 3>;

  DepositContext(const real3_t& dx, real_t fnqs) : deposit(dx, fnqs) {}

  int p;
  real3_t x;
  psc::deposit::code::Deposit1stCc<real_t, dim_t> deposit;
};

template <typename T, typename D>
class Deposit1stCc
{
public:
  using dim_t = D;
  using real_t = T;
  using real3_t = gt::sarray<real_t, 3>;
  using DepositCtx = DepositContext<real_t, dim_t>;

  Deposit1stCc(const Grid_t& grid)
    : dx_{grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
      fnqs_(grid.norm.fnqs)
  {}

  template <typename MF>
  void operator()(DepositCtx& ctx, MF& mflds, int m, real_t val)
  {
    auto ib = mflds.ib();
    ctx.deposit(mflds.storage().view(_all, _all, _all, m, ctx.p), ib, ctx.x,
                val);
  }

  template <typename Mparticles, typename F>
  void process(const Mparticles& mprts, F&& func)
  {
    auto accessor = mprts.accessor();

    DepositCtx ctx{dx_, fnqs_};
    for (ctx.p = 0; ctx.p < mprts.n_patches(); ctx.p++) {
      for (auto prt : accessor[ctx.p]) {
        ctx.x = prt.x();
        func(ctx, prt);
      }
    }
  }

private:
  real3_t dx_;
  real_t fnqs_;
};
