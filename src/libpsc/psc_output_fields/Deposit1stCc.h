
#pragma once

#include "psc_bits.h"
#include <psc/deposit.hxx>

#include <kg/Vec3.h>

// FIXME, this is still too intermingled, both doing the actual deposit as well
// as the particle / patch processing OTOH, there is still opportunity for
// optimization, in particular when operator() gets called back multiple times,
// we don't have to find the IP coefficients again Obviously, the rest of the IP
// macro should be converted, too

template <typename MF, typename D>
struct DepositContext
{
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;
  using real3_t = gt::sarray<real_t, 3>;

  DepositContext(Mfields& mflds, const real3_t& dx, real_t fnqs)
    : mflds{mflds}, ib{mflds.ib()}, deposit{dx, fnqs}
  {}

  void operator()(int m, real_t val)
  {
    deposit(mflds.storage().view(_all, _all, _all, m, p), ib, x, val);
  }

  Mfields& mflds;
  Int3 ib;
  psc::deposit::code::Deposit1stCc<real_t, dim_t> deposit;
  int p;
  real3_t x;
};

template <typename MF, typename D>
class Deposit1stCc
{
public:
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;
  using real3_t = gt::sarray<real_t, 3>;
  using DepositCtx = DepositContext<Mfields, dim_t>;

  template <typename Mparticles, typename F>
  void operator()(const real3_t& dx, real_t fnqs, Mfields& mflds,
                  const Mparticles& mprts, F&& func)
  {
    auto accessor = mprts.accessor();

    DepositCtx ctx{mflds, dx, fnqs};
    for (ctx.p = 0; ctx.p < mprts.n_patches(); ctx.p++) {
      for (auto prt : accessor[ctx.p]) {
        ctx.x = prt.x();
        func(ctx, prt);
      }
    }
  }
};

template <typename D, typename MF, typename MP, typename F>
void deposit1stCc(MF& mflds, const MP& mprts, F&& func)
{
  using real_t = typename MF::real_t;
  using real3_t = gt::sarray<real_t, 3>;

  const auto& domain = mflds.grid().domain;
  real3_t dx = {domain.dx[0], domain.dx[1], domain.dx[2]};
  real_t fnqs = mflds.grid().norm.fnqs;
  Deposit1stCc<MF, D> deposit;
  deposit(dx, fnqs, mflds, mprts, std::forward<F>(func));
}
