
#pragma once

#include "psc_bits.h"
#include <psc/deposit.hxx>

#include <kg/Vec3.h>

namespace psc
{
namespace moment
{

template <typename S, typename D>
struct DepositContext
{
  using storage_type = S;
  using dim_t = D;
  using real_t = typename storage_type::value_type;
  using real3_t = gt::sarray<real_t, 3>;

  DepositContext(storage_type& mflds_gt, const Int3& ib, const real3_t& dx,
                 real_t fnqs)
    : mflds_gt{mflds_gt}, ib{ib}, deposit{dx, fnqs}
  {}

  void operator()(int m, real_t val)
  {
    deposit(mflds_gt.view(_all, _all, _all, m, p), ib, x, val);
  }

  storage_type& mflds_gt;
  Int3 ib;
  psc::deposit::code::Deposit1stCc<real_t, dim_t> deposit;
  int p;
  real3_t x;
};

template <typename S, typename D>
class Deposit1stCc
{
public:
  using storage_type = S;
  using dim_t = D;
  using real_t = typename storage_type::value_type;
  using real3_t = gt::sarray<real_t, 3>;
  using DepositCtx = DepositContext<storage_type, dim_t>;

  template <typename Mparticles, typename F>
  void operator()(const real3_t& dx, real_t fnqs, storage_type& mflds_gt,
                  const Int3& ib, const Mparticles& mprts, F&& func)
  {
    auto accessor = mprts.accessor();

    DepositCtx ctx{mflds_gt, ib, dx, fnqs};
    for (ctx.p = 0; ctx.p < mprts.n_patches(); ctx.p++) {
      for (auto prt : accessor[ctx.p]) {
        ctx.x = prt.x();
        func(ctx, prt);
      }
    }
  }
};

template <typename D, typename MF, typename MP, typename F>
void deposit_1st_cc(MF& mflds, const MP& mprts, F&& func)
{
  using real_t = typename MF::real_t;
  using real3_t = gt::sarray<real_t, 3>;

  const auto& domain = mflds.grid().domain;
  real3_t dx = {domain.dx[0], domain.dx[1], domain.dx[2]};
  real_t fnqs = mflds.grid().norm.fnqs;
  Deposit1stCc<typename MF::Storage, D> deposit;
  deposit(dx, fnqs, mflds.storage(), mflds.ib(), mprts, std::forward<F>(func));
}

} // namespace moment
} // namespace psc
