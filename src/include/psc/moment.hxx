
#pragma once

#include "psc_bits.h"
#include <psc/deposit.hxx>

namespace psc
{
namespace moment
{

template <typename S, typename D,
          template <typename, typename> class DepositCode>
struct DepositParticlesContext
{
  using storage_type = S;
  using dim_t = D;
  using real_t = typename storage_type::value_type;
  using real3_t = gt::sarray<real_t, 3>;
  using Deposit = DepositCode<real_t, dim_t>;

  DepositParticlesContext(storage_type& mflds_gt, const Int3& ib,
                          const real3_t& dx, real_t fnqs)
    : mflds_gt{mflds_gt}, ib{ib}, deposit{dx, fnqs}
  {}

  void operator()(int m, real_t val)
  {
    deposit(mflds_gt.view(_all, _all, _all, m, p), ib, x, val);
  }

  storage_type& mflds_gt;
  Int3 ib;
  Deposit deposit;
  int p;
  real3_t x;
};

template <typename S, typename D,
          template <typename, typename> class DepositCode>
class DepositParticles
{
public:
  using storage_type = S;
  using dim_t = D;
  using real_t = typename storage_type::value_type;
  using real3_t = gt::sarray<real_t, 3>;
  using DepositCtx = DepositParticlesContext<storage_type, dim_t, DepositCode>;

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

template <template <typename, typename> class DepositCode, typename D,
          typename MF, typename MP, typename F>
void deposit(MF& mflds_gt, const Int3& ib, const MP& mprts, F&& func)
{
  using real_t = typename MF::value_type;
  using real3_t = gt::sarray<real_t, 3>;

  const auto& domain = mprts.grid().domain;
  real3_t dx = {domain.dx[0], domain.dx[1], domain.dx[2]};
  real_t fnqs = mprts.grid().norm.fnqs;
  DepositParticles<MF, D, DepositCode> deposit;
  deposit(dx, fnqs, mflds_gt, ib, mprts, std::forward<F>(func));
}

template <typename D, typename MF, typename MP, typename F>
void deposit_1st_cc(MF& mflds_gt, const Int3& ib, const MP& mprts, F&& func)
{
  deposit<psc::deposit::code::Deposit1stCc, D>(mflds_gt, ib, mprts,
                                               std::forward<F>(func));
}

template <typename D, typename MF, typename MP, typename F>
void deposit_1st_nc(MF& mflds_gt, const Int3& ib, const MP& mprts, F&& func)
{
  deposit<psc::deposit::code::Deposit1stNc, D>(mflds_gt, ib, mprts,
                                               std::forward<F>(func));
}

template <typename D, typename MF, typename MP, typename F>
void deposit_2nd_nc(MF& mflds_gt, const Int3& ib, const MP& mprts, F&& func)
{
  deposit<psc::deposit::code::Deposit2ndNc, D>(mflds_gt, ib, mprts,
                                               std::forward<F>(func));
}

} // namespace moment
} // namespace psc
