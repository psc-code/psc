
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

// ===========================================================================
// FIXME _particle_calc_vxi

template <typename Particle>
static inline void _particle_calc_vxi(const Particle& prt,
                                      typename Particle::real_t vxi[3])
{
  typename Particle::real_t root =
    1.f / std::sqrt(1.f + sqr(prt.u()[0]) + sqr(prt.u()[1]) + sqr(prt.u()[2]));
  vxi[0] = prt.u()[0] * root;
  vxi[1] = prt.u()[1] * root;
  vxi[2] = prt.u()[2] * root;
}

// ===========================================================================
// moment_n

template <template <typename, typename> class DepositCode, typename D>
class moment_n
{
public:
  using dim_t = D;

  static std::string name() { return "n" + DepositCode<float, D>::suffix(); }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return addKindSuffix({"n"}, kinds);
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    deposit<DepositCode, dim_t>(mflds_gt, ib, mprts,
                                [&](auto& deposit_one, const auto& prt) {
                                  int m = prt.kind();
                                  deposit_one(m, prt.w());
                                });
  }
};

// ===========================================================================
// moment_rho

template <template <typename, typename> class DepositCode, typename D>
class moment_rho
{
public:
  using dim_t = D;

  static std::string name() { return "rho" + DepositCode<float, D>::suffix(); }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return {"rho"};
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    deposit<DepositCode, dim_t>(mflds_gt, ib, mprts,
                                [&](auto& deposit_one, const auto& prt) {
                                  deposit_one(0, prt.w() * prt.q());
                                });
  }
};

// ===========================================================================
// moment_v

template <template <typename, typename> class DepositCode, typename D>
class moment_v
{
public:
  using dim_t = D;

  static std::string name() { return "v" + DepositCode<float, D>::suffix(); }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return addKindSuffix({"vx", "vy", "vz"}, kinds);
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    using real_t = typename MP::real_t;
    deposit<DepositCode, dim_t>(
      mflds_gt, ib, mprts, [&](auto& deposit_one, const auto& prt) {
        real_t vxi[3];
        _particle_calc_vxi(prt, vxi);
        for (int m = 0; m < 3; m++) {
          deposit_one(m + 3 * prt.kind(), prt.w() * vxi[m]);
        }
      });
  }
};

// ===========================================================================
// moment_p

template <template <typename, typename> class DepositCode, typename D>
class moment_p
{
public:
  using dim_t = D;

  static std::string name() { return "p" + DepositCode<float, D>::suffix(); }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return addKindSuffix({"px", "py", "pz"}, kinds);
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    deposit<DepositCode, dim_t>(
      mflds_gt, ib, mprts, [&](auto& deposit_one, const auto& prt) {
        for (int m = 0; m < 3; m++) {
          deposit_one(m + 3 * prt.kind(), prt.w() * prt.m() * prt.u()[m]);
        }
      });
  }
};

// ===========================================================================
// moment_T

template <template <typename, typename> class DepositCode, typename D>
class moment_T
{
public:
  using dim_t = D;

  static std::string name() { return "T" + DepositCode<float, D>::suffix(); }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return addKindSuffix({"Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz"}, kinds);
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    using real_t = typename MP::real_t;
    deposit<DepositCode, dim_t>(mflds_gt, ib, mprts, [&](const auto& prt) {
      int mm = prt.kind() * 6;
      real_t vxi[3];
      _particle_calc_vxi(prt, vxi);
      deposit_one(mm + 0, prt.w() * prt.m() * prt.u()[0] * vxi[0]);
      deposit_one(mm + 1, prt.w() * prt.m() * prt.u()[1] * vxi[1]);
      deposit_one(mm + 2, prt.w() * prt.m() * prt.u()[2] * vxi[2]);
      deposit_one(mm + 3, prt.w() * prt.m() * prt.u()[0] * vxi[1]);
      deposit_one(mm + 4, prt.w() * prt.m() * prt.u()[0] * vxi[2]);
      deposit_one(mm + 5, prt.w() * prt.m() * prt.u()[1] * vxi[2]);
    });
  }
}; // namespace moment

// ===========================================================================
// moment_all

template <template <typename, typename> class DepositCode, typename D>
class moment_all
{
public:
  using dim_t = D;

  static std::string name() { return "all" + DepositCode<float, D>::suffix(); }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         kinds);
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    using real_t = typename MP::real_t;
    deposit<DepositCode, dim_t>(
      mflds_gt, ib, mprts, [&](auto& deposit_one, const auto& prt) {
        int mm = prt.kind() * 13;
        real_t vxi[3];
        _particle_calc_vxi(prt, vxi);
        deposit_one(mm + 0, prt.w() * prt.q());
        deposit_one(mm + 1, prt.w() * prt.q() * vxi[0]);
        deposit_one(mm + 2, prt.w() * prt.q() * vxi[1]);
        deposit_one(mm + 3, prt.w() * prt.q() * vxi[2]);
        deposit_one(mm + 4, prt.w() * prt.m() * prt.u()[0]);
        deposit_one(mm + 5, prt.w() * prt.m() * prt.u()[1]);
        deposit_one(mm + 6, prt.w() * prt.m() * prt.u()[2]);
        deposit_one(mm + 7, prt.w() * prt.m() * prt.u()[0] * vxi[0]);
        deposit_one(mm + 8, prt.w() * prt.m() * prt.u()[1] * vxi[1]);
        deposit_one(mm + 9, prt.w() * prt.m() * prt.u()[2] * vxi[2]);
        deposit_one(mm + 10, prt.w() * prt.m() * prt.u()[0] * vxi[1]);
        deposit_one(mm + 11, prt.w() * prt.m() * prt.u()[1] * vxi[2]);
        deposit_one(mm + 12, prt.w() * prt.m() * prt.u()[2] * vxi[0]);
      });
  }
};

} // namespace moment
} // namespace psc
