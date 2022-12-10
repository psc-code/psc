
#pragma once

#include <math.h>

#include "common_moments.cxx"

// ======================================================================
// n

template <typename MP, typename MF, typename D>
struct Moment_n_1st_nc
{
  using Mparticles = MP;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;

  constexpr static char const* name = "n_1st_nc";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"n"}; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    psc::deposit::Deposit1stNc<real_t, dim_t> deposit(
      {grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
      grid.norm.fnqs);

    auto ib = mflds.ib();
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds.gt().view(_all, _all, _all, _all, p);
      for (auto prt : accessor[p]) {
        auto fld = flds.view(_all, _all, _all, prt.kind());
        deposit(prt, fld, ib, prt.w());
      }
    }
  }
};

// ======================================================================
// rho

template <typename MP, typename MF, typename D>
struct Moment_rho_1st_nc : ItemMomentCRTP<Moment_rho_1st_nc<MP, MF, D>, MF>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc<MP, MF, D>, MF>;
  using Mparticles = MP;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;

  using Base::n_comps;

  constexpr static char const* name = "rho_1st_nc";
  static int n_comps(const Grid_t& grid) { return 1; }
  static std::vector<std::string> fld_names() { return {"rho"}; }
  constexpr static int flags = 0;

  explicit Moment_rho_1st_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    const Grid_t& grid = mprts.grid();
    psc::deposit::Deposit1stNc<real_t, dim_t> deposit(
      {grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
      grid.norm.fnqs);

    auto ib = Base::mres_.ib();
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto fld = Base::mres_.gt().view(_all, _all, _all, 0, p);
      for (auto prt : accessor[p]) {
        deposit(prt, fld, ib, prt.q());
      }
    }
    Base::bnd_.add_ghosts(Base::mres_);
  }

  auto gt()
  {
    auto bnd = -Base::mres_.ib();
    return Base::mres_.storage().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                                      _s(bnd[2], -bnd[2]));
  }
};
