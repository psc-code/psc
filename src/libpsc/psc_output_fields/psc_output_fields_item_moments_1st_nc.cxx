
#pragma once

#include <math.h>

#include "common_moments.cxx"
#include "Deposit1stCc.h"

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
    psc::moment::deposit_1st_nc<dim_t>(mflds.storage(), mflds.ib(), mprts,
                                       [&](auto& deposit_one, const auto& prt) {
                                         int m = prt.kind();
                                         deposit_one(m, prt.w());
                                       });
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
    Base::mres_.gt().view() = 0.f;
    psc::moment::deposit_1st_nc<dim_t>(
      Base::mres_.storage(), Base::mres_.ib(), mprts,
      [&](auto& deposit_one, const auto& prt) { deposit_one(0, prt.q()); });
    Base::bnd_.add_ghosts(Base::mres_);
  }
};
