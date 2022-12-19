
#pragma once

#include <math.h>

#include <psc/moment.hxx>

// ======================================================================
// n

template <typename MF, typename D>
struct Moment_n_1st_nc
{
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;

  constexpr static char const* name = "n_1st_nc";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"n"}; }
  constexpr static int flags = POFI_BY_KIND;

  template <typename Mparticles>
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

template <typename MF, typename D>
struct Moment_rho_1st_nc : ItemMomentCRTP<Moment_rho_1st_nc<MF, D>, MF>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc<MF, D>, MF>;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;

  constexpr static char const* name = "rho_1st_nc";
  static int n_comps(const Grid_t& grid) { return 1; }
  static std::vector<std::string> fld_names() { return {"rho"}; }
  constexpr static int flags = 0;

  template <typename Mparticles>
  explicit Moment_rho_1st_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    psc::moment::deposit_1st_nc<dim_t>(
      Base::mres_gt_, Base::mres_ib_, mprts,
      [&](auto& deposit_one, const auto& prt) { deposit_one(0, prt.q()); });
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
  }
};
