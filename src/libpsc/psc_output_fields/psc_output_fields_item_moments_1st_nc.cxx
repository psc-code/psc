
#pragma once

#include <math.h>

#include "common_moments.cxx"

template <typename R>
class Deposit1stNc
{
public:
  using real_t = R;
  using real3_t = gt::sarray<real_t, 3>;

  Deposit1stNc(const real3_t& dx, real_t fnqs)
    : dxi_{real_t(1.) / dx}, fnqs_{fnqs}
  {}

  template <typename P, typename F>
  void operator()(const P& prt, F& flds, int m, real_t val, const Grid_t& grid)
  {
    auto xi = prt.x(); /* don't shift back in time */
    Vec3<real_t> x = {xi[0] * dxi_[0], xi[1] * dxi_[1], xi[2] * dxi_[2]};
    real_t value = prt.w() * fnqs_ * val;

    auto fm = flds.storage().view(_all, _all, _all, m);
    auto ib = gt::shape(flds.ib()[0], flds.ib()[1], flds.ib()[2]);
    if (!grid.isInvar(0) && !grid.isInvar(1) && !grid.isInvar(2)) {
      psc::DepositNc<real_t, dim_xyz> deposit;
      deposit(fm, ib, x, value);
    } else if (grid.isInvar(0) && !grid.isInvar(1) && !grid.isInvar(2)) {
      psc::DepositNc<real_t, dim_yz> deposit;
      deposit(fm, ib, x, value);
    }
  }

  real3_t dxi_;
  real_t fnqs_;
};

// ======================================================================
// n

template <typename MP, typename MF>
struct Moment_n_1st_nc
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;

  constexpr static char const* name = "n_1st_nc";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"n"}; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    Deposit1stNc<real_t> deposit(
      {grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
      grid.norm.fnqs);

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt : accessor[p]) {
        int m = prt.kind();
        deposit(prt, flds, m, 1.f, grid);
      }
    }
  }
};

// ======================================================================
// rho

template <typename MP, typename MF>
struct Moment_rho_1st_nc : ItemMomentCRTP<Moment_rho_1st_nc<MP, MF>, MF>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc<MP, MF>, MF>;
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;

  using Base::n_comps;

  constexpr static char const* name = "rho_1st_nc";
  static int n_comps(const Grid_t& grid) { return 1; }
  static std::vector<std::string> fld_names() { return {"rho"}; }
  constexpr static int flags = 0;

  explicit Moment_rho_1st_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    const auto& grid = mprts.grid();
    Deposit1stNc<real_t> deposit(
      {grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
      grid.norm.fnqs);

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto res = Base::mres_[p];
      for (auto prt : accessor[p]) {
        deposit(prt, res, 0, prt.q(), grid);
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

// ======================================================================
// v

template <typename MP, typename MF>
struct Moment_v_1st_nc
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;

  constexpr static char const* name = "v_1st_nc";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"vx", "vy", "vz"}; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    Deposit1stNc<real_t> deposit(
      {grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
      grid.norm.fnqs);

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt : accessor[p]) {
        int mm = prt.kind() * 3;

        real_t vxi[3];
        particle_calc_vxi(prt, vxi);

        for (int m = 0; m < 3; m++) {
          deposit(prt, flds, mm + m, vxi[m], grid);
        }
      }
    }
  }
};

#define MAKE_POFI_OPS(MP, MF, TYPE)                                            \
  FieldsItemMomentOps<Moment_n_1st_nc<MP, MF>>                                 \
    psc_output_fields_item_n_1st_nc_##TYPE##_ops;                              \
  FieldsItemMomentOps<Moment_rho_1st_nc<MP, MF>>                               \
    psc_output_fields_item_rho_1st_nc_##TYPE##_ops;                            \
  FieldsItemMomentOps<Moment_v_1st_nc<MP, MF>>                                 \
    psc_output_fields_item_v_1st_nc_##TYPE##_ops;
