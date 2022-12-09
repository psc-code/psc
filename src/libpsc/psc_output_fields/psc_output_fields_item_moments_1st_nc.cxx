
#pragma once

#include <math.h>

#include "common_moments.cxx"

template <typename R, typename D>
class Deposit1stNc
{
public:
  using real_t = R;
  using dim_t = D;
  using real3_t = gt::sarray<real_t, 3>;

  Deposit1stNc(const real3_t& dx, real_t fnqs)
    : dxi_{real_t(1.) / dx}, fnqs_{fnqs}
  {}

  template <typename P, typename F>
  void operator()(const P& prt, const F& flds, const gt::shape_type<3>& ib,
                  real_t val)
  {
    auto xi = prt.x(); /* don't shift back in time */
    Vec3<real_t> x = {xi[0] * dxi_[0], xi[1] * dxi_[1], xi[2] * dxi_[2]};
    real_t value = fnqs_ * val;

    psc::DepositNc<real_t, dim_t> deposit;
    deposit(flds, ib, x, value);
  }

  real3_t dxi_;
  real_t fnqs_;
};

template <typename R, typename D>
class Moments_n
{
public:
  using real_t = R;
  using dim_t = D;

  Moments_n(const Grid_t& grid)
    : deposit_({grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
               grid.norm.fnqs)
  {}

  template <typename MF, typename MP>
  void operator()(MF& mflds, MP& mprts)
  {
    auto ib = gt::shape(mflds.ib()[0], mflds.ib()[1], mflds.ib()[2]);
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds.gt().view(_all, _all, _all, _all, p);
      for (auto prt : accessor[p]) {
        auto fld = flds.view(_all, _all, _all, prt.kind());
        deposit_(prt, fld, ib, prt.w());
      }
    }
  }

  Deposit1stNc<real_t, dim_t> deposit_;
};

template <typename R, typename D>
class Moments_rho
{
public:
  using real_t = R;
  using dim_t = D;

  Moments_rho(const Grid_t& grid)
    : deposit_({grid.domain.dx[0], grid.domain.dx[1], grid.domain.dx[2]},
               grid.norm.fnqs)
  {}

  template <typename MF, typename MP>
  void operator()(MF& mflds, MP& mprts)
  {
    auto ib = gt::shape(mflds.ib()[0], mflds.ib()[1], mflds.ib()[2]);
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto fld = mflds.gt().view(_all, _all, _all, 0, p);
      for (auto prt : accessor[p]) {
        deposit_(prt, fld, ib, prt.q());
      }
    }
  }

  Deposit1stNc<real_t, dim_t> deposit_;
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
    const Grid_t& grid = mprts.grid();
    if (!grid.isInvar(0) && !grid.isInvar(1) && !grid.isInvar(2)) {
      Moments_n<real_t, dim_xyz> moments{grid};
      moments(mflds, mprts);
    } else if (grid.isInvar(0) && !grid.isInvar(1) && !grid.isInvar(2)) {
      Moments_n<real_t, dim_yz> moments{grid};
      moments(mflds, mprts);
    } else {
      assert(0);
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
    const Grid_t& grid = mprts.grid();
    if (!grid.isInvar(0) && !grid.isInvar(1) && !grid.isInvar(2)) {
      Moments_rho<real_t, dim_xyz> moments{grid};
      moments(Base::mres_, mprts);
    } else if (grid.isInvar(0) && !grid.isInvar(1) && !grid.isInvar(2)) {
      Moments_rho<real_t, dim_yz> moments{grid};
      moments(Base::mres_, mprts);
    } else {
      assert(0);
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

#define MAKE_POFI_OPS(MP, MF, TYPE)                                            \
  FieldsItemMomentOps<Moment_n_1st_nc<MP, MF>>                                 \
    psc_output_fields_item_n_1st_nc_##TYPE##_ops;                              \
  FieldsItemMomentOps<Moment_rho_1st_nc<MP, MF>>                               \
    psc_output_fields_item_rho_1st_nc_##TYPE##_ops;
