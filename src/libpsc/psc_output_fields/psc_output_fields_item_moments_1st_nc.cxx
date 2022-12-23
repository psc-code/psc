
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
  using moment_type =
    psc::moment::moment_n<psc::deposit::code::Deposit1stNc, dim_t>;

  template <typename Mparticles>
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    moment_type{}(mflds.storage(), mflds.ib(), mprts);
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
  using storage_type = typename Base::storage_type;
  using real_t = typename Base::real_t;
  using space_type = typename Base::space_type;
  using moment_type =
    psc::moment::moment_rho<psc::deposit::code::Deposit1stNc, dim_t>;

  explicit Moment_rho_1st_nc(const Grid_t& grid) : Base{grid} {}

  template <typename Mparticles>
  auto operator()(const Mparticles& mprts)
  {
    Int3 ib = -mprts.grid().ibn;
    storage_type mres =
      psc::mflds::zeros<real_t, space_type>(mprts.grid(), 1, ib);
    moment_type{}(mres, ib, mprts);
    Base::bnd_.add_ghosts(mprts.grid(), mres, ib);
    return mres;
  }

  auto storage() = delete;
};
