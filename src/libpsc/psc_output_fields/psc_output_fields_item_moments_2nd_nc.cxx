
#include "fields.hxx"
#include <psc/moment.hxx>

#include <math.h>

// ======================================================================
// n

template <typename MF, typename D>
class Moment_n_2nd_nc : public ItemMomentCRTP<Moment_n_2nd_nc<MF, D>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_n_2nd_nc<MF, D>, MF>;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;
  using moment_type =
    psc::moment::moment_n<psc::deposit::code::Deposit2ndNc, dim_t>;

  static std::string name_impl() { return "n_2nd_nc"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  template <typename Mparticles>
  explicit Moment_n_2nd_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    moment_type{}(Base::mres_gt_, Base::mres_ib_, mprts);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
  }
};

// ======================================================================
// rho

template <typename MF, typename D>
class Moment_rho_2nd_nc : public ItemMomentCRTP<Moment_rho_2nd_nc<MF, D>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_rho_2nd_nc<MF, D>, MF>;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;
  using moment_type =
    psc::moment::moment_rho<psc::deposit::code::Deposit2ndNc, dim_t>;

  static std::string name_impl() { return "rho_2nd_nc"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return {"rho"};
  }

  template <typename Mparticles>
  explicit Moment_rho_2nd_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    moment_type{}(Base::mres_gt_, Base::mres_ib_, mprts);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
  }
};
