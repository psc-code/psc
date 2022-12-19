
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

  constexpr static char const* name = "n_2nd_nc";

  static int n_comps(const Grid_t& grid) { return 1 * grid.kinds.size(); }

  std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  template <typename Mparticles>
  explicit Moment_n_2nd_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    psc::moment::deposit_2nd_nc<dim_t>(Base::mres_gt_, Base::mres_ib_, mprts,
                                       [&](auto& deposit_one, const auto& prt) {
                                         int m = prt.kind();
                                         deposit_one(m, prt.w());
                                       });
    Base::bnd_.add_ghosts(Base::mres_);
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

  static char const* name() { return "rho_2nd_nc"; }
  static int n_comps(const Grid_t& grid) { return 1; }
  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return {"rho"};
  }

  int n_comps() const { return Base::mres_.n_comps(); }

  template <typename Mparticles>
  explicit Moment_rho_2nd_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    psc::moment::deposit_2nd_nc<dim_t>(Base::mres_gt_, Base::mres_ib_, mprts,
                                       [&](auto& deposit_one, const auto& prt) {
                                         deposit_one(0, prt.w() * prt.q());
                                       });
    Base::bnd_.add_ghosts(Base::mres_);
  }

  auto gt()
  {
    auto bnd = -Base::mres_.ib();
    return Base::mres_.storage().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                                      _s(bnd[2], -bnd[2]));
  }
};
