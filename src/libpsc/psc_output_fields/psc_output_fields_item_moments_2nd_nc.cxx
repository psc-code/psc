
#include "fields.hxx"

#include <math.h>

#include "common_moments.cxx"

// ======================================================================
// n

template <typename MP, typename D, typename MF = Mfields<typename MP::real_t>>
class Moment_n_2nd_nc : public ItemMomentCRTP<Moment_n_2nd_nc<MP, D, MF>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_n_2nd_nc<MP, D, MF>, MF>;
  using Mparticles = MP;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mparticles::real_t;

  using Base::n_comps;

  constexpr static char const* name = "n_2nd_nc";

  static int n_comps(const Grid_t& grid) { return 1 * grid.kinds.size(); }

  std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  explicit Moment_n_2nd_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    const auto& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1],
           dzi = 1.f / grid.domain.dx[2];

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = Base::mres_[p];
      for (auto prt : accessor[p]) {
        int m = prt.kind();
        DEPOSIT_TO_GRID_2ND_NC(prt, flds, m, 1.f);
      }
    }
    Base::bnd_.add_ghosts(Base::mres_);
  }
};

// ======================================================================
// rho

template <typename MP, typename D, typename MF = Mfields<typename MP::real_t>>
class Moment_rho_2nd_nc
  : public ItemMomentCRTP<Moment_rho_2nd_nc<MP, D, MF>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_rho_2nd_nc<MP, D, MF>, MF>;
  using Mparticles = MP;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mparticles::real_t;

  static char const* name() { return "rho_2nd_nc"; }
  static int n_comps(const Grid_t& grid) { return 1; }
  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return {"rho"};
  }

  int n_comps() const { return Base::mres_.n_comps(); }

  explicit Moment_rho_2nd_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    const auto& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1],
           dzi = 1.f / grid.domain.dx[2];

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = Base::mres_[p];
      for (auto prt : accessor[p]) {
        int m = prt.kind();
        DEPOSIT_TO_GRID_2ND_NC(prt, flds, 0, prt.q());
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
