
#pragma once

#include "Deposit1stCc.h"

#include "fields_item.hxx"

#include <cmath>

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

template <typename Particle>
static inline void __particle_calc_vxi(const Particle& prt,
                                       typename Particle::real_t vxi[3])
{
  typename Particle::real_t root =
    1.f / std::sqrt(1.f + sqr(prt.u[0]) + sqr(prt.u[1]) + sqr(prt.u[2]));
  vxi[0] = prt.u[0] * root;
  vxi[1] = prt.u[1] * root;
  vxi[2] = prt.u[2] * root;
}

// ======================================================================
// n_1st

template <typename MP, typename MF, typename D>
class Moment_n_1st : public ItemMomentCRTP<Moment_n_1st<MP, MF, D>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_n_1st<MP, MF, D>, MF>;
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mfields::real_t;
  using dim_t = D;

  using Base::n_comps;

  constexpr static char const* name = "n_1st";

  static int n_comps(const Grid_t& grid) { return 1 * grid.kinds.size(); }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  int n_comps() const { return Base::mres_.n_comps(); }
  Int3 ibn() const { return {}; }

  explicit Moment_n_1st(const Grid_t& grid) : Base{grid} {}

  explicit Moment_n_1st(const Mparticles& mprts) : Base{mprts.grid()}
  {
    update(mprts);
  }

  void update(const Mparticles& mprts)
  {
    Base::mres_.storage().view() = 0.f;
    psc::moment::deposit_1st_cc<dim_t>(Base::mres_.storage(), Base::mres_.ib(),
                                       mprts,
                                       [&](auto& deposit_one, const auto& prt) {
                                         int m = prt.kind();
                                         deposit_one(m, 1.f);
                                       });
    Base::bnd_.add_ghosts(Base::mres_);
  }

  auto gt()
  {
    Int3 bnd = Base::mres_.ibn();
    return Base::mres_.gt().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                                 _s(bnd[2], -bnd[2]));
  }
};

// ======================================================================
// v_1st

template <typename MF, typename D>
class Moment_v_1st : public ItemMomentCRTP<Moment_v_1st<MF, D>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_v_1st<MF, D>, MF>;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;

  constexpr static char const* name = "v_1st";

  static int n_comps(const Grid_t& grid) { return 3 * grid.kinds.size(); }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"vx", "vy", "vz"}, grid.kinds);
  }

  template <typename Mparticles>
  explicit Moment_v_1st(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_.storage().view() = 0.f;
    psc::moment::deposit_1st_cc<dim_t>(
      Base::mres_.storage(), Base::mres_.ib(), mprts,
      [&](auto& deposit_one, const auto& prt) {
        real_t vxi[3];
        _particle_calc_vxi(prt, vxi);
        for (int m = 0; m < 3; m++) {
          deposit_one(m + 3 * prt.kind(), vxi[m]);
        }
      });
    Base::bnd_.add_ghosts(Base::mres_);
  }

  auto gt()
  {
    Int3 bnd = Base::mres_.ibn();
    return Base::mres_.gt().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                                 _s(bnd[2], -bnd[2]));
  }
};

// ======================================================================
// p_1st

template <typename MF, typename D>
class Moment_p_1st : public ItemMomentCRTP<Moment_p_1st<MF, D>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_p_1st<MF, D>, MF>;
  using Mfields = MF;
  using dim_t = D;
  using real_t = typename Mfields::real_t;

  constexpr static char const* name = "p_1st";

  static int n_comps(const Grid_t& grid) { return 3 * grid.kinds.size(); }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"px", "py", "pz"}, grid.kinds);
  }

  template <typename Mparticles>
  explicit Moment_p_1st(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_.storage().view() = 0.f;
    psc::moment::deposit_1st_cc<dim_t>(
      Base::mres_.storage(), Base::mres_.ib(), mprts,
      [&](auto& deposit_one, const auto& prt) {
        for (int m = 0; m < 3; m++) {
          deposit_one(m + 3 * prt.kind(), prt.m() * prt.u()[m]);
        }
      });
    Base::bnd_.add_ghosts(Base::mres_);
  }

  auto gt()
  {
    Int3 bnd = Base::mres_.ibn();
    return Base::mres_.gt().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                                 _s(bnd[2], -bnd[2]));
  }
};

// ======================================================================
// T_1st

template <typename MF, typename D>
struct Moment_T_1st
{
  using Mfields = MF;
  using dim_t = D;

  constexpr static char const* name = "T_1st";

  static int n_comps(const Grid_t& grid) { return 3 * grid.kinds.size(); }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz"},
                         grid.kinds);
  }

  template <typename Mparticles>
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    using real_t = typename Mparticles::real_t;

    psc::moment::deposit_1st_cc<dim_t>(mflds, mprts, [&](const auto& prt) {
      int mm = prt.kind() * 6;

      real_t vxi[3];
      _particle_calc_vxi(prt, vxi);
      auto pxi = prt.u();
      deposit(mflds, prt, mm + 0, prt.m() * pxi[0] * vxi[0]);
      deposit(mflds, prt, mm + 1, prt.m() * pxi[1] * vxi[1]);
      deposit(mflds, prt, mm + 2, prt.m() * pxi[2] * vxi[2]);
      deposit(mflds, prt, mm + 3, prt.m() * pxi[0] * vxi[1]);
      deposit(mflds, prt, mm + 4, prt.m() * pxi[0] * vxi[2]);
      deposit(mflds, prt, mm + 5, prt.m() * pxi[1] * vxi[2]);
    });
  }
};

// ======================================================================
// Moments_1st
//
// all moments calculated at once
// FIXME: add KE

template <typename MP, typename MF, typename D>
class Moments_1st : public ItemMomentCRTP<Moments_1st<MP, MF, D>, MF>
{
public:
  using Base = ItemMomentCRTP<Moments_1st<MP, MF, D>, MF>;
  using Mparticles = MP;
  using Mfields = MF;
  using dim_t = D;
  using value_type = typename Mfields::real_t;
  using space = typename Mfields::space;

  using Base::n_comps;

  constexpr static int n_moments = 13;
  static char const* name() { return "all_1st"; }

  static int n_comps(const Grid_t& grid)
  {
    return n_moments * grid.kinds.size();
  }

  std::vector<std::string> comp_names()
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         Base::grid().kinds);
  }

  explicit Moments_1st(const Mparticles& mprts) : Base{mprts.grid()}
  {
    using real_t = typename Mparticles::real_t;

    psc::moment::deposit_1st_cc<dim_t>(
      Base::mres_.storage(), Base::mres_.ib(), mprts,
      [&](auto& deposit_one, const auto& prt) {
        int mm = prt.kind() * n_moments;
        real_t vxi[3];
        _particle_calc_vxi(prt, vxi);
        deposit_one(mm + 0, prt.q());
        deposit_one(mm + 1, prt.q() * vxi[0]);
        deposit_one(mm + 2, prt.q() * vxi[1]);
        deposit_one(mm + 3, prt.q() * vxi[2]);
        deposit_one(mm + 4, prt.m() * prt.u()[0]);
        deposit_one(mm + 5, prt.m() * prt.u()[1]);
        deposit_one(mm + 6, prt.m() * prt.u()[2]);
        deposit_one(mm + 7, prt.m() * prt.u()[0] * vxi[0]);
        deposit_one(mm + 8, prt.m() * prt.u()[1] * vxi[1]);
        deposit_one(mm + 9, prt.m() * prt.u()[2] * vxi[2]);
        deposit_one(mm + 10, prt.m() * prt.u()[0] * vxi[1]);
        deposit_one(mm + 11, prt.m() * prt.u()[1] * vxi[2]);
        deposit_one(mm + 12, prt.m() * prt.u()[2] * vxi[0]);
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

#ifdef USE_CUDA

#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "psc_particles_single.h"

// ======================================================================
// Moments_1st
//
// all moments calculated at once
// FIXME: add KE

template <typename BS, typename D>
class Moments_1st<MparticlesCuda<BS>, MfieldsSingle, D>
  : public ItemMomentCRTP<Moments_1st<MparticlesCuda<BS>, MfieldsSingle, D>,
                          MfieldsSingle>
{
public:
  using Base = ItemMomentCRTP<Moments_1st<MparticlesCuda<BS>, MfieldsSingle, D>,
                              MfieldsSingle>;
  using dim_t = D;
  using Mparticles = MparticlesCuda<BS>;
  using Mfields = MfieldsSingle;

  using Base::n_comps;

  using Sub = Moments_1st<MparticlesSingle, Mfields, D>;

  constexpr static int n_moments = Sub::n_moments;
  static char const* name() { return Sub::name(); }

  static int n_comps(const Grid_t& grid) { return Sub::n_comps(grid); }

  std::vector<std::string> comp_names()
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         Base::grid().kinds);
  }

  explicit Moments_1st(const Mparticles& _mprts) : Base{_mprts.grid()}
  {
    static int pr, pr_A, pr_B, pr_C;
    if (!pr) {
      pr = prof_register("Moments_1st cuda", 1., 0, 0);
      pr_A = prof_register("Moments_1st get", 1., 0, 0);
      pr_B = prof_register("Moments_1st process", 1., 0, 0);
      pr_C = prof_register("Moments_1st addg", 1., 0, 0);
    }

    prof_start(pr);
    prof_start(pr_A);
    auto& mprts = const_cast<Mparticles&>(_mprts);
    auto&& h_mprts = mprts.template get_as<MparticlesSingle>();
    prof_stop(pr_A);

    using Particle = typename MparticlesSingle::ConstAccessor::Particle;
    using Real = typename Particle::real_t;
    using R = Real;
    using real_t = R;

    psc::moment::deposit_1st_cc<dim_t>(
      Base::mres_.storage(), Base::mres_.ib(), h_mprts,
      [&](auto& deposit_one, const auto& prt) {
        int mm = prt.kind() * n_moments;
        real_t vxi[3];
        _particle_calc_vxi(prt, vxi);
        deposit_one(mm + 0, prt.q());
        deposit_one(mm + 1, prt.q() * vxi[0]);
        deposit_one(mm + 2, prt.q() * vxi[1]);
        deposit_one(mm + 3, prt.q() * vxi[2]);
        deposit_one(mm + 4, prt.m() * prt.u()[0]);
        deposit_one(mm + 5, prt.m() * prt.u()[1]);
        deposit_one(mm + 6, prt.m() * prt.u()[2]);
        deposit_one(mm + 7, prt.m() * prt.u()[0] * vxi[0]);
        deposit_one(mm + 8, prt.m() * prt.u()[1] * vxi[1]);
        deposit_one(mm + 9, prt.m() * prt.u()[2] * vxi[2]);
        deposit_one(mm + 10, prt.m() * prt.u()[0] * vxi[1]);
        deposit_one(mm + 11, prt.m() * prt.u()[1] * vxi[2]);
        deposit_one(mm + 12, prt.m() * prt.u()[2] * vxi[0]);
      });

    prof_stop(pr_B);

    prof_start(pr_C);
    Base::bnd_.add_ghosts(Base::mres_);
    prof_stop(pr_C);

    mprts.put_as(h_mprts, MP_DONT_COPY);
    prof_stop(pr);
  }

  auto gt()
  {
    auto bnd = -Base::mres_.ib();
    return Base::mres_.storage().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                                      _s(bnd[2], -bnd[2]));
  }
};

#endif
