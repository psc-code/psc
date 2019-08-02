
#pragma once

#include "Deposit1stCc.h"

#include "fields_item.hxx"

#include <cmath>

// ======================================================================
// n_1st

template <typename MF>
struct Moment_n_1st : ItemMomentCRTP<Moment_n_1st<MF>, MF>
{
  using Base = ItemMomentCRTP<Moment_n_1st<MF>, MF>;
  using Mfields = MF;

  constexpr static char const* name = "n_1st";

  static int n_comps(const Grid_t& grid) { return 1 * grid.kinds.size(); }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  using Base::Base;

  template <typename Mparticles>
  void operator()(Mparticles& mprts)
  {
    using Particle = typename Mparticles::ConstAccessor::Particle;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, Base::mres_};
    deposit.process([&](const Particle& prt) {
      int m = prt.kind();
      deposit(prt, m, 1.f);
    });
    Base::bnd_.add_ghosts(Base::mres_);
  }
};

// ======================================================================
// v_1st

template <typename MF>
struct Moment_v_1st
{
  using Mfields = MF;

  constexpr static char const* name = "v_1st";

  static int n_comps(const Grid_t& grid) { return 3 * grid.kinds.size(); }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"vx", "vy", "vz"}, grid.kinds);
  }

  template <typename Mparticles>
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    using Particle = typename Mparticles::ConstAccessor::Particle;
    using Real = typename Particle::real_t;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      Real vxi[3];
      particle_calc_vxi(prt, vxi);

      int mm = prt.kind() * 3;
      for (int m = 0; m < 3; m++) {
        deposit(prt, mm + m, vxi[m]);
      }
    });
  }
};

// ======================================================================
// p_1st

template <typename MF>
struct Moment_p_1st
{
  using Mfields = MF;

  constexpr static char const* name = "p_1st";

  static int n_comps(const Grid_t& grid) { return 3 * grid.kinds.size(); }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"px", "py", "pz"}, grid.kinds);
  }

  template <typename Mparticles>
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    using Particle = typename Mparticles::ConstAccessor::Particle;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      int mm = prt.kind() * 3;
      auto pxi = prt.u();
      for (int m = 0; m < 3; m++) {
        deposit(prt, mm + m, prt.m() * pxi[m]);
      }
    });
  }
};

// ======================================================================
// T_1st

template <typename MF>
struct Moment_T_1st
{
  using Mfields = MF;

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
    using Particle = typename Mparticles::ConstAccessor::Particle;
    using Real = typename Particle::real_t;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      int mm = prt.kind() * 6;

      Real vxi[3];
      particle_calc_vxi(prt, vxi);
      auto pxi = prt.u();
      deposit(prt, mm + 0, prt.m() * pxi[0] * vxi[0]);
      deposit(prt, mm + 1, prt.m() * pxi[1] * vxi[1]);
      deposit(prt, mm + 2, prt.m() * pxi[2] * vxi[2]);
      deposit(prt, mm + 3, prt.m() * pxi[0] * vxi[1]);
      deposit(prt, mm + 4, prt.m() * pxi[0] * vxi[2]);
      deposit(prt, mm + 5, prt.m() * pxi[1] * vxi[2]);
    });
  }
};

// ======================================================================
// Moments_1st
//
// all moments calculated at once
// FIXME: add KE

template <typename MF>
class Moments_1st : public ItemMomentCRTP<Moments_1st<MF>, MF>
{
public:
  using Base = ItemMomentCRTP<Moments_1st<MF>, MF>;
  using Mfields = MF;

  constexpr static int n_moments = 13;
  static char const* name() { return "all_1st"; }

  static int n_comps(const Grid_t& grid)
  {
    return n_moments * grid.kinds.size();
  }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         grid.kinds);
  }

  Moments_1st(const Grid_t& grid) : Base{grid} {}

  template <typename Mparticles>
  void operator()(Mparticles& mprts)
  {
    using Particle = typename Mparticles::ConstAccessor::Particle;
    using Real = typename Particle::real_t;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, Base::mres_};
    deposit.process([&](const Particle& prt) {
      int mm = prt.kind() * n_moments;
      Real vxi[3];
      particle_calc_vxi(prt, vxi);
      deposit(prt, mm + 0, prt.q());
      deposit(prt, mm + 1, prt.q() * vxi[0]);
      deposit(prt, mm + 2, prt.q() * vxi[1]);
      deposit(prt, mm + 3, prt.q() * vxi[2]);
      deposit(prt, mm + 4, prt.m() * prt.u()[0]);
      deposit(prt, mm + 5, prt.m() * prt.u()[1]);
      deposit(prt, mm + 6, prt.m() * prt.u()[2]);
      deposit(prt, mm + 7, prt.m() * prt.u()[0] * vxi[0]);
      deposit(prt, mm + 8, prt.m() * prt.u()[1] * vxi[1]);
      deposit(prt, mm + 9, prt.m() * prt.u()[2] * vxi[2]);
      deposit(prt, mm + 10, prt.m() * prt.u()[0] * vxi[1]);
      deposit(prt, mm + 11, prt.m() * prt.u()[1] * vxi[2]);
      deposit(prt, mm + 12, prt.m() * prt.u()[2] * vxi[0]);
    });
    Base::bnd_.add_ghosts(Base::mres_);
  }
};
