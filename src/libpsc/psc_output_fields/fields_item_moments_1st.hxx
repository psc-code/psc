
#pragma once

#include "Deposit1stCc.h"

#include "fields_item.hxx"

#include <cmath>

// ======================================================================
// n_1st

template <typename MF>
struct Moment_n_1st
{
  using Mfields = MF;

  constexpr static char const* name = "n_1st";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"n"}; }
  constexpr static int flags = POFI_BY_KIND;

  template <typename Mparticles>
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    using Particle = typename Mparticles::ConstAccessor::Particle;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      int m = prt.kind();
      deposit(prt, m, 1.f);
    });
  }
};

// ======================================================================
// v_1st

template <typename MF>
struct Moment_v_1st
{
  using Mfields = MF;

  constexpr static char const* name = "v_1st";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"vx", "vy", "vz"}; }
  constexpr static int flags = POFI_BY_KIND;

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
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"px", "py", "pz"}; }
  constexpr static int flags = POFI_BY_KIND;

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
// vv_1st

template <typename MF>
struct Moment_vv_1st
{
  using Mfields = MF;

  constexpr static char const* name = "vv_1st";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"vxvx", "vyvy", "vzvz"};
  }
  constexpr static int flags = POFI_BY_KIND;

  template <typename Mparticles>
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    using Particle = typename Mparticles::ConstAccessor::Particle;
    using Real = typename Particle::real_t;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      int mm = prt.kind() * 3;
      Real vxi[3];
      particle_calc_vxi(prt, vxi);

      for (int m = 0; m < 3; m++) {
        deposit(prt, mm + m, vxi[m] * vxi[m]);
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
  constexpr static int n_comps = 6;
  static std::vector<std::string> fld_names()
  {
    return {"Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz"};
  }
  constexpr static int flags = POFI_BY_KIND;

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
// Tvv_1st

template <typename MF>
struct Moment_Tvv_1st
{
  using Mfields = MF;

  constexpr static char const* name = "Tvv_1st";
  constexpr static int n_comps = 6;
  static std::vector<std::string> fld_names()
  {
    return {"vxvx", "vyvy", "vzvz", "vxvy", "vxvz", "vyvz"};
  }
  constexpr static int flags = POFI_BY_KIND;

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
      deposit(prt, mm + 0, prt.m() * vxi[0] * vxi[0]);
      deposit(prt, mm + 1, prt.m() * vxi[1] * vxi[1]);
      deposit(prt, mm + 2, prt.m() * vxi[2] * vxi[2]);
      deposit(prt, mm + 3, prt.m() * vxi[0] * vxi[1]);
      deposit(prt, mm + 4, prt.m() * vxi[0] * vxi[2]);
      deposit(prt, mm + 5, prt.m() * vxi[1] * vxi[2]);
    });
  }
};

// ======================================================================
// Moments_1st
//
// all moments calculated at once
// FIXME: add KE

template <typename MF>
struct Moments_1st
{
  using Mfields = MF;

  constexpr static char const* name = "all_1st";
  constexpr static int n_comps = 13;
  static std::vector<std::string> fld_names()
  {
    return {"rho", "jx",  "jy",  "jz",  "px",  "py", "pz",
            "txx", "tyy", "tzz", "txy", "tyz", "tzx"};
  }
  constexpr static int flags = POFI_BY_KIND;

  template <typename Mparticles>
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    using Particle = typename Mparticles::ConstAccessor::Particle;
    using Real = typename Particle::real_t;

    auto deposit = Deposit1stCc<Mparticles, Mfields>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      int mm = prt.kind() * n_comps;
      Real q = prt.q(), m = prt.m();
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
  }
};
