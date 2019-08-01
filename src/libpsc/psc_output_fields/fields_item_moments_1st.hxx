
#pragma once

#include "Deposit1stCc.h"

#include "fields_item.hxx"

#include <cmath>

// ======================================================================
// n_1st

template <typename MP, typename MF>
struct Moment_n_1st
{
  using Mparticles = MP;
  using Mfields = MF;

  constexpr static char const* name = "n_1st";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return {"n"}; }
  constexpr static int flags = POFI_BY_KIND;

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

template <typename MP, typename MF>
struct Moment_v_1st
{
  using Mparticles = MP;
  using Mfields = MF;

  constexpr static char const* name = "v_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return {"vx", "vy", "vz"}; }
  constexpr static int flags = POFI_BY_KIND;

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

template <typename MP, typename MF>
struct Moment_p_1st
{
  using Mparticles = MP;
  using Mfields = MF;

  constexpr static char const* name = "p_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return {"px", "py", "pz"}; }
  constexpr static int flags = POFI_BY_KIND;

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

template <typename MP, typename MF>
struct Moment_vv_1st
{
  using Mparticles = MP;
  using Mfields = MF;

  constexpr static char const* name = "vv_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return {"vxvx", "vyvy", "vzvz"}; }
  constexpr static int flags = POFI_BY_KIND;

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

template <typename MP, typename MF>
struct Moment_T_1st
{
  using Mparticles = MP;
  using Mfields = MF;

  constexpr static char const* name = "T_1st";
  constexpr static int n_comps = 6;
  constexpr static fld_names_t fld_names()
  {
    return {"Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz"};
  }
  constexpr static int flags = POFI_BY_KIND;

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

template <typename MP, typename MF>
struct Moment_Tvv_1st
{
  using Mparticles = MP;
  using Mfields = MF;

  constexpr static char const* name = "Tvv_1st";
  constexpr static int n_comps = 6;
  constexpr static fld_names_t fld_names()
  {
    return {"vxvx", "vyvy", "vzvz", "vxvy", "vxvz", "vyvz"};
  }
  constexpr static int flags = POFI_BY_KIND;

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
