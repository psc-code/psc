
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
    using Particle = typename Mparticles::ConstAccessor::Particle;

    Base::mres_.storage().view() = 0.;
    auto deposit = Deposit1stCc<Mparticles, Mfields, dim_t>{mprts, Base::mres_};
    deposit.process([&](const Particle& prt) {
      int m = prt.kind();
      deposit(prt, m, 1.f);
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
struct Moment_v_1st
{
  using Mfields = MF;
  using dim_t = D;

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

    auto deposit = Deposit1stCc<Mparticles, Mfields, dim_t>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      Real vxi[3];
      _particle_calc_vxi(prt, vxi);

      int mm = prt.kind() * 3;
      for (int m = 0; m < 3; m++) {
        deposit(prt, mm + m, vxi[m]);
      }
    });
  }
};

// ======================================================================
// p_1st

template <typename MF, typename D>
struct Moment_p_1st
{
  using Mfields = MF;
  using dim_t = D;

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

    auto deposit = Deposit1stCc<Mparticles, Mfields, dim_t>{mprts, mflds};
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
    using Particle = typename Mparticles::ConstAccessor::Particle;
    using Real = typename Particle::real_t;

    auto deposit = Deposit1stCc<Mparticles, Mfields, dim_t>{mprts, mflds};
    deposit.process([&](const Particle& prt) {
      int mm = prt.kind() * 6;

      Real vxi[3];
      _particle_calc_vxi(prt, vxi);
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
    using Particle = typename Mparticles::ConstAccessor::Particle;
    using Real = typename Particle::real_t;

    auto deposit = Deposit1stCc<Mparticles, Mfields, dim_t>{mprts, Base::mres_};
    deposit.process([&](const Particle& prt) {
      int mm = prt.kind() * n_moments;
      Real vxi[3];
      _particle_calc_vxi(prt, vxi);
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

    auto deposit =
      Deposit1stCc<MparticlesSingle, Mfields, dim_t>{h_mprts, Base::mres_};

    prof_start(pr_B);
    printf("np %d\n", h_mprts.size());
    for (int p = 0; p < h_mprts.n_patches(); p++) {
      deposit.flds_ = deposit.mflds_[p];
      auto flds = make_Fields3d<dim_xyz>(deposit.mflds_[p]);
      auto prts = h_mprts[p];
      for (const auto& prt : prts) {
        int mm = prt.kind * n_moments;
        Real vxi[3];
        __particle_calc_vxi(prt, vxi);

        auto xi = prt.x; /* don't shift back in time */
        R u = xi[0] * deposit.dxi_[0] - .5f;
        R v = xi[1] * deposit.dxi_[1] - .5f;
        R w = xi[2] * deposit.dxi_[2] - .5f;

        int jx = fint(u);
        int jy = fint(v);
        int jz = fint(w);
        R h1 = u - jx;
        R h2 = v - jy;
        R h3 = w - jz;

        R g[2][3] = {{1.f - h1, 1.f - h2, 1.f - h3}, {h1, h2, h3}};

        int jxd = 1, jyd = 1, jzd = 1;
        if (deposit.is_invar_[0]) {
          jx = 0;
          g[0][0] = 1.;
          g[1][0] = 0.;
          jxd = 0;
        }
        if (deposit.is_invar_[1]) {
          jy = 0;
          g[0][1] = 1.;
          g[1][1] = 0.;
          jyd = 0;
        }
        if (deposit.is_invar_[2]) {
          jz = 0;
          g[0][2] = 1.;
          g[1][2] = 0.;
          jzd = 0;
        }

        assert(jx >= -1 && jx < deposit.ldims_[0]);
        assert(jy >= -1 && jy < deposit.ldims_[1]);
        assert(jz >= -1 && jz < deposit.ldims_[2]);

        R fnq = h_mprts.prt_w(prt) * deposit.fnqs_;
        R q = h_mprts.prt_q(prt);
        R m = h_mprts.prt_m(prt);

        int stride = &flds(mm + 1, 0, 0, 0) - &flds(mm, 0, 0, 0);
        int d[3];
        for (d[2] = 0; d[2] <= jzd; d[2]++) {
          for (d[1] = 0; d[1] <= jyd; d[1]++) {
            for (d[0] = 0; d[0] <= jxd; d[0]++) {
              R fac = fnq * g[d[0]][0] * g[d[1]][1] * g[d[2]][2];
              R* pp = &flds(mm, jx + d[0], jy + d[1], jz + d[2]);
              pp[0 * stride] += fac * q;
              pp[1 * stride] += fac * q * vxi[0];
              pp[2 * stride] += fac * q * vxi[1];
              pp[3 * stride] += fac * q * vxi[2];
              pp[4 * stride] += fac * m * prt.u[0];
              pp[5 * stride] += fac * m * prt.u[1];
              pp[6 * stride] += fac * m * prt.u[2];
              pp[7 * stride] += fac * m * prt.u[0] * vxi[0];
              pp[8 * stride] += fac * m * prt.u[1] * vxi[1];
              pp[9 * stride] += fac * m * prt.u[2] * vxi[2];
              pp[10 * stride] += fac * m * prt.u[0] * vxi[1];
              pp[11 * stride] += fac * m * prt.u[1] * vxi[2];
              pp[12 * stride] += fac * m * prt.u[2] * vxi[0];
            }
          }
        }
      }
    }
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
