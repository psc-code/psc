
#pragma once

#include <psc/moment.hxx>
#include "fields_item.hxx"

#include <cmath>

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

template <typename MF, typename D>
class Moment_n_1st : public ItemMomentCRTP<Moment_n_1st<MF, D>, MF>
{
public:
  using Base = ItemMomentCRTP<Moment_n_1st<MF, D>, MF>;
  using Mfields = MF;
  using real_t = typename Mfields::real_t;
  using dim_t = D;
  using moment_type =
    psc::moment::moment_n<psc::deposit::code::Deposit1stCc, dim_t>;

  static std::string name_impl() { return "n_1st"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  template <typename MP>
  explicit Moment_n_1st(const MP& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    moment_type{}(Base::mres_gt_, Base::mres_ib_, mprts);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
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
  using moment_type =
    psc::moment::moment_v<psc::deposit::code::Deposit1stCc, dim_t>;

  static std::string name_impl() { return "v_1st"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"vx", "vy", "vz"}, grid.kinds);
  }

  template <typename Mparticles>
  explicit Moment_v_1st(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    moment_type{}(Base::mres_gt_, Base::mres_ib_, mprts);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
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
  using moment_type =
    psc::moment::moment_p<psc::deposit::code::Deposit1stCc, dim_t>;

  static std::string name_impl() { return "p_1st"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"px", "py", "pz"}, grid.kinds);
  }

  template <typename Mparticles>
  explicit Moment_p_1st(const Mparticles& mprts) : Base{mprts.grid()}
  {
    Base::mres_gt_.view() = 0.f;
    moment_type{}(Base::mres_gt_, Base::mres_ib_, mprts);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
  }
};

// ======================================================================
// T_1st

template <typename MF, typename D>
struct Moment_T_1st
{
  using Mfields = MF;
  using dim_t = D;
  using moment_type =
    psc::moment::moment_T<psc::deposit::code::Deposit1stCc, dim_t>;

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
    moment_type{}(mflds.storage(), mflds.ib(), mprts);
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
  using moment_type =
    psc::moment::moment_all<psc::deposit::code::Deposit1stCc, dim_t>;

  static std::string name_impl() { return "all_1st"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         grid.kinds);
  }

  explicit Moments_1st(const Mparticles& mprts) : Base{mprts.grid()}
  {
    moment_type{}(Base::mres_gt_, Base::mres_ib_, mprts);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
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
  using real_t = Mfields::real_t;

  using Sub = Moments_1st<MparticlesSingle, Mfields, D>;

  static std::string name_impl() { return Sub::name_impl(); }

  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return Sub::comp_names_impl(grid);
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

    prof_start(pr_B);
    typename Sub::moment_type()(Base::mres_gt_, Base::mres_ib_, h_mprts);
    prof_stop(pr_B);

    prof_start(pr_C);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
    prof_stop(pr_C);

    mprts.put_as(h_mprts, MP_DONT_COPY);
    prof_stop(pr);
  }
};

#endif
