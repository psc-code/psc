
#pragma once

#include "fields_item.hxx"
#include "bnd_cuda_3_impl.hxx"
#include "psc_fields_cuda.h"
#include "cuda_moments.cuh"

template <typename BS>
struct cuda_mparticles;

// ======================================================================
// Moment_rho_1st_nc_cuda

template <typename _Mparticles, typename dim>
struct Moment_rho_1st_nc_cuda
  : ItemMomentCRTP<Moment_rho_1st_nc_cuda<_Mparticles, dim>, MfieldsCuda>
{
  using Base =
    ItemMomentCRTP<Moment_rho_1st_nc_cuda<_Mparticles, dim>, MfieldsCuda>;
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda3<Mfields>;

  static std::string name_impl() { return "rho_1st_nc"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return {"rho"};
  }

  Moment_rho_1st_nc_cuda(const Grid_t& grid) : Base{grid}, bnd_{grid, grid.ibn}
  {}

  void operator()(Mparticles& mprts)
  {
    auto& cmprts = *mprts.cmprts();

    Base::mres_gt_.view() = 0.;
    CudaMoments1stNcRho<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, Base::mres_gt_, Base::mres_ib_);
    bnd_.add_ghosts(Base::mres_, 0, Base::mres_gt_.shape(3));
  }

private:
  Bnd bnd_;
};

// ======================================================================
// Moment_n_1st_cuda

template <typename _Mparticles, typename dim>
class Moment_n_1st_cuda
  : public ItemMomentCRTP<Moment_n_1st_cuda<_Mparticles, dim>, MfieldsCuda>
{
public:
  using Base = ItemMomentCRTP<Moment_n_1st_cuda<_Mparticles, dim>, MfieldsCuda>;
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda3<Mfields>;

  constexpr static int n_moments = 1;
  static std::string name_impl() { return "n_1st_cuda"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  explicit Moment_n_1st_cuda(const Grid_t& grid)
    : Base{grid}, bnd_{grid, grid.ibn}
  {}

  explicit Moment_n_1st_cuda(const Mparticles& mprts)
    : Base{mprts.grid()}, bnd_{mprts.grid(), mprts.grid().ibn}
  {
    update(mprts);
  }

  void update(const Mparticles& mprts)
  {
    static int pr, pr_1, pr_2;
    if (!pr) {
      pr = prof_register("mom_n_cuda", 1, 0, 0);
      pr_1 = prof_register("mom_n_cuda zero", 1, 0, 0);
      pr_2 = prof_register("mom_n_cuda addg", 1, 0, 0);
    }

    prof_start(pr);
    auto& _mprts = const_cast<Mparticles&>(mprts);
    auto& cmprts = *_mprts.cmprts();

    prof_start(pr_1);
    Base::mres_gt_.view() = 0.;
    prof_stop(pr_1);

    CudaMoments1stN<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, Base::mres_gt_, Base::mres_ib_);

    prof_start(pr_2);
    bnd_.add_ghosts(Base::mres_, 0, Base::mres_.n_comps());
    prof_stop(pr_2);

    prof_stop(pr);
  }

private:
  Bnd bnd_;
};

// ======================================================================
// Moment_1st_cuda

template <typename _Mparticles, typename dim>
class Moment_1st_cuda
  : public ItemMomentCRTP<Moment_1st_cuda<_Mparticles, dim>, MfieldsCuda>
{
public:
  using Base = ItemMomentCRTP<Moment_1st_cuda<_Mparticles, dim>, MfieldsCuda>;
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda3<Mfields>;
  using value_type = typename Mfields::real_t;
  using space = gt::space::device;

  constexpr static int n_moments = 13;
  static std::string name_impl() { return "all_1st"; }
  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         grid.kinds);
  }

  explicit Moment_1st_cuda(const Grid_t& grid)
    : Base{grid}, bnd_{grid, grid.ibn}
  {}

  explicit Moment_1st_cuda(const Mparticles& mprts)
    : Base{mprts.grid()}, bnd_{mprts.grid(), mprts.grid().ibn}
  {
    update(mprts);
  }

  void update(const Mparticles& mprts)
  {
    static int pr, pr_1, pr_2, pr_3;
    if (!pr) {
      pr = prof_register("mom_cuda", 1, 0, 0);
      pr_1 = prof_register("mom_cuda zero", 1, 0, 0);
      pr_2 = prof_register("mom_cuda calc", 1, 0, 0);
      pr_3 = prof_register("mom_cuda addg", 1, 0, 0);
    }

    prof_start(pr);
    auto& _mprts = const_cast<Mparticles&>(mprts);
    auto& cmprts = *_mprts.cmprts();

    prof_start(pr_1);
    Base::mres_gt_.view() = 0.;
    prof_stop(pr_1);

    prof_start(pr_2);
    CudaMoments1stAll<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, Base::mres_gt_, Base::mres_ib_);
    prof_stop(pr_2);

    prof_start(pr_3);
    bnd_.add_ghosts(Base::mres_, 0, Base::mres_.n_comps());
    prof_stop(pr_3);

    prof_stop(pr);
  }

private:
  Bnd bnd_;
};
