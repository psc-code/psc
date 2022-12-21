
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
{
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda3<Mfields>;

  constexpr static const char* name = "rho_1st_nc";
  static int n_comps(const Grid_t&) { return 1; }
  static std::vector<std::string> comp_names_impl()
  {
    return {"rho_nc_cuda"};
  } // FIXME
  constexpr static int flags = 0;

  Moment_rho_1st_nc_cuda(const Grid_t& grid)
    : mres_{grid, n_comps(grid), grid.ibn}, bnd_{grid, grid.ibn}
  {}

  void operator()(Mparticles& mprts)
  {
    auto& cmprts = *mprts.cmprts();

    mres_.gt().view() = 0.;
    CudaMoments1stNcRho<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, mres_);
    bnd_.add_ghosts(mres_, 0, mres_.n_comps());
  }

  Mfields& result() { return mres_; }

  auto gt() { return view_interior(mres_.gt(), mres_.ibn()); }

private:
  Mfields mres_;
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
  static char const* name() { return "n_1st_cuda"; }

  static int n_comps_impl(const Grid_t& grid)
  {
    return n_moments * grid.kinds.size();
  }

  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"n"}, grid.kinds);
  }

  Int3 ibn() const { return {}; }

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
    Base::mres_.gt().view() = 0.;
    prof_stop(pr_1);

    CudaMoments1stN<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, Base::mres_);

    prof_start(pr_2);
    bnd_.add_ghosts(Base::mres_, 0, Base::mres_.n_comps());
    prof_stop(pr_2);

    prof_stop(pr);
  }

  const Mfields& result() const { return Base::mres_; }

  auto gt()
  {
    Int3 bnd = Base::mres_.ibn();
    return std::move(Base::mres_.gt())
      .view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]), _s(bnd[2], -bnd[2]));
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
  static std::string name() { return "all_1st"; }
  static int n_comps_impl(const Grid_t& grid)
  {
    return n_moments * grid.kinds.size();
  }

  static std::vector<std::string> comp_names_impl(const Grid_t& grid)
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         grid.kinds);
  }

  int n_comps() const { return Base::mres_.n_comps(); }
  Int3 ibn() const { return Base::mres_.ibn(); }

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
    Base::mres_.gt().view() = 0.;
    prof_stop(pr_1);

    prof_start(pr_2);
    CudaMoments1stAll<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, Base::mres_);
    prof_stop(pr_2);

    prof_start(pr_3);
    bnd_.add_ghosts(Base::mres_, 0, Base::mres_.n_comps());
    prof_stop(pr_3);

    prof_stop(pr);
  }

  const Mfields& result() const { return Base::mres_; }

  auto gt()
  {
    auto bnd = Base::mres_.ibn();
    return Base::mres_.gt().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                                 _s(bnd[2], -bnd[2]));
  }

private:
  Bnd bnd_;
};
