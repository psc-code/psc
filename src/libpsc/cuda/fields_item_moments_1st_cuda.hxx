
#pragma once

#include "fields_item.hxx"
#include "bnd_cuda_3_impl.hxx"
#include "psc_fields_cuda.h"
#include "cuda_moments.cuh"

template <typename BS>
struct cuda_mparticles;

// ======================================================================
// Moment_rho_1st_nc_cuda

template <typename dim_t>
struct Moment_rho_1st_nc_cuda
  : ItemMomentCRTP<Moment_rho_1st_nc_cuda<dim_t>, MfieldsCuda,
                   BndCuda3<MfieldsCuda>>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc_cuda<dim_t>, MfieldsCuda,
                              BndCuda3<MfieldsCuda>>;
  using Mfields = MfieldsCuda;
  using moment_type =
    psc::moment::moment_rho<psc::deposit::code::Deposit1stNc, dim_t>;

  Moment_rho_1st_nc_cuda(const Grid_t& grid) : Base{grid} {}

  template <typename Mparticles>
  void operator()(Mparticles& mprts)
  {
    auto& cmprts = *mprts.cmprts();

    Base::mres_gt_.view() = 0.;
    CudaMoments1stNcRho<cuda_mparticles<typename Mparticles::BS>, dim_t>
      cmoments;
    cmoments(cmprts, Base::mres_gt_, Base::mres_ib_);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
  }
};

// ======================================================================
// Moment_n_1st_cuda

template <typename dim_t>
class Moment_n_1st_cuda
  : public ItemMomentCRTP<Moment_n_1st_cuda<dim_t>, MfieldsCuda,
                          BndCuda3<MfieldsCuda>>
{
public:
  using Base = ItemMomentCRTP<Moment_n_1st_cuda<dim_t>, MfieldsCuda,
                              BndCuda3<MfieldsCuda>>;
  using Mfields = MfieldsCuda;
  using moment_type =
    psc::moment::moment_n<psc::deposit::code::Deposit1stCc, dim_t>;

  template <typename Mparticles>
  explicit Moment_n_1st_cuda(const Mparticles& mprts) : Base{mprts.grid()}
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

    CudaMoments1stN<cuda_mparticles<typename Mparticles::BS>, dim_t> cmoments;
    cmoments(cmprts, Base::mres_gt_, Base::mres_ib_);

    prof_start(pr_2);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
    prof_stop(pr_2);

    prof_stop(pr);
  }
};

// ======================================================================
// Moment_1st_cuda

template <typename dim_t>
class Moment_1st_cuda
  : public ItemMomentCRTP<Moment_1st_cuda<dim_t>, MfieldsCuda,
                          BndCuda3<MfieldsCuda>>
{
public:
  using Base =
    ItemMomentCRTP<Moment_1st_cuda<dim_t>, MfieldsCuda, BndCuda3<MfieldsCuda>>;
  using Mfields = MfieldsCuda;
  using value_type = typename Mfields::real_t;
  using space = gt::space::device;
  using moment_type =
    psc::moment::moment_all<psc::deposit::code::Deposit1stCc, dim_t>;

  template <typename Mparticles>
  explicit Moment_1st_cuda(const Mparticles& mprts) : Base{mprts.grid()}
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
    CudaMoments1stAll<cuda_mparticles<typename Mparticles::BS>, dim_t> cmoments;
    cmoments(cmprts, Base::mres_gt_, Base::mres_ib_);
    prof_stop(pr_2);

    prof_start(pr_3);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
    prof_stop(pr_3);

    prof_stop(pr);
  }
};
