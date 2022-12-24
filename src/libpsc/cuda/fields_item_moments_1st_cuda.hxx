
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
  using storage_type = typename Base::storage_type;
  using real_t = typename Base::real_t;
  using space_type = typename Base::space_type;
  using moment_type =
    psc::moment::moment_rho<psc::deposit::code::Deposit1stNc, dim_t>;

  using Base::Base;

  template <typename Mparticles>
  auto operator()(Mparticles& mprts)
  {
    auto& cmprts = *mprts.cmprts();

    Int3 ib = -mprts.grid().ibn;
    // FIXME, gt::gtensor and psc::gtensor are slightly different, and ideally
    // zeros() shouldn't actually allocate, but probably it does, so this wastes
    // memory and a copy
    storage_type mres =
      psc::mflds::zeros<real_t, space_type>(mprts.grid(), 1, ib);
    CudaMoments1stNcRho<cuda_mparticles<typename Mparticles::BS>, dim_t>
      cmoments;
    cmoments(cmprts, mres, ib);
    Base::bnd_.add_ghosts(mprts.grid(), mres, ib);
    return mres;
  }

  auto storage() = delete;
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
  using storage_type = typename Base::storage_type;
  using real_t = typename Base::real_t;
  using space_type = typename Base::space_type;
  using moment_type =
    psc::moment::moment_n<psc::deposit::code::Deposit1stCc, dim_t>;

  using Base::Base;

  template <typename Mparticles>
  auto operator()(const Mparticles& mprts)
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
    Int3 ib = -mprts.grid().ibn;
    storage_type mres =
      psc::mflds::zeros<real_t, space_type>(mprts.grid(), 1, ib);
    prof_stop(pr_1);

    CudaMoments1stN<cuda_mparticles<typename Mparticles::BS>, dim_t> cmoments;
    cmoments(cmprts, mres, ib);

    prof_start(pr_2);
    Base::bnd_.add_ghosts(mprts.grid(), mres, ib);
    prof_stop(pr_2);

    prof_stop(pr);
    return mres;
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
  using storage_type = typename Base::storage_type;
  using real_t = typename Base::real_t;
  using space_type = typename Base::space_type;
  using moment_type =
    psc::moment::moment_all<psc::deposit::code::Deposit1stCc, dim_t>;

  using Base::Base;

  template <typename Mparticles>
  auto operator()(const Mparticles& mprts)
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
    Int3 ib = -mprts.grid().ibn;
    storage_type mres =
      psc::mflds::zeros<real_t, space_type>(mprts.grid(), 1, ib);
    prof_stop(pr_1);

    prof_start(pr_2);
    CudaMoments1stAll<cuda_mparticles<typename Mparticles::BS>, dim_t> cmoments;
    cmoments(cmprts, mres, ib);
    prof_stop(pr_2);

    prof_start(pr_3);
    Base::bnd_.add_ghosts(mprts.grid(), mres, ib);
    prof_stop(pr_3);

    prof_stop(pr);
    return mres;
  }
};
