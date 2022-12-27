
#pragma once

#include "fields_item.hxx"
#include "bnd_cuda_3_impl.hxx"
#include "psc_fields_cuda.h"
#include "cuda_moments.cuh"

template <typename BS>
struct cuda_mparticles;

// ======================================================================
// Moment_rho_1st_nc_cuda

template <typename D>
class moment_rho_1st_nc_cuda
{
public:
  using dim_t = D;

  static std::string name() { return "rho_1st_nc"; }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return {"rho"};
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    auto& cmprts = *const_cast<MP&>(mprts).cmprts();

    CudaMoments1stNcRho<cuda_mparticles<typename MP::BS>, dim_t> cmoments;
    cmoments(cmprts, mflds_gt, ib);
  }
};

template <typename dim_t>
struct Moment_rho_1st_nc_cuda
  : ItemMomentCRTP<Moment_rho_1st_nc_cuda<dim_t>, MfieldsCuda::Storage,
                   BndCuda3>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc_cuda<dim_t>,
                              MfieldsCuda::Storage, BndCuda3>;
  using moment_type = moment_rho_1st_nc_cuda<dim_t>;

  using Base::Base;
};

// ======================================================================
// Moment_n_1st_cuda

template <typename D>
class moment_n_1st_cc_cuda
{
public:
  using dim_t = D;

  static std::string name() { return "n_1st_cc"; }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return addKindSuffix({"n"}, kinds);
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    auto& cmprts = *const_cast<MP&>(mprts).cmprts();

    CudaMoments1stN<cuda_mparticles<typename MP::BS>, dim_t> cmoments;
    cmoments(cmprts, mflds_gt, ib);
  }
};

template <typename dim_t>
class Moment_n_1st_cuda
  : public ItemMomentCRTP<Moment_n_1st_cuda<dim_t>, MfieldsCuda::Storage,
                          BndCuda3>
{
public:
  using Base =
    ItemMomentCRTP<Moment_n_1st_cuda<dim_t>, MfieldsCuda::Storage, BndCuda3>;
  using moment_type = moment_n_1st_cc_cuda<dim_t>;

  using Base::Base;
};

// ======================================================================
// Moments_1st_cuda

template <typename D>
class moments_1st_cc_cuda
{
public:
  using dim_t = D;

  static std::string name() { return "all_1st_cc"; }
  static std::vector<std::string> comp_names(const Grid_t::Kinds& kinds)
  {
    return addKindSuffix({"rho", "jx", "jy", "jz", "px", "py", "pz", "txx",
                          "tyy", "tzz", "txy", "tyz", "tzx"},
                         kinds);
  }

  template <typename MFLDS_GT, typename MP>
  void operator()(MFLDS_GT& mflds_gt, const Int3& ib, const MP& mprts)
  {
    auto& cmprts = *const_cast<MP&>(mprts).cmprts();

    CudaMoments1stAll<cuda_mparticles<typename MP::BS>, dim_t> cmoments;
    cmoments(cmprts, mflds_gt, ib);
  }
};

template <typename dim_t>
class Moments_1st_cuda
  : public ItemMomentCRTP<Moments_1st_cuda<dim_t>, MfieldsCuda::Storage,
                          BndCuda3>
{
public:
  using Base =
    ItemMomentCRTP<Moments_1st_cuda<dim_t>, MfieldsCuda::Storage, BndCuda3>;
  using moment_type = moments_1st_cc_cuda<dim_t>;

  using Base::Base;
};
