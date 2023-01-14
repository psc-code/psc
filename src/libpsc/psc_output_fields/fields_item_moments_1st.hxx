
#pragma once

#include <psc/moment.hxx>
#include "fields_item.hxx"

// ======================================================================
// 1st cc

template <typename S, typename D>
using Moment_n_1st =
  ItemMoment<psc::moment::moment_n<psc::deposit::code::Deposit1stCc, D>, S>;

template <typename S, typename D>
using Moment_v_1st =
  ItemMoment<psc::moment::moment_v<psc::deposit::code::Deposit1stCc, D>, S>;

template <typename S, typename D>
using Moment_p_1st =
  ItemMoment<psc::moment::moment_p<psc::deposit::code::Deposit1stCc, D>, S>;

template <typename S, typename D>
using Moment_T_1st =
  ItemMoment<psc::moment::moment_T<psc::deposit::code::Deposit1stCc, D>, S>;

template <typename S, typename D>
using Moments_1st =
  ItemMoment<psc::moment::moment_all<psc::deposit::code::Deposit1stCc, D>, S>;

// ======================================================================
// 1st nc

template <typename S, typename D>
using Moment_rho_1st_nc =
  ItemMoment<psc::moment::moment_rho<psc::deposit::code::Deposit1stNc, D>, S>;

// ======================================================================
// 2nd nc

template <typename S, typename D>
using Moment_rho_2nd_nc =
  ItemMoment<psc::moment::moment_rho<psc::deposit::code::Deposit2ndNc, D>, S>;

#ifdef USE_CUDA

#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "psc_particles_single.h"

// ======================================================================
// Moments_1st_to_host

template <typename D>
class Moments_1st_to_host
{
public:
  using dim_t = D;
  using storage_type = MfieldsSingle::Storage;
  using value_type = storage_type::value_type;
  using space_type = storage_type::space_type;
  using Sub = Moments_1st<storage_type, D>;

  static std::string name() { return Sub::name(); }
  int n_comps() { return sub_.n_comps(); }
  const std::vector<std::string>& comp_names() { return sub_.comp_names(); }

  explicit Moments_1st_to_host(const Grid_t& grid) : sub_{grid} {}

  template <typename Mparticles>
  auto operator()(const Mparticles& _mprts)
  {
    static int pr, pr_A, pr_B;
    if (!pr) {
      pr = prof_register("Moments_1st cuda", 1., 0, 0);
      pr_A = prof_register("Moments_1st get", 1., 0, 0);
      pr_B = prof_register("Moments_1st process", 1., 0, 0);
    }

    prof_start(pr);
    prof_start(pr_A);
    auto& d_mprts = const_cast<Mparticles&>(_mprts);
    auto&& h_mprts = d_mprts.template get_as<MparticlesSingle>();
    prof_stop(pr_A);

    prof_start(pr_B);
    auto mres = sub_(h_mprts);
    prof_stop(pr_B);

    d_mprts.put_as(h_mprts, MP_DONT_COPY);
    prof_stop(pr);
    return mres;
  }

private:
  Sub sub_;
};

#endif
