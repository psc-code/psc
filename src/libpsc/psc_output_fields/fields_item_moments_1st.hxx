
#pragma once

#include <psc/moment.hxx>
#include "fields_item.hxx"

#include <cmath>

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

// ======================================================================
// Moments_1st
//
// all moments calculated at once
// FIXME: add KE

template <typename MP, typename S, typename D>
class Moments_1st : public ItemMomentCRTP<Moments_1st<MP, S, D>, S>
{
public:
  using Base = ItemMomentCRTP<Moments_1st<MP, S, D>, S>;
  using moment_type =
    psc::moment::moment_all<psc::deposit::code::Deposit1stCc, D>;

  using Base::Base;
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
class Moments_1st<MparticlesCuda<BS>, MfieldsSingle::Storage, D>
{
public:
  using Sub = Moments_1st<MparticlesSingle, MfieldsSingle::Storage, D>;
  using dim_t = D;
  using Mparticles = MparticlesCuda<BS>;
  using storage_type = MfieldsSingle::Storage;
  using value_type = storage_type::value_type;
  using space_type = storage_type::space_type;
  using moment_type =
    psc::moment::moment_all<psc::deposit::code::Deposit1stCc, dim_t>;

  static std::string name() { return Sub::name(); }
  int n_comps() { return sub_.n_comps(); }
  const std::vector<std::string>& comp_names() { return sub_.comp_names(); }

  explicit Moments_1st(const Grid_t& grid) : sub_{grid} {}

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
