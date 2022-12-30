
#pragma once

#include "../libpsc/psc_push_fields/marder_impl.hxx"
#include "psc_particles_single.h"
#include "mparticles_cuda.hxx"
//#include "fields_item_dive_cuda.hxx"
#include "fields_item_moments_1st_cuda.hxx"
#include "cuda_bits.h"

#include <gtensor/reductions.h>

#include <mrc_io.h>

template <typename BS, typename D>
struct MarderCuda
  : public MarderCommon<MparticlesCuda<BS>, MfieldsStateCuda, MfieldsCuda, D,
                        Moment_rho_1st_nc_cuda<D>, BndCuda3>
{
  using Base = MarderCommon<MparticlesCuda<BS>, MfieldsStateCuda, MfieldsCuda,
                            D, Moment_rho_1st_nc_cuda<D>, BndCuda3>;
  using Mparticles = typename Base::Mparticles;
  using MfieldsState = typename Base::MfieldsState;
  using Mfields = typename Base::Mfields;
  using dim_t = typename Base::dim_t;
  using Bnd = typename Base::Bnd;
  using real_t = typename Mfields::real_t;
  using Moment_t = typename Base::Item_rho_t;
  using Base::bnd_;
  using Base::diffusion_;
  using Base::dump_;
  using Base::io_;
  using Base::loop_;
  using Base::res_;
  using Base::rho_;

  using Base::Base;

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MfieldsStateCuda& mflds, MparticlesCuda<BS>& mprts)
  {
    const auto& grid = mprts.grid();
    static int pr;
    if (!pr) {
      pr = prof_register("marder", 1., 0, 0);
    }

    prof_start(pr);
    // need to fill ghost cells first (should be unnecessary with only variant
    // 1) FIXME
    bnd_.fill_ghosts(mflds, EX, EX + 3);

    Moment_t item_rho{grid};
    auto&& rho = psc::mflds::interior(grid, item_rho(mprts));

    for (int i = 0; i < loop_; i++) {
      Base::calc_aid_fields(mflds, rho);
      Base::print_max(res_);
      Base::correct(mflds);
      bnd_.fill_ghosts(mflds, EX, EX + 3);
    }
    prof_stop(pr);
  }
};
