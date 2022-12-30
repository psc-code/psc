
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
};
