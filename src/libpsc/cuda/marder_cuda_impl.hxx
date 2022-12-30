
#pragma once

#include "../libpsc/psc_push_fields/marder_impl.hxx"
#include "psc_particles_single.h"
#include "mparticles_cuda.hxx"
//#include "fields_item_dive_cuda.hxx"
#include "fields_item_moments_1st_cuda.hxx"

template <typename BS, typename D>
using MarderCuda =
  MarderCommon<MparticlesCuda<BS>, MfieldsStateCuda, MfieldsCuda, D,
               Moment_rho_1st_nc_cuda<D>, BndCuda3>;
