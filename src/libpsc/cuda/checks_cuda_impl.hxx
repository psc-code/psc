
#pragma once

#include "checks.hxx"

template<typename Mparticles>
struct ChecksCuda : ChecksBase, ChecksParams
{
  ChecksCuda(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params)
  {}

  void continuity_before_particle_push(Mparticles& mprts)
  {}
 
  void continuity_after_particle_push(Mparticles& mprts, MfieldsStateCuda& mflds)
  {}

  void gauss(Mparticles& mprts, MfieldsStateCuda& mflds)
  {}
};
