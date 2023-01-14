
#pragma once

template <typename Mparticles, typename MfieldsState>
struct ChecksVpic : ChecksParams
{
  ChecksVpic(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params)
  {}

  void continuity_before_particle_push(Mparticles& mprts) {}
  void continuity_after_particle_push(Mparticles& mprts, MfieldsState& mflds) {}
  void gauss(Mparticles& mprts, MfieldsState& mflds) {}
};
