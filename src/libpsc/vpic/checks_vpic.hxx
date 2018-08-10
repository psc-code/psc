
#pragma once

struct ChecksVpic : ChecksParams, ChecksBase
{
  ChecksVpic(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params)
  {}
  
  void continuity_before_particle_push(psc* psc, MparticlesVpic& mprts) {}
  void continuity_after_particle_push(psc* psc, MparticlesVpic& mprts, MfieldsStateVpic& mflds) {}
  void gauss(psc* psc, MparticlesVpic& mprts, MfieldsStateVpic& mflds) {}

  void continuity_before_particle_push(psc* psc, MparticlesBase& mprts) override {}
  void continuity_after_particle_push(psc* psc, MparticlesBase& mprts, MfieldsStateBase& mflds) override {}
  void gauss(psc* psc, MparticlesBase& mprts, MfieldsStateBase& mflds) override {}
};
