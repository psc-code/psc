
#pragma once

struct ChecksParams
{
  int continuity_every_step;   // check charge continuity eqn every so many steps
  double continuity_threshold; // acceptable error in continuity eqn
  bool continuity_verbose;     // always print continuity error, even if acceptable
  bool continuity_dump_always; // always dump d_rho, div_j, even if acceptable

  int gauss_every_step;   // check Gauss's Law every so many steps
  double gauss_threshold; // acceptable error in Gauss's Law
  bool gauss_verbose;     // always print Gauss's Law error, even if acceptable
  bool gauss_dump_always; // always dump E, div_rho, even if acceptable
};

// ======================================================================
// class ChecksBase

class ChecksBase
{
public:
  virtual void continuity_before_particle_push(MparticlesBase& mprts) = 0;
  virtual void continuity_after_particle_push(MparticlesBase& mprts, MfieldsStateBase& mflds) = 0;
  virtual void gauss(MparticlesBase& mprts, MfieldsStateBase& mflds) = 0;
};

