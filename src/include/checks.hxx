
#pragma once

struct ChecksParams
{
  int continuity_every_step = 0;       // check charge continuity eqn every so many steps
  double continuity_threshold = 1e-13; // acceptable error in continuity eqn
  bool continuity_verbose = true;      // always print continuity error, even if acceptable
  bool continuity_dump_always = false; // always dump d_rho, div_j, even if acceptable

  int gauss_every_step = 0;       // check Gauss's Law every so many steps
  double gauss_threshold = 1e-13; // acceptable error in Gauss's Law
  bool gauss_verbose = true;      // always print Gauss's Law error, even if acceptable
  bool gauss_dump_always = false; // always dump E, div_rho, even if acceptable
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

