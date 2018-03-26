
#pragma once

#include "psc_output_particles_private.h"
#include "particles.hxx"

// ======================================================================
// class OutputParticlesBase

class OutputParticlesBase
{
public:
  virtual void run(PscMparticlesBase mprts_base) = 0;
};

// ======================================================================
// PscOutputParticles

template<typename S>
struct PscOutputParticles
{
  using sub_t = S;

  static_assert(std::is_convertible<sub_t*, OutputParticlesBase*>::value,
  		"sub classes used in PscOutputParticlesParticles must derive from OutputParticlesBase");
  
  explicit PscOutputParticles(psc_output_particles *outp)
    : outp_(outp)
  {}

  void run(PscMparticlesBase mprts_base)
  {
    struct psc_output_particles_ops *ops = psc_output_particles_ops(outp_);
    assert(ops->run);
    psc_stats_start(st_time_output);
    ops->run(outp_, mprts_base.mprts());
    psc_stats_stop(st_time_output);
  }

  psc_push_fields *outp() { return outp_; }
  
  sub_t* operator->() { return sub(); }

  sub_t* sub() { return mrc_to_subobj(outp_, sub_t); }

private:
  psc_output_particles* outp_;
};

using PscOutputParticlesBase = PscOutputParticles<OutputParticlesBase>;
