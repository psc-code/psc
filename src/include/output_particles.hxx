
#pragma once

#include "psc_output_particles_private.h"
#include "particles.hxx"
#include "psc.h"

// ======================================================================
// class OutputParticlesBase

class OutputParticlesBase
{
public:
  virtual void run(MparticlesBase& mprts_base) = 0;

public:
  bool inited = true; // FIXME hack to avoid dtor call when not yet constructed
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

  void run(MparticlesBase& mprts_base)
  {
    psc_stats_start(st_time_output);
    sub()->run(mprts_base);
    psc_stats_stop(st_time_output);
  }

  psc_output_particles* outp() { return outp_; }
  
  sub_t* operator->() { return sub(); }

  sub_t* sub() { return mrc_to_subobj(outp_, sub_t); }

private:
  psc_output_particles* outp_;
};

using PscOutputParticlesBase = PscOutputParticles<OutputParticlesBase>;

// ======================================================================
// OutputParticlesWrapper

template<typename OutputParticles_t>
class OutputParticlesWrapper
{
public:
  const static size_t size = sizeof(OutputParticles_t);
  
  static void setup(struct psc_output_particles *_outp)
  {
    PscOutputParticles<OutputParticles_t> outp(_outp);
    new(outp.sub()) OutputParticles_t{*ggrid, _outp->params};
  }

  static void destroy(struct psc_output_particles *_outp)
  {
    if (!mrc_to_subobj(_outp, OutputParticlesBase)->inited) return; // FIXME
    PscOutputParticles<OutputParticles_t> outp(_outp);
    outp.sub()->~OutputParticles_t();
  }
};

