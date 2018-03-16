
#pragma once

#include "particles.hxx"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// InjectBase

struct InjectBase
{
  InjectBase(bool do_inject, int every_step, int tau, int kind_n)
    : do_inject(do_inject),
      every_step(every_step), tau(tau), kind_n(kind_n)
  {}

  virtual void run(PscMparticlesBase mprts_base, PscMfieldsBase mflds_base) = 0;
  
  // param
  const bool do_inject; // whether to inject particles at all
  const int every_step; // inject every so many steps
  const int tau; // in steps
  const int kind_n; // the number of particles to inject are based on this kind's density
};

