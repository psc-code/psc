
#pragma once

#include "particles.hxx"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// InjectBase

struct InjectBase
{
  InjectBase(int interval, int tau, int kind_n)
    : interval(interval), tau(tau), kind_n(kind_n)
  {}

  virtual void run(MparticlesBase& mprts_base, MfieldsBase& mflds_base) = 0;
  
  // param
  const int interval; // inject every so many steps
  const int tau; // in steps
  const int kind_n; // the number of particles to inject are based on this kind's density
};

