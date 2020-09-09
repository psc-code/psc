
#pragma once

#include "particles.hxx"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// InjectBase

struct InjectBase
{
  InjectBase(int interval, int tau)
    : interval(interval), tau(tau)
  {}

  // param
  const int interval; // inject every so many steps
  const int tau; // in steps
};

