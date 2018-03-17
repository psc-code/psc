
#pragma once

#include "particles.hxx"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// HeatingBase

struct HeatingBase
{
  virtual void run(PscMparticlesBase mprts_base) = 0;
};

