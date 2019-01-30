
#pragma once

#include <fields3d.hxx>

// ======================================================================
// MarderBase

struct MarderBase
{
  virtual void run(MfieldsStateBase& mflds_base, MparticlesBase& mprts_base) = 0;
};

