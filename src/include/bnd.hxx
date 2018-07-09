
#pragma once

#include "fields3d.hxx"

// ======================================================================
// BndBase

struct BndBase
{
  virtual void fill_ghosts(MfieldsBase& mflds, int mb, int me) = 0;
  virtual void add_ghosts(MfieldsBase& mflds, int mb, int me) = 0;
};

