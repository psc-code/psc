
#pragma once

#include "fields3d.hxx"

// ======================================================================
// BndFieldsBase

struct BndFieldsBase
{
  virtual void fill_ghosts_E(MfieldsBase& mflds_base) = 0;
  virtual void fill_ghosts_H(MfieldsBase& mflds_base) = 0;
  virtual void add_ghosts_J(MfieldsBase& mflds_base) = 0;
};

