
#pragma once

#include "bnd_fields.hxx"

struct BndFieldsVpic : BndFieldsBase
{
  void fill_ghosts_E(MfieldsState& mflds) {}
  void fill_ghosts_H(MfieldsState& mflds) {}
  void add_ghosts_J(MfieldsState& mflds) {}

  void fill_ghosts_E(MfieldsBase& mflds_base) override {}
  void fill_ghosts_H(MfieldsBase& mflds_base) override {}
  void add_ghosts_J(MfieldsBase& mflds_base) override {}
};

