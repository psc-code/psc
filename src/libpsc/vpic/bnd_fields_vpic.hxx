
#pragma once

#include "bnd_fields.hxx"

template<typename MfieldsState>
struct BndFieldsVpic : BndFieldsBase
{
  void fill_ghosts_E(MfieldsState& mflds) {}
  void fill_ghosts_H(MfieldsState& mflds) {}
  void add_ghosts_J(MfieldsState& mflds) {}
};

