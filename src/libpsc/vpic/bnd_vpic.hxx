
#pragma once

#include "bnd.hxx"

template <typename MfieldsState>
struct BndVpic : BndBase
{
  BndVpic() {}

  void fill_ghosts(MfieldsState& mflds, int mb, int me) {}
  void add_ghosts(MfieldsState& mflds, int mb, int me) {}
};
