
#pragma once

#include "bnd.hxx"

template<typename MfieldsState>
struct BndVpic : BndBase
{
  BndVpic(const Grid_t& grid, int ibn[3])
  {}

  void fill_ghosts(MfieldsState& mflds, int mb, int me) {}
  void add_ghosts(MfieldsState& mflds, int mb, int me) {}
};

