
#pragma once

#include "bnd.hxx"

struct BndVpic : BndBase
{
  BndVpic(const Grid_t& grid, mrc_domain *domain, int ibn[3])
  {}

  void fill_ghosts(MfieldsStateVpic& mflds, int mb, int me) {}
  void add_ghosts(MfieldsStateVpic& mflds, int mb, int me) {}

  void fill_ghosts(MfieldsBase& mflds_base, int mb, int me) override { assert(0); }
  void add_ghosts(MfieldsBase& mflds_base, int mb, int me) override { assert(0); }
};

