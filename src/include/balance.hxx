
#pragma once

#include "psc.h"
#include "particles.hxx"

// ======================================================================
// BalanceBase

struct BalanceBase
{
  virtual void initial(Grid_t*& grid, std::vector<uint>& n_prts_by_patch) = 0;
  virtual void operator()(Grid_t*& grid, MparticlesBase& mp) = 0;
};

extern int psc_balance_generation_cnt;
