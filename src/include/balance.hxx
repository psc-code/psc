
#pragma once

#include "psc.h"
#include "particles.hxx"

// ======================================================================
// BalanceBase

struct BalanceBase
{
  virtual std::vector<uint> initial(psc* psc, const std::vector<uint>& n_prts_by_patch) = 0;
  virtual void operator()(psc* psc, MparticlesBase& mp) = 0;
};

extern int psc_balance_generation_cnt;

