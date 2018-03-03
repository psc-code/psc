
#pragma once

#include "psc_sort_private.h"

// ======================================================================
// PscSort

template<typename S>
struct PscSort
{
  using sub_t = S;
  
  explicit PscSort(psc_sort *sort)
    : sort_(sort)
  {}
  
  void operator()(PscMparticlesBase mprts)
  {
    psc_sort_run(sort_, mprts.mprts());
    // sub()->run(mprts.mprts());
  }
  
  sub_t* sub() { return mrc_to_subobj(sort_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_sort *sort_;
};

// ======================================================================
// SortBase

class SortBase
{
public:
  virtual void run(struct psc_mparticles *mprts_base) = 0;
};

using PscSortBase = PscSort<SortBase>;
