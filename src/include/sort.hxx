
#pragma once

#include "particles.hxx"
#include "psc_sort_private.h"

#include <mrc_profile.h>

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

// ======================================================================
// PscSort_
//
// wraps get_as / put_as

template<class Sort>
class PscSort_ : SortBase
{
public:
  using mparticles_t = typename Sort::mparticles_t;
  
  void run(struct psc_mparticles *mprts_base) override
  {
    mparticles_t mprts = mprts_base->get_as<mparticles_t>();
    
    sort_(mprts);
    
    mprts.put_as(mprts_base);
  }

private:
  Sort sort_{};
};

template<typename Sort>
class PscSortWrapper
{
public:
  const static size_t size = sizeof(Sort);
  
  static void setup(psc_sort* _sort)
  {
    PscSort<Sort> sort(_sort);
    new(sort.sub()) Sort;
  }

  static void destroy(psc_sort* _sort)
  {
    PscSort<Sort> sort(_sort);
    sort->~Sort();
  }

  static void run(struct psc_sort* _sort, struct psc_mparticles* mprts)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("sort", 1., 0, 0);
    }

    PscSort<Sort> sort(_sort);
    prof_start(pr);
    sort->run(mprts);
    prof_stop(pr);
  }
};

