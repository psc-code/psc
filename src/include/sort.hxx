
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
    if (ppsc->timestep % sort_->every != 0)
      return;
    
    static int st_time_sort, pr;
    if (!st_time_sort) {
      st_time_sort = psc_stats_register("time sort");
      pr = prof_register("sort", 1., 0, 0);
    }
    
    psc_stats_start(st_time_sort);
    prof_start(pr);
    
    (*this)->run(mprts.mprts());
    
    psc_stats_stop(st_time_sort);
    prof_stop(pr);
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
};

