
#pragma once

#include "particles.hxx"
#include "psc_sort_private.h"

#include <mrc_profile.h>

// ======================================================================
// PscSort

struct SortBase;

template<typename S>
struct PscSort
{
  using sub_t = S;
  
  // static_assert(std::is_convertible<sub_t*, SortBase*>::value,
  // 		"sub classes used in PscSort must derive from SortBase");
  
  explicit PscSort(psc_sort *sort)
    : sort_(sort)
  {}
  
  void operator()(PscMparticlesBase mprts)
  {
    if (sort_->every <= 0 || ppsc->timestep % sort_->every != 0)
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

struct SortBase
{
  virtual void run(struct psc_mparticles *mprts_base) = 0;
};

using PscSortBase = PscSort<SortBase>;

template<typename Derived>
struct SortCRTP : SortBase
{
  void run(struct psc_mparticles *mprts_base) override
  {
    using Mparticles = typename Derived::Mparticles;
    auto mprts = mprts_base->get_as<PscMparticles<Mparticles>>();
    auto& derived = *static_cast<Derived*>(this);
    derived.sort(*mprts.sub());
    mprts.put_as(mprts_base);
  }
};

// ======================================================================
// PscSortWrapper

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

