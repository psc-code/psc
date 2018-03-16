
#pragma once

#include "psc_sort_private.h"
#include "particles.hxx"

#include "psc_stats.h"
#include <mrc_profile.h>

// ======================================================================
// SortBase

struct SortBase
{
  virtual void run(PscMparticlesBase mprts_base) = 0;
};

// ======================================================================
// PscSort

template<typename S>
struct PscSort
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, SortBase*>::value,
  		"sub classes used in PscSort must derive from SortBase");
  
  explicit PscSort(psc_sort *sort)
    : sort_(sort)
  {}
  
  void operator()(PscMparticlesBase mprts)
  {
    static int st_time_sort;
    if (!st_time_sort) {
      st_time_sort = psc_stats_register("time sort");
    }
    
    psc_stats_start(st_time_sort);
    sub()->run(mprts);
    psc_stats_stop(st_time_sort);
  }
  
  sub_t* sub() { return mrc_to_subobj(sort_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_sort *sort_;
};

using PscSortBase = PscSort<SortBase>;

// ======================================================================
// SortConvert

template<typename Sort_t>
struct SortConvert : Sort_t, SortBase
{
  using Base = Sort_t;
  using Mparticles = typename Sort_t::Mparticles;
  using Base::Base;

  void run(PscMparticlesBase mprts_base) override
  {
    auto mprts = mprts_base.get_as<PscMparticles<Mparticles>>();
    (*this)(*mprts.sub());
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
    new(sort.sub()) Sort{};
  }

  static void destroy(psc_sort* _sort)
  {
    PscSort<Sort> sort(_sort);
    sort->~Sort();
  }
};

