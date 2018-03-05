
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
    sub()->run(mprts);
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
  virtual void run(PscMparticlesBase mprts_base) = 0;
};

using PscSortBase = PscSort<SortBase>;

template<typename Derived, typename MP>
struct SortCRTP : SortBase
{
  using Mparticles = MP;
  
  SortCRTP(int interval)
    : interval_(interval)
  {}

  // FIXME, the 2x replicated funcs here aren't nice to start with.
  // There should be a way to have a get_as that's specialized for known types,
  // so that in case of the two known types being equal, nothing gets done..
  
  void operator()(Mparticles& mprts)
  {
    if (interval_ <= 0 || ppsc->timestep % interval_ != 0)
      return;
    
    static int st_time_sort, pr;
    if (!st_time_sort) {
      st_time_sort = psc_stats_register("time sort");
      pr = prof_register("sort", 1., 0, 0);
    }
    
    psc_stats_start(st_time_sort);

    MPI_Comm comm = MPI_COMM_WORLD; // FIXME
    mpi_printf(comm, "***** Sorting...\n");
    prof_start(pr);
    static_cast<Derived*>(this)->sort(mprts);
    prof_stop(pr);

    psc_stats_stop(st_time_sort);
  }
  
  void run(PscMparticlesBase mprts_base) override
  {
    if (interval_ <= 0 || ppsc->timestep % interval_ != 0)
      return;
    
    static int st_time_sort, pr;
    if (!st_time_sort) {
      st_time_sort = psc_stats_register("time sort");
      pr = prof_register("sort", 1., 0, 0);
    }
    
    psc_stats_start(st_time_sort);
    
    auto mprts = mprts_base.get_as<PscMparticles<Mparticles>>();
    MPI_Comm comm = MPI_COMM_WORLD; // FIXME
    mpi_printf(comm, "***** Sorting...\n");
    prof_start(pr);
    static_cast<Derived*>(this)->sort(*mprts.sub());
    prof_stop(pr);
    mprts.put_as(mprts_base);

    psc_stats_stop(st_time_sort);
  }

private:
  int interval_;
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
    new(sort.sub()) Sort(_sort->every);
  }

  static void destroy(psc_sort* _sort)
  {
    PscSort<Sort> sort(_sort);
    sort->~Sort();
  }
};

