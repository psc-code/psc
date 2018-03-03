
#include "psc_sort_impl.hxx"
#include "sort.hxx"

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

#include "psc_particles_single.h"
#include "psc_particles_double.h"

// ======================================================================
// psc_sort: subclass "countsort_single"

struct psc_sort_ops_countsort_single : psc_sort_ops {
  using PscSort = PscSortWrapper<PscSort_<SortCountsort<PscMparticlesSingle>>>;
  psc_sort_ops_countsort_single() {
    name                  = "countsort_single";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
    run                   = PscSort::run;
  }
} psc_sort_countsort_single_ops;

// ======================================================================
// psc_sort: subclass "countsort2_single"

struct psc_sort_ops_countsort2_single : psc_sort_ops {
  using PscSort = PscSortWrapper<PscSort_<SortCountsort2<PscMparticlesSingle>>>;
  psc_sort_ops_countsort2_single() {
    name                  = "countsort2_single";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
    run                   = PscSort::run;
  }
} psc_sort_countsort2_single_ops;

// ======================================================================
// psc_sort: subclass "countsort_double"

struct psc_sort_ops_countsort_double : psc_sort_ops {
  using PscSort = PscSortWrapper<PscSort_<SortCountsort<PscMparticlesDouble>>>;
  psc_sort_ops_countsort_double() {
    name                  = "countsort_double";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
    run                   = PscSort::run;
  }
} psc_sort_countsort_double_ops;

// ======================================================================
// psc_sort: subclass "countsort2_double"

struct psc_sort_ops_countsort2_double : psc_sort_ops {
  using PscSort = PscSortWrapper<PscSort_<SortCountsort2<PscMparticlesDouble>>>;
  psc_sort_ops_countsort2_double() {
    name                  = "countsort2_double";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
    run                   = PscSort::run;
  }
} psc_sort_countsort2_double_ops;



