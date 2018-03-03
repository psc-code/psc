
#include "psc_sort_impl.hxx"
#include "sort.hxx"

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
    PscSort<Sort> sort(_sort);
    sort->run(_sort, mprts);
  }
};

#include "psc_particles_single.h"
#include "psc_particles_double.h"

// ======================================================================
// psc_sort: subclass "countsort_single"

struct psc_sort_ops_countsort_single : psc_sort_ops {
  using PscSort = PscSortWrapper<psc_sort_countsort<PscMparticlesSingle>>;
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
  using PscSort = PscSortWrapper<psc_sort_countsort2<PscMparticlesSingle>>;
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
  using PscSort = PscSortWrapper<psc_sort_countsort<PscMparticlesDouble>>;
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
  using PscSort = PscSortWrapper<psc_sort_countsort2<PscMparticlesDouble>>;
  psc_sort_ops_countsort2_double() {
    name                  = "countsort2_double";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
    run                   = PscSort::run;
  }
} psc_sort_countsort2_double_ops;



