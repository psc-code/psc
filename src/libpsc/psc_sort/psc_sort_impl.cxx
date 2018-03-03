
#include "psc_sort_impl.hxx"
#include "sort.hxx"

template<typename Sort>
static void psc_sort_setup(psc_sort* _sort)
{
  PscSort<Sort> sort(_sort);
  new(sort.sub()) Sort;
}

template<typename Sort>
static void psc_sort_destroy(psc_sort* _sort)
{
  PscSort<Sort> sort(_sort);
  sort->~Sort();
}

#include "psc_particles_single.h"
#include "psc_particles_double.h"

// ======================================================================
// psc_sort: subclass "countsort_single"

struct psc_sort_ops_countsort_single : psc_sort_ops {
  using Sort = psc_sort_countsort<PscMparticlesSingle>;
  psc_sort_ops_countsort_single() {
    name                  = "countsort_single";
    size                  = sizeof(Sort);
    run                   = Sort::run;
    setup                 = psc_sort_setup<Sort>;
    destroy               = psc_sort_destroy<Sort>;
  }
} psc_sort_countsort_single_ops;

// ======================================================================
// psc_sort: subclass "countsort2_single"

struct psc_sort_ops_countsort2_single : psc_sort_ops {
  using Sort = psc_sort_countsort2<PscMparticlesSingle>;
  psc_sort_ops_countsort2_single() {
    name                  = "countsort2_single";
    size                  = sizeof(Sort);
    run                   = Sort::run;
    setup                 = psc_sort_setup<Sort>;
    destroy               = psc_sort_destroy<Sort>;
  }
} psc_sort_countsort2_single_ops;

// ======================================================================
// psc_sort: subclass "countsort_double"

struct psc_sort_ops_countsort_double : psc_sort_ops {
  using Sort = psc_sort_countsort<PscMparticlesDouble>;
  psc_sort_ops_countsort_double() {
    name                  = "countsort_double";
    size                  = sizeof(Sort);
    run                   = Sort::run;
    setup                 = psc_sort_setup<Sort>;
    destroy               = psc_sort_destroy<Sort>;
  }
} psc_sort_countsort_double_ops;

// ======================================================================
// psc_sort: subclass "countsort2_double"

struct psc_sort_ops_countsort2_double : psc_sort_ops {
  using Sort = psc_sort_countsort2<PscMparticlesDouble>;
  psc_sort_ops_countsort2_double() {
    name                  = "countsort2_double";
    size                  = sizeof(Sort);
    run                   = Sort::run;
    setup                 = psc_sort_setup<Sort>;
    destroy               = psc_sort_destroy<Sort>;
  }
} psc_sort_countsort2_double_ops;



