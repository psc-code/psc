
#include "psc_sort_impl.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"

// ======================================================================
// psc_sort: subclass "countsort_single"

struct psc_sort_ops_countsort_single : psc_sort_ops {
  psc_sort_ops_countsort_single() {
    name                  = "countsort_single";
    size                  = sizeof(psc_sort_countsort<PscMparticlesSingle>);
    run                   = psc_sort_countsort<PscMparticlesSingle>::run;
  }
} psc_sort_countsort_single_ops;

// ======================================================================
// psc_sort: subclass "countsort2_single"

struct psc_sort_ops_countsort2_single : psc_sort_ops {
  psc_sort_ops_countsort2_single() {
    name                  = "countsort2_single";
    size                  = sizeof(psc_sort_countsort2<PscMparticlesSingle>);
    run                   = psc_sort_countsort2<PscMparticlesSingle>::run;
  }
} psc_sort_countsort2_single_ops;

// ======================================================================
// psc_sort: subclass "countsort_double"

struct psc_sort_ops_countsort_double : psc_sort_ops {
  psc_sort_ops_countsort_double() {
    name                  = "countsort_double";
    size                  = sizeof(psc_sort_countsort<PscMparticlesDouble>);
    run                   = psc_sort_countsort<PscMparticlesDouble>::run;
  }
} psc_sort_countsort_double_ops;

// ======================================================================
// psc_sort: subclass "countsort2_double"

struct psc_sort_ops_countsort2_double : psc_sort_ops {
  psc_sort_ops_countsort2_double() {
    name                  = "countsort2_double";
    size                  = sizeof(psc_sort_countsort2<PscMparticlesDouble>);
    run                   = psc_sort_countsort2<PscMparticlesDouble>::run;
  }
} psc_sort_countsort2_double_ops;



