
#include "psc_sort_private.h"
#include "psc_particles_as_single.h"

#include <mrc_profile.h>
#include <mrc_params.h>

#include <string.h>

#include "psc_sort_common.cxx"

// ======================================================================
// psc_sort: subclass "qsort_single"

struct psc_sort_ops_qsort_single : psc_sort_ops {
  psc_sort_ops_qsort_single() {
    name                  = "qsort_single";
    run                   = psc_sort_qsort_run;
  }
} psc_sort_qsort_single_ops;

// ======================================================================
// psc_sort: subclass "countsort_single"

struct psc_sort_ops_countsort_single : psc_sort_ops {
  psc_sort_ops_countsort_single() {
    name                  = "countsort_single";
    run                   = psc_sort_countsort_run;
  }
} psc_sort_countsort_single_ops;

// ======================================================================
// psc_sort: subclass "countsort2_single"

struct psc_sort_ops_countsort2_single : psc_sort_ops {
  psc_sort_ops_countsort2_single() {
    name                  = "countsort2_single";
    size                  = sizeof(psc_sort_countsort2<mparticles_t>);
    run                   = psc_sort_countsort2<mparticles_t>::run;
  }
} psc_sort_countsort2_single_ops;


