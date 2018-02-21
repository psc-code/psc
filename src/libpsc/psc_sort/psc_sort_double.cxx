
#include "psc_sort_private.h"
#include "psc_particles_as_double.h"

#include <mrc_profile.h>
#include <mrc_params.h>

#include <string.h>

#include "psc_sort_common.cxx"

// ======================================================================
// psc_sort: subclass "countsort_double"

struct psc_sort_ops_countsort_double : psc_sort_ops {
  psc_sort_ops_countsort_double() {
    name                  = "countsort_double";
    run                   = psc_sort_countsort_run;
  }
} psc_sort_countsort_double_ops;

// ======================================================================
// psc_sort: subclass "countsort2_double"

struct psc_sort_ops_countsort2_double : psc_sort_ops {
  psc_sort_ops_countsort2_double() {
    name                  = "countsort2_double";
    size                  = sizeof(psc_sort_countsort2<mparticles_t>);
    run                   = psc_sort_countsort2<mparticles_t>::run;
  }
} psc_sort_countsort2_double_ops;


