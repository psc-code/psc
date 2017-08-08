
#include "psc_sort_private.h"
#include "psc_particles_as_c.h"

#include <mrc_profile.h>
#include <mrc_params.h>

#include <string.h>
#include <psc_sort_common.c>

// ======================================================================
// psc_sort: subclass "qsort"

struct psc_sort_ops psc_sort_qsort_ops = {
  .name                  = "qsort",
  .run                   = psc_sort_qsort_run,
};

// ======================================================================
// psc_sort: subclass "countsort"

struct psc_sort_ops psc_sort_countsort_ops = {
  .name                  = "countsort",
  .run                   = psc_sort_countsort_run,
};

// ======================================================================
// psc_sort: subclass "countsort2"

struct psc_sort_ops psc_sort_countsort2_ops = {
  .name                  = "countsort2",
  .size                  = sizeof(struct psc_sort_countsort2),
  .param_descr           = psc_sort_countsort2_descr,
  .run                   = psc_sort_countsort2_run,
};

