
#include "psc_sort_private.h"
#include "psc_particles_as_single.h"

#include <mrc_profile.h>
#include <mrc_params.h>

#include <string.h>
#include <psc_sort_common.c>

// ======================================================================
// psc_sort: subclass "qsort_single"

struct psc_sort_ops psc_sort_qsort_single_ops = {
  .name                  = "qsort_single",
  .run                   = psc_sort_qsort_run,
};

// ======================================================================
// psc_sort: subclass "countsort_single"

struct psc_sort_ops psc_sort_countsort_single_ops = {
  .name                  = "countsort_single",
  .run                   = psc_sort_countsort_run,
};

// ======================================================================
// psc_sort: subclass "countsort2_single"

struct psc_sort_ops psc_sort_countsort2_single_ops = {
  .name                  = "countsort2_single",
  .size                  = sizeof(struct psc_sort_countsort2),
  .param_descr           = psc_sort_countsort2_descr,
  .run                   = psc_sort_countsort2_run,
};

