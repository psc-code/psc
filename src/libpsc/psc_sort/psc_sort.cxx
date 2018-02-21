
#include "psc_sort_private.h"

#include <mrc_params.h>

// ======================================================================
// forward to subclass

void
psc_sort_run(struct psc_sort *sort, struct psc_mparticles *mprts)
{
  if (ppsc->timestep % sort->every != 0)
    return;

  static int st_time_sort;
  if (!st_time_sort) {
    st_time_sort = psc_stats_register("time sort");
  }

  psc_stats_start(st_time_sort);

  struct psc_sort_ops *ops = psc_sort_ops(sort);
  assert(ops->run);
  ops->run(sort, mprts);

  psc_stats_stop(st_time_sort);
}

// ======================================================================
// psc_sort_init

extern struct psc_sort_ops psc_sort_none_ops;
extern struct psc_sort_ops psc_sort_countsort_single_ops;
extern struct psc_sort_ops psc_sort_countsort2_single_ops;
extern struct psc_sort_ops psc_sort_countsort_double_ops;
extern struct psc_sort_ops psc_sort_countsort2_double_ops;
extern struct psc_sort_ops psc_sort_vpic_ops;

static void
psc_sort_init()
{
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_none_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_countsort_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_countsort2_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_countsort_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_countsort2_double_ops);
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_vpic_ops);
#endif
}

// ======================================================================
// psc_sort class

#define VAR(x) (void *)offsetof(struct psc_sort, x)
static struct param psc_sort_descr[] = {
  { "every"             , VAR(every)               , PARAM_INT(100)       },
  {},
};
#undef VAR

struct mrc_class_psc_sort_ : mrc_class_psc_sort {
  mrc_class_psc_sort_() {
    name             = "psc_sort";
    size             = sizeof(struct psc_sort);
    param_descr      = psc_sort_descr;
    init             = psc_sort_init;
  }
} mrc_class_psc_sort;

