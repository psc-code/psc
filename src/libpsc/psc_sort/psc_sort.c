
#include "psc_sort_private.h"

#include <mrc_params.h>

// ======================================================================
// forward to subclass

void
psc_sort_run(struct psc_sort *sort, mparticles_base_t *particles)
{
  if (ppsc->timestep % sort->every != 0)
    return;

  struct psc_sort_ops *ops = psc_sort_ops(sort);
  assert(ops->run);
  ops->run(sort, particles);
}

// ======================================================================
// psc_sort_init

static void
psc_sort_init()
{
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_none_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_fortran_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_qsort_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_countsort_ops);
  mrc_class_register_subclass(&mrc_class_psc_sort, &psc_sort_countsort2_ops);
}

// ======================================================================
// psc_sort class

#define VAR(x) (void *)offsetof(struct psc_sort, x)
static struct param psc_sort_descr[] = {
  { "every"             , VAR(every)               , PARAM_INT(100)       },
  {},
};
#undef VAR

struct mrc_class_psc_sort mrc_class_psc_sort = {
  .name             = "psc_sort",
  .size             = sizeof(struct psc_sort),
  .param_descr      = psc_sort_descr,
  .init             = psc_sort_init,
};

