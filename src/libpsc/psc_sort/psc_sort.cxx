
#include "psc_sort_private.h"

#include <mrc_params.h>

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
  { "every"             , VAR(every)               , PARAM_INT(0)         },
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

