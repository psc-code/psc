
#include "psc.h"
#include "util/profile.h"

static int
compare(const void *_a, const void *_b)
{
  const struct f_particle *a = _a, *b = _b;

  if (a->cni < b->cni) {
    return -1;
  } else if (a->cni == b->cni) {
    return 0;
  } else {
    return 1;
  }
}

static void
qsort_sort()
{
  static int pr;
  if (!pr) {
    pr = prof_register("qsort_sort", 1., 0, 0);
  }
  prof_start(pr);
  qsort(psc.f_part, psc.n_part, sizeof(*psc.f_part), compare);
  prof_stop(pr);
}

struct psc_sort_ops psc_sort_ops_qsort = {
  .name = "qsort",
  .sort = qsort_sort,
};
