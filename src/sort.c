
#include "psc.h"
#include "util/profile.h"

#include <string.h>

// ======================================================================
// quicksort

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

// ======================================================================
// counting sort

static void
countsort_sort()
{
  static int pr;
  if (!pr) {
    pr = prof_register("countsort_sort", 1., 0, 0);
  }
  prof_start(pr);
  
  int N = psc.fld_size;
  unsigned int *cnts = malloc(N * sizeof(*cnts));
  memset(cnts, 0, N * sizeof(*cnts));

  // count
  for (int i = 0; i < psc.n_part; i++) {
    unsigned int cni = psc.f_part[i].cni;
    assert(cni < N);
    cnts[cni]++;
  }

  // calc offsets
  int cur = 0;
  for (int i = 1; i < N; i++) {
    int n = cnts[i];
    cnts[i] = cur;
    cur += n;
  }
  assert(cur == psc.n_part);

  // move into new position
  struct f_particle *f_part2 = malloc(psc.n_part * sizeof(*f_part2));
  for (int i = 0; i < psc.n_part; i++) {
    unsigned int cni = psc.f_part[i].cni;
    memcpy(&f_part2[cnts[cni]], &psc.f_part[i], sizeof(*f_part2));
    cnts[cni]++;
  }

  // back to in-place
  memcpy(psc.f_part, f_part2, psc.n_part * sizeof(*f_part2));

  free(f_part2);
  free(cnts);

  prof_stop(pr);
}

struct psc_sort_ops psc_sort_ops_countsort = {
  .name = "countsort",
  .sort = countsort_sort,
};
