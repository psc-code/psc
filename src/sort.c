
#include "psc.h"
#include "util/profile.h"

#include <string.h>

// ======================================================================
// quicksort

#if PARTICLES_BASE == PARTICLES_FORTRAN

static void
find_cell_indices(particles_base_t *pp)
{
  PIC_find_cell_indices(pp);
}

static inline int
get_cell_index(const particle_base_t *p)
{
  return p->cni;
}

#else

static void
find_cell_indices(particles_base_t *pp)
{
  assert(0);
}

static inline int
get_cell_index(const particle_base_t *p)
{
  assert(0);
}

#endif

static int
compare(const void *_a, const void *_b)
{
  const particle_base_t *a = _a, *b = _b;

  if (get_cell_index(a) < get_cell_index(b)) {
    return -1;
  } else if (get_cell_index(a) == get_cell_index(b)) {
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
  find_cell_indices(&psc.pp);
  qsort(psc.pp.particles, psc.pp.n_part, sizeof(*psc.pp.particles), compare);
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
  
  find_cell_indices(&psc.pp);

  int N = psc.fld_size;
  unsigned int *cnts = malloc(N * sizeof(*cnts));
  memset(cnts, 0, N * sizeof(*cnts));

  // count
  for (int i = 0; i < psc.pp.n_part; i++) {
    unsigned int cni = get_cell_index(&psc.pp.particles[i]);
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
  assert(cur == psc.pp.n_part);

  // move into new position
  particle_base_t *particles2 = malloc(psc.pp.n_part * sizeof(*particles2));
  for (int i = 0; i < psc.pp.n_part; i++) {
    unsigned int cni = get_cell_index(&psc.pp.particles[i]);
    memcpy(&particles2[cnts[cni]], &psc.pp.particles[i], sizeof(*particles2));
    cnts[cni]++;
  }

  // back to in-place
  memcpy(psc.pp.particles, particles2, psc.pp.n_part * sizeof(*particles2));

  free(particles2);
  free(cnts);

  prof_stop(pr);
}

struct psc_sort_ops psc_sort_ops_countsort = {
  .name = "countsort",
  .sort = countsort_sort,
};

// ======================================================================
// counting sort 2
// use a separate array of cell indices 

static void
countsort2_sort()
{
  static int pr;
  if (!pr) {
    pr = prof_register("countsort2_sort", 1., 0, 0);
  }

  find_cell_indices(&psc.pp);

  int N = psc.fld_size;

  unsigned int *cnis = malloc(psc.pp.n_part * sizeof(*cnis));
  for (int i = 0; i < psc.pp.n_part; i++) {
    cnis[i] = get_cell_index(&psc.pp.particles[i]);
    assert(cnis[i] < N);
  }

  prof_start(pr);
  
  unsigned int *cnts = malloc(N * sizeof(*cnts));
  memset(cnts, 0, N * sizeof(*cnts));

  // count
  for (int i = 0; i < psc.pp.n_part; i++) {
    unsigned int cni = cnis[i];
    cnts[cni]++;
  }

  // calc offsets
  int cur = 0;
  for (int i = 1; i < N; i++) {
    int n = cnts[i];
    cnts[i] = cur;
    cur += n;
  }
  assert(cur == psc.pp.n_part);

  // move into new position
  particle_base_t *particles2 = malloc(psc.pp.n_part * sizeof(*particles2));
  for (int i = 0; i < psc.pp.n_part; i++) {
    unsigned int cni = cnis[i];
    int n = 1;
    while (i+n < psc.pp.n_part && cnis[i+n] == cni) {
      n++;
    }
    memcpy(&particles2[cnts[cni]], &psc.pp.particles[i], n * sizeof(*particles2));
    cnts[cni] += n;
    i += n - 1;
  }

  // back to in-place
  memcpy(psc.pp.particles, particles2, psc.pp.n_part * sizeof(*particles2));

  free(particles2);
  free(cnts);

  prof_stop(pr);
}

struct psc_sort_ops psc_sort_ops_countsort2 = {
  .name = "countsort2",
  .sort = countsort2_sort,
};

// ======================================================================
// none sort
//
// This doesn't actually sort, it's used to turn sorting off.

static void
none_sort()
{
}

struct psc_sort_ops psc_sort_ops_none = {
  .name = "none",
  .sort = none_sort,
};
