
#include "psc.h"
#include <mrc_profile.h>

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
  foreach_patch(p) {
    particles_base_t *pp = &psc.particles.p[p];
    find_cell_indices(pp);
    qsort(pp->particles, pp->n_part, sizeof(*pp->particles), compare);
  }
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
 
  struct psc_patch *patch = &psc.patch[0];
  particles_base_t *pp = &psc.particles.p[0];
  find_cell_indices(pp);

  int N = 1;
  for (int d = 0; d < 3; d++) {
    N *= patch->ldims[d] + 2 * psc.ibn[d];
  }
  unsigned int *cnts = malloc(N * sizeof(*cnts));
  memset(cnts, 0, N * sizeof(*cnts));

  // count
  for (int i = 0; i < pp->n_part; i++) {
    unsigned int cni = get_cell_index(&pp->particles[i]);
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
  assert(cur == pp->n_part);

  // move into new position
  particle_base_t *particles2 = malloc(pp->n_part * sizeof(*particles2));
  for (int i = 0; i < pp->n_part; i++) {
    unsigned int cni = get_cell_index(&pp->particles[i]);
    memcpy(&particles2[cnts[cni]], &pp->particles[i], sizeof(*particles2));
    cnts[cni]++;
  }

  // back to in-place
  memcpy(pp->particles, particles2, pp->n_part * sizeof(*particles2));

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

  struct psc_patch *patch = &psc.patch[0];
  particles_base_t *pp = &psc.particles.p[0];
  find_cell_indices(pp);

  int N = 1;
  for (int d = 0; d < 3; d++) {
    N *= patch->ldims[d] + 2 * psc.ibn[d];
  }

  unsigned int *cnis = malloc(pp->n_part * sizeof(*cnis));
  for (int i = 0; i < pp->n_part; i++) {
    cnis[i] = get_cell_index(&pp->particles[i]);
    assert(cnis[i] < N);
  }

  prof_start(pr);
  
  unsigned int *cnts = malloc(N * sizeof(*cnts));
  memset(cnts, 0, N * sizeof(*cnts));

  // count
  for (int i = 0; i < pp->n_part; i++) {
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
  assert(cur == pp->n_part);

  // move into new position
  particle_base_t *particles2 = malloc(pp->n_part * sizeof(*particles2));
  for (int i = 0; i < pp->n_part; i++) {
    unsigned int cni = cnis[i];
    int n = 1;
    while (i+n < pp->n_part && cnis[i+n] == cni) {
      n++;
    }
    memcpy(&particles2[cnts[cni]], &pp->particles[i], n * sizeof(*particles2));
    cnts[cni] += n;
    i += n - 1;
  }

  // back to in-place
  memcpy(pp->particles, particles2, pp->n_part * sizeof(*particles2));

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
