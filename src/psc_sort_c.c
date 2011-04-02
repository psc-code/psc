
#include "psc_sort_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

#include <string.h>

// ======================================================================
// quicksort

#if PARTICLES_BASE == PARTICLES_FORTRAN

static void
find_cell_indices(int p, particles_base_t *pp)
{
  PIC_find_cell_indices(pp);
}

static inline int
get_cell_index(int p, const particle_base_t *p)
{
  return p->cni;
}

#else

static void
find_cell_indices(int p, particles_base_t *pp)
{
}

static inline int
get_cell_index(int p, const particle_base_t *part)
{
  struct psc_patch *patch = &psc.patch[p];
  particle_base_real_t dxi = 1.f / psc.dx[0];
  particle_base_real_t dyi = 1.f / psc.dx[1];
  particle_base_real_t dzi = 1.f / psc.dx[2];
  int *ldims = patch->ldims;
  
  particle_base_real_t u = (part->xi - patch->xb[0]) * dxi;
  particle_base_real_t v = (part->yi - patch->xb[1]) * dyi;
  particle_base_real_t w = (part->zi - patch->xb[2]) * dzi;
  int j0 = particle_base_real_nint(u);
  int j1 = particle_base_real_nint(v);
  int j2 = particle_base_real_nint(w);
    
  assert(j0 >= 0 && j0 < ldims[0]);
  assert(j1 >= 0 && j1 < ldims[1]);
  assert(j2 >= 0 && j2 < ldims[2]);
  return ((j2) * ldims[1] + j1) * ldims[0] + j0;
}

#endif

static int
compare(const void *_a, const void *_b)
{
  const particle_base_t *a = _a, *b = _b;

  if (get_cell_index(0, a) < get_cell_index(0, b)) {
    return -1;
  } else if (get_cell_index(0, a) == get_cell_index(0, b)) {
    return 0;
  } else {
    return 1;
  }
}

static void
psc_sort_qsort_run(struct psc_sort *sort, mparticles_base_t *particles)
{

  static int pr;
  if (!pr) {
    pr = prof_register("qsort_sort", 1., 0, 0);
  }
  prof_start(pr);
  assert(psc.nr_patches == 1);
  foreach_patch(p) {
    particles_base_t *pp = &particles->p[p];
    find_cell_indices(p, pp);
    qsort(pp->particles, pp->n_part, sizeof(*pp->particles), compare);
  }
  prof_stop(pr);
}

// ======================================================================
// psc_sort: subclass "qsort"

struct psc_sort_ops psc_sort_qsort_ops = {
  .name                  = "qsort",
  .run                   = psc_sort_qsort_run,
};

// ======================================================================
// counting sort

static void
psc_sort_countsort_run(struct psc_sort *sort, mparticles_base_t *particles)
{
  static int pr;
  if (!pr) {
    pr = prof_register("countsort_sort", 1., 0, 0);
  }
  prof_start(pr);
 
  struct psc_patch *patch = &psc.patch[0];
  particles_base_t *pp = &particles->p[0];
  find_cell_indices(0, pp);

  int N = 1;
  for (int d = 0; d < 3; d++) {
    N *= patch->ldims[d] + 2 * psc.ibn[d];
  }
  unsigned int *cnts = malloc(N * sizeof(*cnts));
  memset(cnts, 0, N * sizeof(*cnts));

  // count
  for (int i = 0; i < pp->n_part; i++) {
    unsigned int cni = get_cell_index(0, &pp->particles[i]);
    assert(cni < N);
    cnts[cni]++;
  }

  // calc offsets
  int cur = 0;
  for (int i = 0; i < N; i++) {
    int n = cnts[i];
    cnts[i] = cur;
    cur += n;
  }
  assert(cur == pp->n_part);

  // move into new position
  particle_base_t *particles2 = malloc(pp->n_part * sizeof(*particles2));
  for (int i = 0; i < pp->n_part; i++) {
    unsigned int cni = get_cell_index(0, &pp->particles[i]);
    memcpy(&particles2[cnts[cni]], &pp->particles[i], sizeof(*particles2));
    cnts[cni]++;
  }

  // back to in-place
  memcpy(pp->particles, particles2, pp->n_part * sizeof(*particles2));

  free(particles2);
  free(cnts);

  prof_stop(pr);
}

// ======================================================================
// psc_sort: subclass "countsort"

struct psc_sort_ops psc_sort_countsort_ops = {
  .name                  = "countsort",
  .run                   = psc_sort_countsort_run,
};

// ======================================================================
// counting sort 2
// use a separate array of cell indices 

static void
psc_sort_countsort2_run(struct psc_sort *sort, mparticles_base_t *particles)
{
  static int pr;
  if (!pr) {
    pr = prof_register("countsort2_sort", 1., 0, 0);
  }

  struct psc_patch *patch = &psc.patch[0];
  particles_base_t *pp = &particles->p[0];
  find_cell_indices(0, pp);

  int N = 1;
  for (int d = 0; d < 3; d++) {
    N *= patch->ldims[d] + 2 * psc.ibn[d];
  }

  unsigned int *cnis = malloc(pp->n_part * sizeof(*cnis));
  for (int i = 0; i < pp->n_part; i++) {
    cnis[i] = get_cell_index(0, &pp->particles[i]);
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
  for (int i = 0; i < N; i++) {
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

// ======================================================================
// psc_sort: subclass "countsort2"

struct psc_sort_ops psc_sort_countsort2_ops = {
  .name                  = "countsort2",
  .run                   = psc_sort_countsort2_run,
};

