
#include "psc_sort_private.h"

#include <mrc_profile.h>

// ======================================================================
// cuda sort

static unsigned int *
create_indices(struct cell_map *map, particles_base_t *pp, struct psc_patch *patch)
{
  unsigned int *cnis = malloc(pp->n_part * sizeof(*cnis));
  for (int i = 0; i < pp->n_part; i++) {
    particle_base_t *p = &pp->particles[i];
    particle_base_real_t dxi = 1.f / ppsc->dx[0];
    particle_base_real_t dyi = 1.f / ppsc->dx[1];
    particle_base_real_t dzi = 1.f / ppsc->dx[2];
    particle_base_real_t xi[3] = {
      (p->xi - patch->xb[0]) * dxi,
      (p->yi - patch->xb[1]) * dyi,
      (p->zi - patch->xb[2]) * dzi };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = particle_base_real_nint(xi[d]);
    }
    
    cnis[i] = cell_map_3to1(map, pos);
    assert(cnis[i] < map->N);
  }
  return cnis;
}

static void
sort_pairs(unsigned int *keys, unsigned int *vals, int n, int n_max)
{
  unsigned int *cnts = calloc(n_max, sizeof(*cnts));
  
  // count
  for (int i = 0; i < n; i++) {
    unsigned int key = keys[i];
    cnts[key]++;
  }
  
  // calc offsets
  int cur = 0;
  for (int i = 0; i < n_max; i++) {
    int cnt = cnts[i];
    cnts[i] = cur;
    cur += cnt;
  }
  assert(cur == n);

  // move
  unsigned int *vals2 = malloc(n * sizeof(*vals2));
  for (int i = 0; i < n; i++) {
    unsigned int key = keys[i];
    vals2[cnts[key]++] = vals[i];
  }
  free(cnts);
  memcpy(vals, vals2, n * sizeof(*vals));
  free(vals2);
}

static void
sort_patch(int p, particles_base_t *pp)
{
  struct psc_patch *patch = &ppsc->patch[p];
  struct cell_map map;

  int N = cell_map_init(&map, patch->ldims, (int[3]) { 1, 8, 8 });
  unsigned int *cnis = create_indices(&map, pp, patch);
  cell_map_free(&map);
  unsigned int *ids = malloc(pp->n_part * sizeof(*ids));
  for (int i = 0; i < pp->n_part; i++) {
    ids[i] = i;
  }
  
  sort_pairs(cnis, ids, pp->n_part, N);
  free(cnis);

  // move into new position
  particle_base_t *particles2 = malloc(pp->n_part * sizeof(*particles2));
  for (int i = 0; i < pp->n_part; i++) {
    particles2[i] = pp->particles[ids[i]];
  }
  free(ids);

  // back to in-place
  memcpy(pp->particles, particles2, pp->n_part * sizeof(*particles2));
  
  free(particles2);
}

static void
psc_sort_cuda_run(struct psc_sort *sort, mparticles_base_t *particles)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort", 1., 0, 0);
  }

  prof_start(pr);
    
  psc_foreach_patch(ppsc, p) {
    sort_patch(p, &particles->p[p]);
  }

  prof_stop(pr);
}

// ======================================================================
// psc_sort: subclass "cuda"

struct psc_sort_ops psc_sort_cuda_ops = {
  .name                  = "cuda",
  .run                   = psc_sort_cuda_run,
};

