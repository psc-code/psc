
#include "psc_sort_private.h"

#include <mrc_profile.h>

// FIXME -> header
EXTERN_C void sort_pairs_host(unsigned int *keys, unsigned int *vals, int n);
EXTERN_C void create_indices_host(unsigned int *cnis, struct cell_map *map,
				  particles_cuda_t *pp, struct psc_patch *patch);
EXTERN_C void particles_cuda_copy_to_device(particles_cuda_t *pp);


// ======================================================================
// cuda sort

#if 0
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
#endif

void
find_cell_indices_host(particles_cuda_t *pp, struct psc_patch *patch,
		       unsigned int *cnis, unsigned int *ids)
{
  struct cell_map map;

  cell_map_init(&map, patch->ldims, (int[3]) { 1, 8, 8 });
  create_indices_host(cnis, &map, pp, patch);
  cell_map_free(&map);
  for (int i = 0; i < pp->n_part; i++) {
    ids[i] = i;
  }
}

static void
sort_patch(int p, particles_cuda_t *pp)
{
  struct psc_patch *patch = &ppsc->patch[p];

  unsigned int *cnis = malloc(pp->n_part * sizeof(*cnis));
  unsigned int *ids = malloc(pp->n_part * sizeof(*ids));

  find_cell_indices_host(pp, patch, cnis, ids);

  sort_pairs_host(cnis, ids, pp->n_part);

  // move into new position
  float4 *xi4 = malloc(pp->n_part * sizeof(*xi4));
  float4 *pxi4 = malloc(pp->n_part * sizeof(*pxi4));
  for (int i = 0; i < pp->n_part; i++) {
    xi4[i] = pp->h_part.xi4[ids[i]];
    pxi4[i] = pp->h_part.pxi4[ids[i]];
  }
  free(cnis);
  free(ids);

  // back to in-place
  memcpy(pp->h_part.xi4, xi4, pp->n_part * sizeof(*xi4));
  memcpy(pp->h_part.pxi4, pxi4, pp->n_part * sizeof(*pxi4));
  
  free(xi4);
  free(pxi4);

  particles_cuda_copy_to_device(pp); 
}

static void
psc_sort_cuda_run(struct psc_sort *sort, mparticles_base_t *particles_base)
{
  mparticles_cuda_t particles;
  psc_mparticles_cuda_get_from(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort", 1., 0, 0);
  }

  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    sort_patch(p, &particles.p[p]);
  }

  prof_stop(pr);
  psc_mparticles_cuda_put_to(&particles, particles_base);
}

// ======================================================================
// psc_sort: subclass "cuda"

struct psc_sort_ops psc_sort_cuda_ops = {
  .name                  = "cuda",
  .run                   = psc_sort_cuda_run,
};

