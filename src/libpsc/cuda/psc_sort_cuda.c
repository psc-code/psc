
#include "psc_sort_private.h"
#include "psc_cuda.h"

#include <mrc_profile.h>

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

static void
sort_patch(int p, particles_cuda_t *pp)
{
  cuda_find_block_indices_ids(pp, pp->d_part.bidx, pp->d_part.ids);
  sort_pairs_device(pp->d_part.bidx, pp->d_part.ids, pp->n_part);
  cuda_reorder_and_offsets(pp, pp->d_part.bidx, pp->d_part.ids);
}

void
cuda_sort_patch(int p, particles_cuda_t *pp)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort_patch", 1., 0, 0);
  }

  prof_start(pr);
  sort_patch(p, pp);
  prof_stop(pr);
}

void
cuda_sort_patch_by_cell(int p, particles_cuda_t *pp)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort_patch", 1., 0, 0);
  }

  prof_start(pr);
  sort_patch_by_cell(p, pp);
  prof_stop(pr);
}

static void
psc_sort_cuda_run(struct psc_sort *sort, mparticles_base_t *particles_base)
{
  mparticles_cuda_t *particles = psc_mparticles_get_cuda(particles_base, 0);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort", 1., 0, 0);
  }

  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    sort_patch(p, psc_mparticles_get_patch_cuda(particles, p));
  }
  prof_stop(pr);

  psc_mparticles_put_cuda(particles, particles_base, 0);
}

// ======================================================================
// psc_sort: subclass "cuda"

struct psc_sort_ops psc_sort_cuda_ops = {
  .name                  = "cuda",
  .run                   = psc_sort_cuda_run,
};

