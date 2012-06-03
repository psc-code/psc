
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
sort_patch(int p, struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  cuda_find_block_indices_ids(prts, cuda->h_dev->bidx, cuda->h_dev->ids);
  sort_pairs_device(cuda->h_dev->bidx, cuda->h_dev->ids, prts->n_part);
  cuda_reorder_and_offsets(prts, cuda->h_dev->bidx, cuda->h_dev->ids);
}

void
cuda_sort_patch(int p, struct psc_particles *prts)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort_patch", 1., 0, 0);
  }

  prof_start(pr);
  sort_patch(p, prts);
  prof_stop(pr);
}

void
cuda_sort_patch_by_cell(int p, struct psc_particles *prts)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort_patch", 1., 0, 0);
  }

  prof_start(pr);
  sort_patch_by_cell(p, prts);
  prof_stop(pr);
}

static void
psc_sort_cuda_run(struct psc_sort *sort, struct psc_particles *prts_base)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, "cuda", 0);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_sort", 1., 0, 0);
  }

  prof_start(pr);
  sort_patch(prts->p, prts);
  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, 0);
}

// ======================================================================
// psc_sort: subclass "cuda"

struct psc_sort_ops psc_sort_cuda_ops = {
  .name                  = "cuda",
  .run                   = psc_sort_cuda_run,
};

