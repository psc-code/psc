
#include "psc_sort_private.h"
#include "psc_particles_as_c.h"

#include "psc_glue.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <string.h>

// ======================================================================
// quicksort

static void
find_cell_indices(int p, particles_t *pp)
{
}

static inline int
get_cell_index(int p, const particle_t *part)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t dxi = 1.f / ppsc->dx[0];
  particle_real_t dyi = 1.f / ppsc->dx[1];
  particle_real_t dzi = 1.f / ppsc->dx[2];
  int *ldims = patch->ldims;
  int *ibn = ppsc->ibn;
  
  particle_real_t u = (part->xi - patch->xb[0]) * dxi;
  particle_real_t v = (part->yi - patch->xb[1]) * dyi;
  particle_real_t w = (part->zi - patch->xb[2]) * dzi;
  int j0 = particle_real_nint(u) + ibn[0];
  int j1 = particle_real_nint(v) + ibn[1];
  int j2 = particle_real_nint(w) + ibn[2];
    
  return ((j2) * (ldims[1] + 2*ibn[1]) + j1) * (ldims[0] + 2*ibn[0]) + j0;
}

static inline int
get_cell_index_2x2x2(int p, const particle_t *part)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t dxi = 1.f / ppsc->dx[0];
  particle_real_t dyi = 1.f / ppsc->dx[1];
  particle_real_t dzi = 1.f / ppsc->dx[2];
  int *ldims = patch->ldims;
  int ibn[3] = { 2, 2, 2 }; // must be divisible by 2
  
  particle_real_t u = (part->xi - patch->xb[0]) * dxi;
  particle_real_t v = (part->yi - patch->xb[1]) * dyi;
  particle_real_t w = (part->zi - patch->xb[2]) * dzi;
  int j0 = particle_real_nint(u) + ibn[0];
  int j1 = particle_real_nint(v) + ibn[1];
  int j2 = particle_real_nint(w) + ibn[2];
    
  return ((((j2 >> 1) * (ldims[1] + 2*ibn[1]) +
	   (j1 >> 1)) * (ldims[0] + 2*ibn[0]) +
	   (j0 >> 1)) << 0);// | ((j2 & 1) << 2) | ((j1 & 1) << 1) | (j0 & 1);
}

// ======================================================================

int
cell_map_init(struct cell_map *map, const int dims[3], const int blocksize[3])
{
  map->b_bits_max = 0;
  for (int d = 0; d < 3; d++) {
    assert(dims[d] % blocksize[d] == 0);
    map->dims[d] = dims[d] / blocksize[d];
    int bits = 0;
    while (blocksize[d] > (1 << bits))
      bits++;
    assert(blocksize[d] == (1 << bits));
    map->b_bits[d] = bits;
    if (bits > map->b_bits_max) {
      map->b_bits_max = bits;
    }
  }
  map->N = dims[0] * dims[1] * dims[2];
  return map->N;
}

void
cell_map_free(struct cell_map *map)
{
}

int
cell_map_3to1(struct cell_map *map, int i[3])
{
  // FIXME, don't change i[]
  for (int d = 0; d < 3; d++) {
    assert(i[d] < (map->dims[d] << map->b_bits[d]));
  }
  int cnt = 0;
  unsigned int idx = 0;
  for (int b = 0; b < map->b_bits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b < map->b_bits[d]) {
	if (i[d] & 1) {
	  idx |= (1 << cnt);
	}
	i[d] >>= 1;
	cnt++;
      }
    }
  }
  idx |= (((i[2]) * (map->dims[1]) + i[1]) * map->dims[0] + i[0]) << cnt;
#if 0
  if (idx >= map->N) {
    printf("idx %d N %d %d %d %d dims %d %d %d\n", idx, map->N, i[0], i[1], i[2],
	   map->dims[0], map->dims[1], map->dims[2]);
  }
#endif
  assert(idx < map->N);
  return idx;
}

void
cell_map_1to3(struct cell_map *map, int idx, int i[3])
{
  for (int d = 0; d < 3; d++) {
    i[d] = 0;
  }
  for (int b = 0; b < map->b_bits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b < map->b_bits[d]) {
	if (idx & 1) {
	  i[d] |= (1 << b);
	}
	idx >>= 1;
      }
    }
  }

  i[0] |= (idx % map->dims[0]) << map->b_bits[0];
  idx /= map->dims[0];
  i[1] |= (idx % map->dims[1]) << map->b_bits[1];
  idx /= map->dims[1];
  i[2] |= idx << map->b_bits[2];
  for (int d = 0; d < 3; d++) {
    assert(i[d] < (map->dims[d] << map->b_bits[d]));
  }
}


static int
compare(const void *_a, const void *_b)
{
  const particle_t *a = _a, *b = _b;

  if (get_cell_index(0, a) < get_cell_index(0, b)) {
    return -1;
  } else if (get_cell_index(0, a) == get_cell_index(0, b)) {
    return 0;
  } else {
    return 1;
  }
}

static void
psc_sort_qsort_run(struct psc_sort *sort, mparticles_base_t *particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("qsort_sort", 1., 0, 0);
  }
  mparticles_t *particles = psc_mparticles_base_get_cf(particles_base);

  prof_start(pr);
  assert(ppsc->nr_patches == 1);
  psc_foreach_patch(ppsc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, 0);
    find_cell_indices(p, pp);
    qsort(pp->particles, pp->n_part, sizeof(*pp->particles), compare);
  }
  prof_stop(pr);

  psc_mparticles_base_put_cf(particles, particles_base);
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
psc_sort_countsort_run(struct psc_sort *sort, mparticles_base_t *particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("countsort_sort", 1., 0, 0);
  }

  mparticles_t *particles = psc_mparticles_base_get_cf(particles_base);

  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particles_t *pp = psc_mparticles_get_patch(particles, 0);
    find_cell_indices(p, pp);
    
    int N = 1;
    for (int d = 0; d < 3; d++) {
      N *= patch->ldims[d] + 2 * ppsc->ibn[d];
    }
    unsigned int *cnts = malloc(N * sizeof(*cnts));
    memset(cnts, 0, N * sizeof(*cnts));
    
    // count
    for (int i = 0; i < pp->n_part; i++) {
      unsigned int cni = get_cell_index(p, &pp->particles[i]);
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
    particle_t *particles2 = malloc(pp->n_part * sizeof(*particles2));
    for (int i = 0; i < pp->n_part; i++) {
      unsigned int cni = get_cell_index(0, &pp->particles[i]);
      memcpy(&particles2[cnts[cni]], &pp->particles[i], sizeof(*particles2));
      cnts[cni]++;
    }
    
    // back to in-place
    memcpy(pp->particles, particles2, pp->n_part * sizeof(*particles2));
    
    free(particles2);
    free(cnts);
  }

  prof_stop(pr);

  psc_mparticles_base_put_cf(particles, particles_base);
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

struct psc_sort_countsort2 {
  int blocksize[3];
  int mask;
};

#define VAR(x) (void *)offsetof(struct psc_sort_countsort2, x)
static struct param psc_sort_countsort2_descr[] = {
  { "blocksize"     , VAR(blocksize)       , PARAM_INT3(1, 1, 1)   },
  { "mask"          , VAR(mask)            , PARAM_INT(0)          },
  {},
};
#undef VAR


static void
psc_sort_countsort2_run(struct psc_sort *sort, mparticles_base_t *particles_base)
{
  struct psc_sort_countsort2 *cs2 = mrc_to_subobj(sort, struct psc_sort_countsort2);

  static int pr;
  if (!pr) {
    pr = prof_register("countsort2_sort", 1., 0, 0);
  }

  mparticles_t *particles = psc_mparticles_base_get_cf(particles_base);

  prof_start(pr);
  unsigned int mask = cs2->mask;
  psc_foreach_patch(ppsc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, 0);
    struct psc_patch *patch = &ppsc->patch[p];
    struct cell_map map;
    int N = cell_map_init(&map, patch->ldims, cs2->blocksize);
    
    unsigned int *cnis = malloc(pp->n_part * sizeof(*cnis));
    for (int i = 0; i < pp->n_part; i++) {
      particle_t *p = &pp->particles[i];
      particle_real_t dxi = 1.f / ppsc->dx[0];
      particle_real_t dyi = 1.f / ppsc->dx[1];
      particle_real_t dzi = 1.f / ppsc->dx[2];
      particle_real_t xi[3] = {
	(p->xi - patch->xb[0]) * dxi,
	(p->yi - patch->xb[1]) * dyi,
	(p->zi - patch->xb[2]) * dzi };
      int pos[3];
      for (int d = 0; d < 3; d++) {
	pos[d] = particle_real_fint(xi[d]);
#if 0
	if (pos[d] < 0 || pos[d] >= patch->ldims[d]) {
	  printf("i %d d %d pos %d // %d xi %g dxi %g\n",
		 i, d, pos[d], patch->ldims[d],
		 (&p->xi)[d], 1. / ppsc->dx[d]);
	}
#endif
      }

      cnis[i] = cell_map_3to1(&map, pos) & ~mask;
      assert(cnis[i] < N);
    }
    cell_map_free(&map);
    
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
    particle_t *particles2 = malloc(pp->n_part * sizeof(*particles2));
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
    free(cnis);
    free(cnts);
  }

  prof_stop(pr);

  psc_mparticles_base_put_cf(particles, particles_base);
}

// ======================================================================
// psc_sort: subclass "countsort2"

struct psc_sort_ops psc_sort_countsort2_ops = {
  .name                  = "countsort2",
  .size                  = sizeof(struct psc_sort_countsort2),
  .param_descr           = psc_sort_countsort2_descr,
  .run                   = psc_sort_countsort2_run,
};

