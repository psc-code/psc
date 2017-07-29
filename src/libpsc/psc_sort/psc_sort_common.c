

// ======================================================================
// quicksort

static void
find_cell_indices(int p, struct psc_particles *pp)
{
}

static inline int
get_cell_index(int p, const particle_t *part)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t dxi = 1.f / patch->dx[0];
  particle_real_t dyi = 1.f / patch->dx[1];
  particle_real_t dzi = 1.f / patch->dx[2];
  int *ldims = patch->ldims;
  int *ibn = ppsc->ibn;
  
  particle_real_t u = part->xi * dxi;
  particle_real_t v = part->yi * dyi;
  particle_real_t w = part->zi * dzi;
  int j0 = particle_real_nint(u) + ibn[0];
  int j1 = particle_real_nint(v) + ibn[1];
  int j2 = particle_real_nint(w) + ibn[2];
    
  return ((j2) * (ldims[1] + 2*ibn[1]) + j1) * (ldims[0] + 2*ibn[0]) + j0;
}

static inline int
get_cell_index_2x2x2(int p, const particle_t *part)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t dxi = 1.f / patch->dx[0];
  particle_real_t dyi = 1.f / patch->dx[1];
  particle_real_t dzi = 1.f / patch->dx[2];
  int *ldims = patch->ldims;
  int ibn[3] = { 2, 2, 2 }; // must be divisible by 2
  
  particle_real_t u = part->xi * dxi;
  particle_real_t v = part->yi * dyi;
  particle_real_t w = part->zi * dzi;
  int j0 = particle_real_nint(u) + ibn[0];
  int j1 = particle_real_nint(v) + ibn[1];
  int j2 = particle_real_nint(w) + ibn[2];
    
  return ((((j2 >> 1) * (ldims[1] + 2*ibn[1]) +
	   (j1 >> 1)) * (ldims[0] + 2*ibn[0]) +
	   (j0 >> 1)) << 0);// | ((j2 & 1) << 2) | ((j1 & 1) << 1) | (j0 & 1);
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
psc_sort_qsort_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("qsort_sort", 1., 0, 0);
  }
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    find_cell_indices(p, prts);
    qsort(particles_get_one(prts, 0), prts->n_part,
	  sizeof(*particles_get_one(prts, 0)), compare);
  }
  prof_stop(pr);

  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ======================================================================
// counting sort

static void
psc_sort_countsort_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("countsort_sort", 1., 0, 0);
  }

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_patch *patch = &ppsc->patch[p];
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    find_cell_indices(p, prts);
    
    int N = 1;
    for (int d = 0; d < 3; d++) {
      N *= patch->ldims[d] + 2 * ppsc->ibn[d];
    }
    unsigned int *cnts = malloc(N * sizeof(*cnts));
    memset(cnts, 0, N * sizeof(*cnts));
    
    // count
    for (int i = 0; i < prts->n_part; i++) {
      unsigned int cni = get_cell_index(prts->p, particles_get_one(prts, i));
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
    assert(cur == prts->n_part);
    
    // move into new position
    particle_t *particles2 = malloc(prts->n_part * sizeof(*particles2));
    for (int i = 0; i < prts->n_part; i++) {
      unsigned int cni = get_cell_index(0, particles_get_one(prts, i));
      memcpy(&particles2[cnts[cni]], particles_get_one(prts, i), sizeof(*particles2));
      cnts[cni]++;
    }
    
    // back to in-place
    memcpy(particles_get_one(prts, 0), particles2, prts->n_part * sizeof(*particles2));
    
    free(particles2);
    free(cnts);
  }
  
  prof_stop(pr);

  psc_mparticles_put_as(mprts, mprts_base, 0);
}


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
psc_sort_countsort2_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  struct psc_sort_countsort2 *cs2 = mrc_to_subobj(sort, struct psc_sort_countsort2);

  static int pr;
  if (!pr) {
    pr = prof_register("countsort2_sort", 1., 0, 0);
  }

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_patch *patch = &ppsc->patch[p];
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);

    unsigned int mask = cs2->mask;
    struct cell_map map;
    int N = cell_map_init(&map, patch->ldims, cs2->blocksize);
    
    unsigned int *cnis = malloc(prts->n_part * sizeof(*cnis));
    for (int i = 0; i < prts->n_part; i++) {
      particle_t *p = particles_get_one(prts, i);
      particle_real_t dxi = 1.f / patch->dx[0];
      particle_real_t dyi = 1.f / patch->dx[1];
      particle_real_t dzi = 1.f / patch->dx[2];
      particle_real_t xi[3] = { p->xi * dxi, p->yi * dyi, p->zi * dzi };
      int pos[3];
      for (int d = 0; d < 3; d++) {
	pos[d] = particle_real_fint(xi[d]);
#if 1
	if (pos[d] < 0 || pos[d] >= patch->ldims[d]) {
	  printf("i %d d %d pos %d // %d xi %g dxi %g\n",
		 i, d, pos[d], patch->ldims[d],
		 (&p->xi)[d], 1. / patch->dx[d]);
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
    for (int i = 0; i < prts->n_part; i++) {
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
    assert(cur == prts->n_part);
    
    // move into new position
    particle_t *particles2 = malloc(prts->n_part * sizeof(*particles2));
    for (int i = 0; i < prts->n_part; i++) {
      unsigned int cni = cnis[i];
      int n = 1;
      while (i+n < prts->n_part && cnis[i+n] == cni) {
	n++;
      }
      memcpy(&particles2[cnts[cni]], particles_get_one(prts, i), n * sizeof(*particles2));
      cnts[cni] += n;
      i += n - 1;
    }
    
    // back to in-place
    memcpy(particles_get_one(prts, 0), particles2, prts->n_part * sizeof(*particles2));
    
    free(particles2);
    free(cnis);
    free(cnts);
  }

  prof_stop(pr);

  psc_mparticles_put_as(mprts, mprts_base, 0);
}


