
#include "psc_collision_private.h"

#include <psc_fields_c.h>
#include <mrc_profile.h>
#include <mrc_params.h>

struct psc_collision_c {
  // parameters
  int every;
  double nudt0;

  // internal
  struct psc_mfields *mflds;
};

#define psc_collision_c(o) mrc_to_subobj(o, struct psc_collision_c)

#define VAR(x) (void *)offsetof(struct psc_collision_c, x)
static struct param psc_collision_c_descr[] = {
  { "every"         , VAR(every)       , PARAM_INT(1)     },
  { "nudt0"         , VAR(nudt0)       , PARAM_DOUBLE(-1.) },
  {},
};
#undef VAR

enum {
  STATS_MIN,
  STATS_MAX,
  STATS_MED,
  STATS_NLARGE,
  STATS_NCOLL,
  NR_STATS,
};

struct psc_collision_stats {
  particle_c_real_t s[NR_STATS];
};


// ----------------------------------------------------------------------
// calc_stats

static int
compare(const void *_a, const void *_b)
{
  const particle_c_real_t *a = _a, *b = _b;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return 1;
  } else {
    return 0;
  }
}

static void
calc_stats(struct psc_collision_stats *stats, particle_c_real_t *nudts, int cnt)
{
  stats->s[STATS_NLARGE] = 0;
  for (int n = cnt - 1; n >= 0; n--) {
    if (nudts[n] < 1.) {
      break;
    }
    stats->s[STATS_NLARGE]++;
  }
  qsort(nudts, cnt, sizeof(*nudts), compare);
  stats->s[STATS_MIN] = nudts[0];
  stats->s[STATS_MAX] = nudts[cnt-1];
  stats->s[STATS_MED] = nudts[cnt/2];
  stats->s[STATS_NCOLL] = cnt;
}

// ----------------------------------------------------------------------
// find_cell_index

static inline int
find_cell_index(particle_c_t *prt, particle_c_real_t *dxi, int ldims[3])
{
  int pos[3];
  particle_c_real_t *xi = &prt->xi;
  for (int d = 0; d < 3; d++) {
    pos[d] = particle_c_real_fint(xi[d] * dxi[d]);
    assert(pos[d] >= 0 && pos[d] < ldims[d]);
  }
  return (pos[2] * ldims[1] + pos[1]) * ldims[0] + pos[0];
}

// ----------------------------------------------------------------------
// find_cell_offsets

static void
find_cell_offsets(int offsets[], struct psc_particles *prts)
{
  particle_c_real_t dxi[3] = { 1.f / ppsc->dx[0], 1.f / ppsc->dx[1], 1.f / ppsc->dx[2] };
  int *ldims = ppsc->patch[prts->p].ldims;
  int last = 0;
  offsets[last] = 0;
  for (int n = 0; n < prts->n_part; n++) {
    particle_c_t *prt = particles_c_get_one(prts, n);
    int cell_index = find_cell_index(prt, dxi, ldims);
    assert(cell_index >= last);
    while (last < cell_index) {
      offsets[++last] = n;
    }
  }
  while (last < ldims[0] * ldims[1] * ldims[2]) {
    offsets[++last] = prts->n_part;
  }
}

// ----------------------------------------------------------------------
// randomize_in_cell

static void
randomize_in_cell(struct psc_particles *prts, int n_start, int n_end)
{
  int nn = n_end - n_start;
  for (int n = 0; n < nn - 1; n++) {
    int n_partner = random() % (nn - n);
    if (n != n_partner) {
      // swap n, n_partner
      particle_c_t tmp = *particles_c_get_one(prts, n_start + n);
      *particles_c_get_one(prts, n_start + n) = *particles_c_get_one(prts, n_start + n_partner);    
      *particles_c_get_one(prts, n_start + n_partner) = tmp;
    }
  }
}

// ----------------------------------------------------------------------
// bc

static particle_c_real_t
bc(struct psc_particles *prts, int nudt1, int n1, int n2)
{
  return 99.;
}

// ----------------------------------------------------------------------
// collide_in_cell

static void
collide_in_cell(struct psc_collision *collision,
		struct psc_particles *prts, int n_start, int n_end,
		struct psc_collision_stats *stats)
{
  struct psc_collision_c *coll = psc_collision_c(collision);

  int nn = n_end - n_start;
  
  int n = 0;
  if (nn < 2) { // can't collide only one (or zero) particles
    return;
  }

  // all particles need to have same weight!
  particle_c_real_t wni = particles_c_get_one(prts, n_start)->wni;
  particle_c_real_t nudt1 = wni * nn * coll->every * coll->nudt0;

  particle_c_real_t *nudts = malloc((nn / 2 + 2) * sizeof(*nudts));
  int cnt = 0;

  if (nn % 2 == 1) { // odd # of particles: do 3-collision
    nudts[cnt++] = bc(prts, .5 * nudt1, n_start    , n_start + 1);
    nudts[cnt++] = bc(prts, .5 * nudt1, n_start    , n_start + 2);
    nudts[cnt++] = bc(prts, .5 * nudt1, n_start + 1, n_start + 2);
    n = 3;
  }
  for (; n < nn;  n += 2) { // do remaining particles as pair
    nudts[cnt++] = bc(prts, nudt1, n_start + n, n_start + n + 1);
  }

  calc_stats(stats, nudts, cnt);
}

// ----------------------------------------------------------------------
// psc_collision_c_setup

static void
psc_collision_c_setup(struct psc_collision *collision)
{
  struct psc_collision_c *coll = psc_collision_c(collision);

  coll->mflds = psc_mfields_create(psc_collision_comm(collision));
  psc_mfields_set_type(coll->mflds, "c");
  psc_mfields_set_domain(coll->mflds, ppsc->mrc_domain);
  psc_mfields_set_param_int(coll->mflds, "nr_fields", 5);
  psc_mfields_set_param_int3(coll->mflds, "ibn", ppsc->ibn);
  psc_mfields_setup(coll->mflds);
  psc_mfields_set_comp_name(coll->mflds, 0, "coll_nudt_min");
  psc_mfields_set_comp_name(coll->mflds, 1, "coll_nudt_med");
  psc_mfields_set_comp_name(coll->mflds, 2, "coll_nudt_max");
  psc_mfields_set_comp_name(coll->mflds, 3, "coll_nudt_nlarge");
  psc_mfields_set_comp_name(coll->mflds, 4, "coll_nudt_ncoll");
  // FIXME, needs to be registered for rebalancing
}

// ----------------------------------------------------------------------
// psc_collision_c_destroy

static void
psc_collision_c_destroy(struct psc_collision *collision)
{
  struct psc_collision_c *coll = psc_collision_c(collision);

  psc_mfields_destroy(coll->mflds);
}

// ----------------------------------------------------------------------
// psc_collision_c_run

static void
psc_collision_c_run(struct psc_collision *collision,
		    struct psc_particles *prts_base)
{
  struct psc_collision_c *coll = psc_collision_c(collision);

  static int pr;
  if (!pr) {
    pr = prof_register("collision", 1., 0, 0);
  }

  assert(coll->nudt0 > 0.);

  struct psc_particles *prts = psc_particles_get_as(prts_base, "c", 0);

  if (ppsc->timestep % coll->every != 0) {
    return;
  }

  prof_start(pr);
  
  int *ldims = ppsc->patch[prts->p].ldims;
  int nr_cells = ldims[0] * ldims[1] * ldims[2];
  int *offsets = calloc(nr_cells + 1, sizeof(*offsets));
  struct psc_collision_stats stats_total = {};
  
  find_cell_offsets(offsets, prts);

  struct psc_fields *flds = psc_mfields_get_patch(coll->mflds, prts->p);
  psc_foreach_3d(ppsc, prts->p, ix, iy, iz, 0, 0) {
    int c = (iz * ldims[1] + iy) * ldims[0] + ix;
    randomize_in_cell(prts, offsets[c], offsets[c+1]);

    struct psc_collision_stats stats = {};
    collide_in_cell(collision, prts, offsets[c], offsets[c+1], &stats);

    for (int s = 0; s < NR_STATS; s++) {
      F3_C(flds, s, ix,iy,iz) = stats.s[s];
      stats_total.s[s] += stats.s[s];
    }
  } psc_foreach_3d_end;

  mprintf("p%d: min %g med %g max %g nlarge %g ncoll %g\n", prts->p,
	  stats_total.s[0] / nr_cells,
	  stats_total.s[1] / nr_cells,
	  stats_total.s[2] / nr_cells,
	  stats_total.s[3] / nr_cells,
	  stats_total.s[4] / nr_cells);

  free(offsets);

  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, 0);
}

// ======================================================================
// psc_collision: subclass "c"

struct psc_collision_ops psc_collision_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_collision_c),
  .param_descr           = psc_collision_c_descr,
  .setup                 = psc_collision_c_setup,
  .destroy               = psc_collision_c_destroy,
  .run                   = psc_collision_c_run,
};
