
#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#if PSC_PARTICLES_AS_DOUBLE

#define PFX(x) psc_particles_double_ ## x

#elif PSC_PARTICLES_AS_SINGLE

#define PFX(x) psc_particles_single_ ## x

#endif

// ======================================================================
// psc_particles "single" / "double"

// ----------------------------------------------------------------------
// psc_particles_sub_setup

#if PSC_PARTICLES_AS_DOUBLE

static void
PFX(setup)(struct psc_particles *prts)
{
  struct psc_particles_double *sub = psc_particles_double(prts);

  int n_alloced = psc_particles_size(prts) * 1.2 + 1000000;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  sub->particles = calloc(n_alloced, sizeof(*sub->particles));
}

#elif PSC_PARTICLES_AS_SINGLE

static void
PFX(setup)(struct psc_particles *prts)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  int n_alloced = psc_particles_size(prts) * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  sub->particles = calloc(n_alloced, sizeof(*sub->particles));
  sub->particles_alt = calloc(n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
}

#endif

// ----------------------------------------------------------------------
// psc_particles_sub_destroy

#if PSC_PARTICLES_AS_DOUBLE

static void
PFX(destroy)(struct psc_particles *prts)
{
  struct psc_particles_double *sub = psc_particles_double(prts);

  free(sub->particles);
}

#elif PSC_PARTICLES_AS_SINGLE

static void
PFX(destroy)(struct psc_particles *prts)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  free(sub->particles);
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
}

#endif
