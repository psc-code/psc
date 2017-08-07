
#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#if PSC_PARTICLES_AS_DOUBLE

#define PFX(x) psc_particles_double_ ## x
#define psc_particles_sub psc_particles_double

#elif PSC_PARTICLES_AS_SINGLE

#define PFX(x) psc_particles_single_ ## x
#define psc_particles_sub psc_particles_single

#elif PSC_PARTICLES_AS_SINGLE_BY_BLOCK

#define PFX(x) psc_particles_single_by_block_ ## x
#define psc_particles_sub psc_particles_single_by_block

#elif PSC_PARTICLES_AS_C

#define PFX(x) psc_particles_c_ ## x
#define psc_particles_sub psc_particles_c

#elif PSC_PARTICLES_AS_FORTRAN

#define PFX(x) psc_particles_fortran_ ## x
#define psc_particles_sub psc_particles_fortran

#endif

// ======================================================================
// psc_particles "single" / "double" / "c"

// ----------------------------------------------------------------------
// psc_particles_sub_setup

static void
PFX(setup)(struct psc_particles *prts)
{
  struct psc_particles_sub *sub = psc_particles_sub(prts);

  int n_alloced = psc_particles_size(prts) * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  sub->particles = calloc(n_alloced, sizeof(*sub->particles));

#if PSC_PARTICLES_AS_SINGLE
  sub->particles_alt = calloc(n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  sub->particles_alt = calloc(n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
  sub->b_off = calloc(sub->nr_blocks + 2, sizeof(*sub->b_off));
#endif
}

// ----------------------------------------------------------------------
// psc_particles_sub_destroy

static void
PFX(destroy)(struct psc_particles *prts)
{
  struct psc_particles_sub *sub = psc_particles_sub(prts);

  free(sub->particles);

#if PSC_PARTICLES_AS_SINGLE
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
  free(sub->b_off);
#endif
}

// ----------------------------------------------------------------------
// psc_particles_sub_realloc

void
PFX(realloc)(struct psc_particles *prts, int new_n_part)
{
  struct psc_particles_sub *sub = psc_particles_sub(prts);

  if (new_n_part <= psc_mparticles_n_alloced(prts->mprts, prts->p))
    return;

  int n_alloced = new_n_part * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  sub->particles = realloc(sub->particles, n_alloced * sizeof(*sub->particles));

#if PSC_PARTICLES_AS_SINGLE
  sub->b_idx = realloc(sub->b_idx, n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(n_alloced * sizeof(*sub->particles_alt));
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  sub->b_idx = realloc(sub->b_idx, n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(n_alloced * sizeof(*sub->particles_alt));
#endif
}

// ----------------------------------------------------------------------
// psc_particles: subclass ops

struct psc_particles_ops PFX(ops) = {
  .name                    = PARTICLE_TYPE,
  .size                    = sizeof(struct psc_particles_sub),
  .setup                   = PFX(setup),
  .destroy                 = PFX(destroy),
};

