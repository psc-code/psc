
#include "psc_output_fields_item_private.h"

#include <math.h>

#include "common_moments.c"

// ======================================================================
// boundary stuff FIXME, should go elsewhere...

static void
add_ghosts_reflecting_lo(struct psc_fields *pf, int d, int mb, int me)
{
  struct psc_patch *patch = ppsc->patch + pf->p;

  if (d == 1) {
    for (int iz = 0; iz < patch->ldims[2]; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = 0; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy-1,iz);
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0; iy < patch->ldims[1]; iy++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iz = 0; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy,iz-1);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_reflecting_hi(struct psc_fields *pf, int d, int mb, int me)
{
  struct psc_patch *patch = ppsc->patch + pf->p;

  if (d == 1) {
    for (int iz = 0; iz < patch->ldims[2]; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = patch->ldims[1] - 1; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy+1,iz);
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0; iy < patch->ldims[1]; iy++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iz = patch->ldims[2] - 1; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy,iz+1);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_boundary(struct psc_fields *res, int mb, int me)
{
  // lo
  for (int d = 0; d < 3; d++) {
    if (ppsc->patch[res->p].off[d] == 0) {
      if (ppsc->domain.bnd_part_lo[d] == BND_PART_REFLECTING) {
	add_ghosts_reflecting_lo(res, d, mb, me);
      }
    }
  }
  // hi
  for (int d = 0; d < 3; d++) {
    if (ppsc->patch[res->p].off[d] + ppsc->patch[res->p].ldims[d] == ppsc->domain.gdims[d]) {
      if (ppsc->domain.bnd_part_hi[d] == BND_PART_REFLECTING) {
	add_ghosts_reflecting_hi(res, d, mb, me);
      }
    }
  }
}

// ======================================================================
// n_1st

static void
do_n_1st_run(int p, fields_t *pf, struct psc_particles *prts)
{
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / ppsc->dx[0], dyi = 1.f / ppsc->dx[1], dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    DEPOSIT_TO_GRID_1ST_CC(part, pf, m, 1.f);
  }
}

static void
n_1st_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	  struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_particles_reorder(prts); // FIXME
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_n_1st_run(res->p, res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  add_ghosts_boundary(res, 0, res->nr_comp);
}

static int
n_1st_get_nr_components(struct psc_output_fields_item *item)
{
  return ppsc->nr_kinds;
}

static const char *
n_1st_get_component_name(struct psc_output_fields_item *item, int m)
{
  static char s[100];
  sprintf(s, "n_%s", ppsc->kinds[m].name);
  return s;
}

// ======================================================================
// v_1st

static void
do_v_1st_run(int p, fields_t *pf, struct psc_particles *prts)
{
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / ppsc->dx[0], dyi = 1.f / ppsc->dx[1], dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts
, n);
    int mm = particle_kind(part) * 3;

    particle_real_t vxi[3];
    particle_calc_vxi(part, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + m, vxi[m]);
    }
  }
}

static void
v_1st_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	  struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_v_1st_run(res->p, res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  add_ghosts_boundary(res, 0, res->nr_comp);
}

static int
v_1st_get_nr_components(struct psc_output_fields_item *item)
{
  return 3 * ppsc->nr_kinds;
}

static const char *
v_1st_get_component_name(struct psc_output_fields_item *item, int m)
{
  static char s[100];
  sprintf(s, "v%c_%s", 'x' + m % 3, ppsc->kinds[m / 3].name);
  return s;
}

// ======================================================================
// p_1st

static void
do_p_1st_run(int p, fields_t *pf, struct psc_particles *prts)
{
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / ppsc->dx[0], dyi = 1.f / ppsc->dx[1], dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts
, n);
    int mm = particle_kind(part) * 3;
    particle_real_t *pxi = &part->pxi;

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + m, particle_mni(part) * pxi[m]);
    }
  }
}

static void
p_1st_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	  struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_p_1st_run(res->p, res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  add_ghosts_boundary(res, 0, res->nr_comp);
}

static int
p_1st_get_nr_components(struct psc_output_fields_item *item)
{
  return 3 * ppsc->nr_kinds;
}

static const char *
p_1st_get_component_name(struct psc_output_fields_item *item, int m)
{
  static char s[100];
  sprintf(s, "p%c_%s", 'x' + m % 3, ppsc->kinds[m / 3].name);
  return s;
}

// ======================================================================
// vv_1st

static void
do_vv_1st_run(int p, fields_t *pf, struct psc_particles *prts)
{
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / ppsc->dx[0], dyi = 1.f / ppsc->dx[1], dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int mm = particle_kind(part) * 3;

    particle_real_t vxi[3];
    particle_calc_vxi(part, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + m, vxi[m] * vxi[m]);
    }
  }
}

static void
vv_1st_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	   struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_vv_1st_run(res->p, res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  add_ghosts_boundary(res, 0, res->nr_comp);
}

static int
vv_1st_get_nr_components(struct psc_output_fields_item *item)
{
  return 3 * ppsc->nr_kinds;
}

static const char *
vv_1st_get_component_name(struct psc_output_fields_item *item, int m)
{
  static char s[100];
  sprintf(s, "v%cv%c_%s", 'x' + m % 3, 'x' + m % 3, ppsc->kinds[m / 3].name);
  return s;
}

#define MAKE_POFI_OPS(WHAT, TYPE)					\
struct psc_output_fields_item_ops psc_output_fields_item_##WHAT##_##TYPE##_ops = { \
  .name               = #WHAT "_" #TYPE,				\
  .get_component_name = WHAT##_get_component_name,			\
  .get_nr_components  = WHAT##_get_nr_components,			\
  .run                = WHAT##_run,					\
  .flags              = POFI_ADD_GHOSTS,				\
}

