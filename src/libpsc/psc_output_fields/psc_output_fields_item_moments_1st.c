
#include "psc_output_fields_item_private.h"
#include "psc_bnd.h"

#include <math.h>

#include "common_moments.c"

// ======================================================================
// boundary stuff FIXME, should go elsewhere...

static void
add_ghosts_reflecting_lo(mfields_c_t *res, int p, int d, int mb, int me)
{
  fields_t *pf = psc_mfields_get_patch(res, p);
  struct psc_patch *patch = ppsc->patch + p;

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
add_ghosts_reflecting_hi(mfields_c_t *res, int p, int d, int mb, int me)
{
  fields_t *pf = psc_mfields_get_patch(res, p);
  struct psc_patch *patch = ppsc->patch + p;

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
add_ghosts_boundary(mfields_c_t *res, int mb, int me)
{
  psc_foreach_patch(ppsc, p) {
    // lo
    for (int d = 0; d < 3; d++) {
      if (ppsc->patch[p].off[d] == 0) {
	if (ppsc->domain.bnd_part_lo[d] == BND_PART_REFLECTING) {
	  add_ghosts_reflecting_lo(res, p, d, mb, me);
	}
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	if (ppsc->domain.bnd_part_hi[d] == BND_PART_REFLECTING) {
	  add_ghosts_reflecting_hi(res, p, d, mb, me);
	}
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
n_1st_run(struct psc_output_fields_item *item, mfields_base_t *flds,
	  mparticles_base_t *particles_base, mfields_c_t *res)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  psc_mfields_zero_range(res, 0, res->nr_fields);
  
  psc_foreach_patch(ppsc, p) {
    do_n_1st_run(p, psc_mfields_get_patch_c(res, p),
		 psc_mparticles_get_patch(particles, p));
  }

  psc_mparticles_put_cf(particles, particles_base, MP_DONT_COPY);

  psc_bnd_add_ghosts(item->bnd, res, 0, res->nr_fields);
  add_ghosts_boundary(res, 0, res->nr_fields);
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
v_1st_run(struct psc_output_fields_item *item, mfields_base_t *flds,
	  mparticles_base_t *particles_base, mfields_c_t *res)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  psc_mfields_zero_range(res, 0, res->nr_fields);
  
  psc_foreach_patch(ppsc, p) {
    do_v_1st_run(p, psc_mfields_get_patch_c(res, p),
		 psc_mparticles_get_patch(particles, p));
  }

  psc_mparticles_put_cf(particles, particles_base, MP_DONT_COPY);

  psc_bnd_add_ghosts(item->bnd, res, 0, res->nr_fields);
  add_ghosts_boundary(res, 0, res->nr_fields);
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
vv_1st_run(struct psc_output_fields_item *item, mfields_base_t *flds,
	  mparticles_base_t *particles_base, mfields_c_t *res)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  psc_mfields_zero_range(res, 0, res->nr_fields);
  
  psc_foreach_patch(ppsc, p) {
    do_vv_1st_run(p, psc_mfields_get_patch_c(res, p),
		  psc_mparticles_get_patch(particles, p));
  }

  psc_mparticles_put_cf(particles, particles_base, MP_DONT_COPY);

  psc_bnd_add_ghosts(item->bnd, res, 0, res->nr_fields);
  add_ghosts_boundary(res, 0, res->nr_fields);
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

