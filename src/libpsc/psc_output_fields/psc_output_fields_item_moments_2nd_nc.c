
#include "psc_output_fields_item_private.h"

#include <math.h>

#include "common_moments.c"

// ======================================================================
// boundary stuff FIXME, should go elsewhere, and FIXME, duplicated, kinda,
// from cell-centered 1st order

static void
add_ghosts_reflecting_lo(struct psc_fields *pf, int d, int mb, int me)
{
  struct psc_patch *patch = ppsc->patch + pf->p;

  if (d == 1) {
    for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = 0; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy+1,iz) += F3(pf, m, ix,iy-1,iz);
	    F3(pf, m, ix,iy-1,iz) = 0.;
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0; iy < patch->ldims[1]; iy++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iz = 0; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz+1) += F3(pf, m, ix,iy,iz-1);
	    F3(pf, m, ix,iy,iz-1) = 0.;
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
    for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = patch->ldims[1]; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy-1,iz) += F3(pf, m, ix,iy+1,iz);
	    F3(pf, m, ix,iy+1,iz) = 0.;
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0; iy < patch->ldims[1]; iy++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iz = patch->ldims[2]; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz-1) += F3(pf, m, ix,iy,iz+1);
	    F3(pf, m, ix,iy,iz+1) += 0.;
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
// n

static void
do_n_run(int p, fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    DEPOSIT_TO_GRID_2ND_NC(part, pf, m, 1.f);
  }
}

static void
n_run(struct psc_output_fields_item *item, struct psc_fields *flds,
      struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_n_run(res->p, res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  add_ghosts_boundary(res, 0, res->nr_comp);
}

// ======================================================================
// rho

static void
do_rho_run(int p, fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    DEPOSIT_TO_GRID_2ND_NC(part, pf, 0, ppsc->kinds[m].q);
  }
}

static void
rho_run(struct psc_output_fields_item *item, struct psc_fields *flds,
      struct psc_particles *prts_base, struct psc_fields *res_base)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  struct psc_fields *res = psc_fields_get_as(res_base, FIELDS_TYPE, 0, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_rho_run(res->p, res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  psc_fields_put_as(res, res_base, 0, 1);
  add_ghosts_boundary(res, 0, 1);
}

