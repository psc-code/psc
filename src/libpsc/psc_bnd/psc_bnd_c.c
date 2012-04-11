
#include "psc_bnd_c.h"
#include "psc_bnd_private.h"
#include "ddc_particles.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include <mrc_domain.h>
#include <mrc_profile.h>
#include <string.h>

struct psc_bnd_c {
  struct mrc_ddc *ddc;
  struct ddc_particles *ddcp;
  struct ddc_particles *ddcp_photons;
};

#define to_psc_bnd_c(bnd) ((struct psc_bnd_c *)((bnd)->obj.subctx))

static void
ddcp_particles_realloc(void *_particles, int p, int new_n_particles)
{
  mparticles_t *particles = _particles;
  particles_t *pp = psc_mparticles_get_patch(particles, p);
  particles_realloc(pp, new_n_particles);
}

static void *
ddcp_particles_get_addr(void *_particles, int p, int n)
{
  mparticles_t *particles = _particles;
  particles_t *pp = psc_mparticles_get_patch(particles, p);
  return &pp->particles[n];
}

static void
ddcp_photons_realloc(void *_particles, int p, int new_n_particles)
{
  mphotons_t *particles = _particles;
  photons_t *pp = &particles->p[p];
  photons_realloc(pp, new_n_particles);
}

static void *
ddcp_photons_get_addr(void *_particles, int p, int n)
{
  mphotons_t *mphotons = _particles;
  photons_t *photons = &mphotons->p[p];
  return &photons->photons[n];
}

// ----------------------------------------------------------------------
// psc_bnd_c_setup

static void
psc_bnd_c_setup(struct psc_bnd *bnd)
{
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);
  struct psc *psc = bnd->psc;

  bnd_c->ddc = psc_bnd_lib_create_ddc(psc);

  bnd_c->ddcp = ddc_particles_create(bnd_c->ddc, sizeof(particle_t),
				     sizeof(particle_real_t),
				     MPI_PARTICLES_REAL,
				     ddcp_particles_realloc,
				     ddcp_particles_get_addr);

  bnd_c->ddcp_photons = ddc_particles_create(bnd_c->ddc, sizeof(photon_t),
					     sizeof(photon_real_t),
					     MPI_PHOTONS_REAL,
					     ddcp_photons_realloc,
					     ddcp_photons_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_c_unsetup

static void
psc_bnd_c_unsetup(struct psc_bnd *bnd)
{
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);

  mrc_ddc_destroy(bnd_c->ddc);
  ddc_particles_destroy(bnd_c->ddcp);
  ddc_particles_destroy(bnd_c->ddcp_photons);
}

// ----------------------------------------------------------------------
// psc_bnd_c_destroy

static void
psc_bnd_c_destroy(struct psc_bnd *bnd)
{
  psc_bnd_c_unsetup(bnd);
}

// ----------------------------------------------------------------------
// check_domain
//
// check if the underlying mrc_domain changed since setup(),
// which might happen, e.g., through rebalancing.
// In this case, do setup() over.

static void
check_domain(struct psc_bnd *bnd)
{
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);
  struct psc *psc = bnd->psc;

  struct mrc_domain *domain = mrc_ddc_get_domain(bnd_c->ddc);
  if (domain != psc->mrc_domain) {
    psc_bnd_c_unsetup(bnd);
    psc_bnd_c_setup(bnd);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_c_add_ghosts

static void
psc_bnd_c_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);
  check_domain(bnd);

  psc_bnd_lib_add_ghosts(bnd_c->ddc, flds_base, mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_c_fill_ghosts

static void
psc_bnd_c_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);
  check_domain(bnd);

  psc_bnd_lib_fill_ghosts(bnd_c->ddc, flds_base, mb, me);
}

// ----------------------------------------------------------------------
// calc_domain_bounds
//
// calculate bounds of local patch, and global domain

static void
calc_domain_bounds(struct psc *psc, int p, double xb[3], double xe[3],
		   double xgb[3], double xge[3], double xgl[3])
{
  struct psc_patch *psc_patch = &psc->patch[p];

  for (int d = 0; d < 3; d++) {
    xb[d] = psc_patch->off[d] * psc->dx[d];
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_UPML) {
      xgb[d] = psc->pml.size * psc->dx[d];
    } else {
      xgb[d] = 0.;
    }
    
    xe[d] = (psc_patch->off[d] + psc_patch->ldims[d]) * psc->dx[d];
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_UPML) {
      xge[d] = (psc->domain.gdims[d] - psc->pml.size) * psc->dx[d];
    } else {
      xge[d] = psc->domain.gdims[d] * psc->dx[d];
    }
    
    xgl[d] = xge[d] - xgb[d];
  }
  for (int d = 0; d < 3; d++) {
    xb[d]  += ppsc->domain.corner[d] / ppsc->coeff.ld;
    xe[d]  += ppsc->domain.corner[d] / ppsc->coeff.ld;
    xgb[d] += ppsc->domain.corner[d] / ppsc->coeff.ld;
    xge[d] += ppsc->domain.corner[d] / ppsc->coeff.ld;
  }
}

// ----------------------------------------------------------------------
// psc_bnd_c_exchange_particles

static void
psc_bnd_c_exchange_particles(struct psc_bnd *bnd, mparticles_base_t *particles_base)
{
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);
  struct psc *psc = bnd->psc;
  check_domain(bnd);

  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  static int pr_A, pr_B;
  if (!pr_A) {
    pr_A = prof_register("c_xchg_part_prep", 1., 0, 0);
    pr_B = prof_register("c_xchg_part_comm", 1., 0, 0);
  }
  prof_start(pr_A);

  struct ddc_particles *ddcp = bnd_c->ddcp;

  double xb[3], xe[3], xgb[3], xge[3], xgl[3];

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.
  // FIXME, calculate once => But then please recalculate whenever the dynamic window changes

  psc_foreach_patch(psc, p) {
    calc_domain_bounds(psc, p, xb, xe, xgb, xge, xgl);

    particles_t *pp = psc_mparticles_get_patch(particles, p);
    struct ddcp_patch *patch = &ddcp->patches[p];
    patch->head = 0;
    for (int dir1 = 0; dir1 < N_DIR; dir1++) {
      patch->nei[dir1].n_send = 0;
    }
    for (int i = 0; i < pp->n_part; i++) {
      particle_t *part = particles_get_one(pp, i);
      particle_real_t *xi = &part->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
      particle_real_t *pxi = &part->pxi;
      if (xi[0] >= xb[0] && xi[0] < xe[0] &&
	  xi[1] >= xb[1] && xi[1] < xe[1] &&
	  xi[2] >= xb[2] && xi[2] < xe[2]) {
	// fast path
	// inside domain: move into right position
	pp->particles[patch->head++] = *part;
      } else {
	// slow path
	int dir[3];
	for (int d = 0; d < 3; d++) {
	  if (xi[d] < xb[d]) {
	    if (xi[d] < xgb[d]) {
	      switch (psc->domain.bnd_part[d]) {
	      case BND_PART_REFLECTING:
		xi[d] = 2.f * xgb[d] - xi[d];
		pxi[d] = -pxi[d];
		dir[d] = 0;
		break;
	      case BND_PART_PERIODIC:
		xi[d] += xgl[d];
		dir[d] = -1;
		break;
	      default:
		assert(0);
	      }
	    } else {
	      // computational bnd
	      dir[d] = -1;
	    }
	  } else if (xi[d] >= xe[d]) {
	    if (xi[d] >= xge[d]) {
	      switch (psc->domain.bnd_part[d]) {
	      case BND_PART_REFLECTING:
		xi[d] = 2.f * xge[d] - xi[d];
		pxi[d] = -pxi[d];
		dir[d] = 0;
		break;
	      case BND_PART_PERIODIC:
		xi[d] -= xgl[d];
		dir[d] = +1;
		break;
	      default:
		assert(0);
	      }
	    } else {
	      dir[d] = +1;
	    }
	  } else {
	    // computational bnd
	    dir[d] = 0;
	  }
	}
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	  pp->particles[patch->head++] = *part;
	} else {
	  ddc_particles_queue(ddcp, patch, dir, part);
	}
      }
    }
  }
  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);

  psc_foreach_patch(psc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    struct ddcp_patch *patch = &ddcp->patches[p];
    pp->n_part = patch->head;
  }
  prof_stop(pr_B);

  psc_mparticles_put_cf(particles, particles_base);
}

// ----------------------------------------------------------------------
// psc_bnd_c_exchange_photons

static void
psc_bnd_c_exchange_photons(struct psc_bnd *bnd, mphotons_t *mphotons)
{
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);
  struct psc *psc = bnd->psc;
  check_domain(bnd);

  static int pr;
  if (!pr) {
    pr = prof_register("c_xchg_photon", 1., 0, 0);
  }
  prof_start(pr);

  struct ddc_particles *ddcp = bnd_c->ddcp_photons;

  double xb[3], xe[3], xgb[3], xge[3], xgl[3];

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.
  // FIXME, calculate once

  psc_foreach_patch(psc, p) {
    calc_domain_bounds(psc, p, xb, xe, xgb, xge, xgl);

    photons_t *photons = &mphotons->p[p];
    struct ddcp_patch *patch = &ddcp->patches[p];
    patch->head = 0;
    for (int dir1 = 0; dir1 < N_DIR; dir1++) {
      patch->nei[dir1].n_send = 0;
    }
    for (int i = 0; i < photons->nr; i++) {
      photon_t *ph = photons_get_one(photons, i);
      photon_real_t *xi = ph->x;
      photon_real_t *pxi = ph->p;
      if (xi[0] >= xb[0] && xi[0] <= xe[0] &&
	  xi[1] >= xb[1] && xi[1] <= xe[1] &&
	  xi[2] >= xb[2] && xi[2] <= xe[2]) {
	// fast path
	// inside domain: move into right position
	photons->photons[patch->head++] = *ph;
      } else {
	// slow path
	int dir[3];
	for (int d = 0; d < 3; d++) {
	  if (xi[d] < xb[d]) {
	    if (xi[d] < xgb[d]) {
	      switch (psc->domain.bnd_part[d]) {
	      case BND_PART_REFLECTING:
		xi[d] = 2.f * xgb[d] - xi[d];
		pxi[d] = -pxi[d];
		dir[d] = 0;
		break;
	      case BND_PART_PERIODIC:
		xi[d] += xgl[d];
		dir[d] = -1;
		break;
	      default:
		assert(0);
	      }
	    } else {
	      // computational bnd
	      dir[d] = -1;
	    }
	  } else if (xi[d] > xe[d]) {
	    if (xi[d] > xge[d]) {
	      switch (psc->domain.bnd_part[d]) {
	      case BND_PART_REFLECTING:
		xi[d] = 2.f * xge[d] - xi[d];
		pxi[d] = -pxi[d];
		dir[d] = 0;
		break;
	      case BND_PART_PERIODIC:
		xi[d] -= xgl[d];
		dir[d] = +1;
		break;
	      default:
		assert(0);
	      }
	    } else {
	      dir[d] = +1;
	    }
	  } else {
	    // computational bnd
	    dir[d] = 0;
	  }
	}
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	  photons->photons[patch->head++] = *ph;
	} else {
	  ddc_particles_queue(ddcp, patch, dir, ph);
	}
      }
    }
  }

  ddc_particles_comm(ddcp, mphotons);
  psc_foreach_patch(psc, p) {
    photons_t *photons = &mphotons->p[p];
    struct ddcp_patch *patch = &ddcp->patches[p];
    photons->nr = patch->head;
  }

  prof_stop(pr);
}


// ======================================================================
// psc_bnd: subclass "c"

struct psc_bnd_ops psc_bnd_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_bnd_c),
  .setup                 = psc_bnd_c_setup,
  .destroy               = psc_bnd_c_destroy,
  .add_ghosts            = psc_bnd_c_add_ghosts,
  .fill_ghosts           = psc_bnd_c_fill_ghosts,
  .exchange_particles    = psc_bnd_c_exchange_particles,
  .exchange_photons      = psc_bnd_c_exchange_photons,
};
