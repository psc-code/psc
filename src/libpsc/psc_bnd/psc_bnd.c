
#include "psc_bnd_private.h"
#include "ddc_particles.h"

#include <mrc_io.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>

// ======================================================================
// psc_bnd photons

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

static void
psc_bnd_photons_setup(struct psc_bnd *bnd)
{
  bnd->ddcp_photons = ddc_particles_create(bnd->ddc, sizeof(photon_t),
					   sizeof(photon_real_t),
					   MPI_PHOTONS_REAL,
					   ddcp_photons_realloc,
					   ddcp_photons_get_addr);
}

static void
psc_bnd_photons_unsetup(struct psc_bnd *bnd)
{
  ddc_particles_destroy(bnd->ddcp_photons);
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
// psc_bnd_photons_exchange

static void
psc_bnd_photons_exchange(struct psc_bnd *bnd, mphotons_t *mphotons)
{
  struct psc *psc = bnd->psc;

  static int pr;
  if (!pr) {
    pr = prof_register("xchg_photon", 1., 0, 0);
  }
  prof_start(pr);

  struct ddc_particles *ddcp = bnd->ddcp_photons;

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
	      switch (psc->domain.bnd_part_lo[d]) {
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
	      switch (psc->domain.bnd_part_hi[d]) {
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
// psc_bnd

// ----------------------------------------------------------------------
// psc_set_psc

void
psc_bnd_set_psc(struct psc_bnd *bnd, struct psc *psc)
{
  bnd->psc = psc;
}

// ----------------------------------------------------------------------
// psc_bnd_setup

static void
_psc_bnd_setup(struct psc_bnd *bnd)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->create_ddc);
  ops->create_ddc(bnd);
  psc_bnd_photons_setup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_destroy

static void
_psc_bnd_destroy(struct psc_bnd *bnd)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->unsetup);
  ops->unsetup(bnd);
  psc_bnd_photons_unsetup(bnd);
  mrc_ddc_destroy(bnd->ddc);
}

// ----------------------------------------------------------------------
// psc_bnd_write

static void
_psc_bnd_write(struct psc_bnd *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_name(bnd);
  mrc_io_write_obj_ref(io, path, "psc", (struct mrc_obj *) bnd->psc);
}

// ----------------------------------------------------------------------
// psc_bnd_read

static void
_psc_bnd_read(struct psc_bnd *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_name(bnd);
  bnd->psc = (struct psc *)
    mrc_io_read_obj_ref(io, path, "psc", &mrc_class_psc);

  psc_bnd_setup(bnd);
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
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);

  struct mrc_domain *domain = mrc_ddc_get_domain(bnd->ddc);
  if (domain != bnd->psc->mrc_domain) {
    ops->unsetup(bnd);
    psc_bnd_photons_unsetup(bnd);
    mrc_ddc_destroy(bnd->ddc);

    ops->create_ddc(bnd);
    ops->setup(bnd);
    psc_bnd_photons_setup(bnd);
  }
}

// ======================================================================
// forward to subclass

void
psc_bnd_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  check_domain(bnd);

  psc_stats_start(st_time_comm);
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->add_ghosts);
  ops->add_ghosts(bnd, flds, mb, me);
  psc_stats_stop(st_time_comm);
}

void
psc_bnd_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  check_domain(bnd);

  psc_stats_start(st_time_comm);
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->fill_ghosts);
  ops->fill_ghosts(bnd, flds, mb, me);
  psc_stats_stop(st_time_comm);
}

void
psc_bnd_exchange_particles(struct psc_bnd *bnd, mparticles_base_t *particles)
{
  static int pr;
  if (!pr) {
    pr = prof_register("xchg_prts", 1., 0, 0);
  }

  check_domain(bnd);

  prof_start(pr);
  psc_stats_start(st_time_comm);
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->exchange_particles);
  ops->exchange_particles(bnd, particles);
  psc_stats_stop(st_time_comm);
  prof_stop(pr);
}

void
psc_bnd_exchange_photons(struct psc_bnd *bnd, mphotons_t *mphotons)
{
  check_domain(bnd);

  psc_stats_start(st_time_comm);
  int n_total = 0;
  psc_foreach_patch(bnd->psc, p) {
    n_total += mphotons->p[p].nr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &n_total, 1, MPI_INT, MPI_SUM, 
		psc_bnd_comm(bnd));
  if (n_total > 0) {
    psc_bnd_photons_exchange(bnd, mphotons);
  }
  psc_stats_stop(st_time_comm);
}

// ======================================================================
// psc_bnd_init

static void
psc_bnd_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_auto_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_single2_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_mix_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_cuda_ops);
#endif
}

// ======================================================================
// psc_bnd class

struct mrc_class_psc_bnd mrc_class_psc_bnd = {
  .name             = "psc_bnd",
  .size             = sizeof(struct psc_bnd),
  .init             = psc_bnd_init,
  .setup            = _psc_bnd_setup,
  .destroy          = _psc_bnd_destroy,
  .write            = _psc_bnd_write,
  .read             = _psc_bnd_read,
};

