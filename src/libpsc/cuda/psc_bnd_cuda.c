
#include "psc_cuda.h"
#include "psc_bnd_private.h"
#include "../psc_bnd/ddc_particles.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>
#include <string.h>

// FIXME header
EXTERN_C void cuda_fill_ghosts(int p, fields_cuda_t *pf, int mb, int me);
EXTERN_C void cuda_add_ghosts(int p, fields_cuda_t *pf, int mb, int me);

struct psc_bnd_cuda {
  struct mrc_ddc *ddc;
  struct ddc_particles *ddcp;
  struct ddc_particles *ddcp_photons;
};

#define to_psc_bnd_cuda(bnd) ((struct psc_bnd_cuda *)((bnd)->obj.subctx))

// ======================================================================
// ddc funcs

static void
copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_t *flds = ctx;
  fields_t *pf = psc_mfields_get_patch(flds, p);
  fields_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf, m - mb, ix,iy,iz) = F3(pf, m, ix,iy,iz);
	}
      }
    }
  }
}

static void
add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_t *flds = ctx;
  fields_t *pf = psc_mfields_get_patch(flds, p);
  fields_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3(pf, m, ix,iy,iz) += MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

static void
copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_t *flds = ctx;
  fields_t *pf = psc_mfields_get_patch(flds, p);
  fields_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3(pf, m, ix,iy,iz) = MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

static struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = copy_to_buf,
  .copy_from_buf = copy_from_buf,
  .add_from_buf  = add_from_buf,
};

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
// psc_bnd_cuda_setup

static void
psc_bnd_cuda_setup(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);
  struct psc *psc = bnd->psc;

  bnd_cuda->ddc = mrc_domain_create_ddc(psc->mrc_domain);
  mrc_ddc_set_funcs(bnd_cuda->ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(bnd_cuda->ddc, "ibn", psc->ibn);
  mrc_ddc_set_param_int(bnd_cuda->ddc, "max_n_fields", 6);
  mrc_ddc_set_param_int(bnd_cuda->ddc, "size_of_type", sizeof(fields_real_t));
  mrc_ddc_setup(bnd_cuda->ddc);

  bnd_cuda->ddcp = ddc_particles_create(bnd_cuda->ddc, sizeof(particle_t),
				     sizeof(particle_real_t),
				     MPI_PARTICLES_REAL,
				     ddcp_particles_realloc,
				     ddcp_particles_get_addr);

  bnd_cuda->ddcp_photons = ddc_particles_create(bnd_cuda->ddc, sizeof(photon_t),
					     sizeof(photon_real_t),
					     MPI_PHOTONS_REAL,
					     ddcp_photons_realloc,
					     ddcp_photons_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_unsetup

static void
psc_bnd_cuda_unsetup(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);

  mrc_ddc_destroy(bnd_cuda->ddc);
  ddc_particles_destroy(bnd_cuda->ddcp);
  ddc_particles_destroy(bnd_cuda->ddcp_photons);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_destroy

static void
psc_bnd_cuda_destroy(struct psc_bnd *bnd)
{
  psc_bnd_cuda_unsetup(bnd);
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
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);
  struct psc *psc = bnd->psc;

  struct mrc_domain *domain = mrc_ddc_get_domain(bnd_cuda->ddc);
  if (domain != psc->mrc_domain) {
    psc_bnd_cuda_unsetup(bnd);
    psc_bnd_setup(bnd);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_add_ghosts

static void
psc_bnd_cuda_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  check_domain(bnd);

  assert(ppsc->nr_patches == 1);
  mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);
  cuda_add_ghosts(0, psc_mfields_get_patch_cuda(flds, 0), mb, me);
  prof_stop(pr);

  psc_mfields_put_cuda(flds, flds_base, mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_fill_ghosts

static void
psc_bnd_cuda_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  check_domain(bnd);

  assert(ppsc->nr_patches == 1);
  mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_fill_ghosts", 1., 0, 0);
  }
  prof_start(pr);
  cuda_fill_ghosts(0, psc_mfields_get_patch_cuda(flds, 0), mb, me);
  prof_stop(pr);

  psc_mfields_put_cuda(flds, flds_base, mb, me);
}

// ----------------------------------------------------------------------
// calc_domain_bounds
//
// calculate bounds of local patch, and global domain

#if 0
static void
calc_domain_bounds(struct psc *psc, int p, double xb[3], double xe[3],
		   double xgb[3], double xge[3], double xgl[3])
{
  struct psc_patch *psc_patch = &psc->patch[p];

  for (int d = 0; d < 3; d++) {
    xb[d] = psc_patch->off[d] * psc->dx[d];
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_PERIODIC) {
      xgb[d] = 0.;
    } else {
      if (psc->domain.bnd_fld_lo[d] == BND_FLD_UPML) {
	xgb[d] = psc->pml.size * psc->dx[d];
      } else {
	xgb[d] = 0.;
      }
      if (psc_patch->off[d] == 0) {
	xb[d] = xgb[d];
      }
    }
    
    xe[d] = (psc_patch->off[d] + psc_patch->ldims[d]) * psc->dx[d];
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_PERIODIC) {
      xge[d] = (psc->domain.gdims[d]) * psc->dx[d];
    } else {
      if (psc->domain.bnd_fld_lo[d] == BND_FLD_UPML) {
	xge[d] = (psc->domain.gdims[d]-1 - psc->pml.size) * psc->dx[d];
      } else {
	xge[d] = (psc->domain.gdims[d]-1) * psc->dx[d];
      }
      if (psc_patch->off[d] + psc_patch->ldims[d] == psc->domain.gdims[d]) {
	  xe[d] = xge[d];
      }
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
#endif

// ----------------------------------------------------------------------
// psc_bnd_cuda_exchange_particles

static void
psc_bnd_cuda_exchange_particles(struct psc_bnd *psc_bnd, mparticles_base_t *particles_base)
{
  assert(ppsc->nr_patches == 1);
  int size;
  MPI_Comm_size(psc_bnd_comm(psc_bnd), &size);
  assert(size == 1);
  // all periodic only, too
  
  mparticles_cuda_t *particles =
	  psc_mparticles_base_get_cuda(particles_base, 0);
  cuda_exchange_particles(0, psc_mparticles_get_patch_cuda(particles, 0));
  psc_mparticles_base_put_cuda(particles, particles_base);
}

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops psc_bnd_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_bnd_cuda),
  .setup                 = psc_bnd_cuda_setup,
  .destroy               = psc_bnd_cuda_destroy,
  .add_ghosts            = psc_bnd_cuda_add_ghosts,
  .fill_ghosts           = psc_bnd_cuda_fill_ghosts,
  .exchange_particles    = psc_bnd_cuda_exchange_particles,
};

