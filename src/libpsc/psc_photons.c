
#include "psc.h"
#include "psc_particles_c.h"

#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>

// ----------------------------------------------------------------------
// photons_alloc

void
photons_alloc(photons_t *pp, int nr)
{
  pp->nr_alloced = nr * 1.2;
  pp->photons = calloc(pp->nr_alloced, sizeof(*pp->photons));
}

// ----------------------------------------------------------------------
// photons_realloc

void
photons_realloc(photons_t *pp, int new_nr)
{
  if (new_nr <= pp->nr)
    return;

  pp->nr_alloced = new_nr * 1.5;
  pp->photons = realloc(pp->photons, pp->nr_alloced * sizeof(*pp->photons));
}

// ----------------------------------------------------------------------
// photons_free

void
photons_free(photons_t *pp)
{
  free(pp->photons);
  pp->nr_alloced = 0;
  pp->photons = NULL;
}

// ----------------------------------------------------------------------
// psc_mphotons_set_domain

void
psc_mphotons_set_domain(mphotons_t *mphotons, struct mrc_domain *domain)
{
  mphotons->domain = domain;
}

// ----------------------------------------------------------------------
// psc_mphotons_setup

static void
_psc_mphotons_setup(mphotons_t *mphotons)
{
  mrc_domain_get_patches(mphotons->domain, &mphotons->nr_patches);
  mphotons->p = calloc(mphotons->nr_patches, sizeof(*mphotons->p));
}

// ----------------------------------------------------------------------
// psc_mphotons_destroy

static void
_psc_mphotons_destroy(mphotons_t *mphotons)
{
  for (int p = 0; p < mphotons->nr_patches; p++) { 
    photons_free(&mphotons->p[p]);
  }
  free(mphotons->p);
  mphotons->nr_patches = -1;
}

static void
_psc_mphotons_write(mphotons_t *mphotons, struct mrc_io *io)
{
  mrc_io_write_ref(io, mphotons, "domain", mphotons->domain);
}

static void
_psc_mphotons_read(mphotons_t *mphotons, struct mrc_io *io)
{
  mphotons->domain = mrc_io_read_ref(io, mphotons, "domain", mrc_domain);

  // FIXME, need to actually write / read the data, too
  psc_mphotons_setup(mphotons);
}

// ======================================================================
// psc_mphotons

struct mrc_class_psc_mphotons mrc_class_psc_mphotons = {
  .name             = "psc_mphotons",
  .size             = sizeof(struct psc_mphotons),
  .setup            = _psc_mphotons_setup,
  .destroy          = _psc_mphotons_destroy,
  .write            = _psc_mphotons_write,
  .read             = _psc_mphotons_read,
};

