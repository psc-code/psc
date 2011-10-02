
#include "psc_bnd_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_bnd_fortran_add_ghosts

static void
psc_bnd_fortran_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			   int mb, int me)
{
  assert(ppsc->nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  mfields_fortran_t *flds = psc_mfields_fortran_get_from(mb, me, flds_base);

  for (int m = mb; m < me; m++) {
    if (ppsc->domain.gdims[0] > 1) {
      PIC_fax(&flds->f[0], m);
    }
    if (ppsc->domain.gdims[1] > 1) {
      PIC_fay(&flds->f[0], m);
    }
    if (ppsc->domain.gdims[2] > 1) {
      PIC_faz(&flds->f[0], m);
    }
  }

  psc_mfields_fortran_put_to(flds, mb, me, flds_base);

  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_fortran_fill_ghosts

static void
psc_bnd_fortran_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			    int mb, int me)
{
  assert(ppsc->nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_fill_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  mfields_fortran_t *flds = psc_mfields_fortran_get_from(mb, me, flds_base);

  for (int m = mb; m < me; m++) {
    if (ppsc->domain.gdims[0] > 1) {
      PIC_fex(&flds->f[0], m);
    }
    if (ppsc->domain.gdims[1] > 1) {
      PIC_fey(&flds->f[0], m);
    }
    if (ppsc->domain.gdims[2] > 1) {
      PIC_fez(&flds->f[0], m);
    }
  }

  psc_mfields_fortran_put_to(flds, mb, me, flds_base);

  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_fortran_exchange_particles

static void
psc_bnd_fortran_exchange_particles(struct psc_bnd *bnd,
				   mparticles_base_t *particles_base)
{
  assert(ppsc->nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_xchg_part", 1., 0, 0);
  }
  prof_start(pr);

  mparticles_fortran_t particles;
  psc_mparticles_fortran_get_from(&particles, particles_base);
  assert(ppsc->nr_patches == 1);
  particles_fortran_t *pp = &particles.p[0];

  if (ppsc->domain.gdims[0] > 1) {
    PIC_pex(pp);
  }
  if (ppsc->domain.gdims[1] > 1) {
    PIC_pey(pp);
  }
  if (ppsc->domain.gdims[2] > 1) {
    PIC_pez(pp);
  }

  psc_mparticles_fortran_put_to(&particles, particles_base);

  prof_stop(pr);
}


// ======================================================================
// psc_bnd: subclass "fortran"

struct psc_bnd_ops psc_bnd_fortran_ops = {
  .name                  = "fortran",
  .add_ghosts            = psc_bnd_fortran_add_ghosts,
  .fill_ghosts           = psc_bnd_fortran_fill_ghosts,
  .exchange_particles    = psc_bnd_fortran_exchange_particles,
};
