
#include "psc_bnd_private.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_bnd_fortran_add_ghosts

static void
psc_bnd_fortran_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			   int mb, int me)
{
  assert(psc.nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  mfields_fortran_t flds;
  fields_fortran_get(&flds, mb, me, flds_base);

  for (int m = mb; m < me; m++) {
    if (psc.domain.gdims[0] > 1) {
      PIC_fax(&flds.f[0], m);
    }
    if (psc.domain.gdims[1] > 1) {
      PIC_fay(&flds.f[0], m);
    }
    if (psc.domain.gdims[2] > 1) {
      PIC_faz(&flds.f[0], m);
    }
  }

  fields_fortran_put(&flds, mb, me, flds_base);

  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_fortran_fill_ghosts

static void
psc_bnd_fortran_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			    int mb, int me)
{
  assert(psc.nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_fill_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  mfields_fortran_t flds;
  fields_fortran_get(&flds, mb, me, flds_base);

  for (int m = mb; m < me; m++) {
    if (psc.domain.gdims[0] > 1) {
      PIC_fex(&flds.f[0], m);
    }
    if (psc.domain.gdims[1] > 1) {
      PIC_fey(&flds.f[0], m);
    }
    if (psc.domain.gdims[2] > 1) {
      PIC_fez(&flds.f[0], m);
    }
  }

  fields_fortran_put(&flds, mb, me, flds_base);

  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_fortran_exchange_particles

static void
psc_bnd_fortran_exchange_particles(struct psc_bnd *bnd,
				   mparticles_base_t *particles_base)
{
  assert(psc.nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_xchg_part", 1., 0, 0);
  }
  prof_start(pr);

  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  
  SET_param_coeff();
  SET_niloc(particles.p[0].n_part);

  if (psc.domain.gdims[0] > 1) {
    PIC_pex();
  }
  if (psc.domain.gdims[1] > 1) {
    PIC_pey();
  }
  if (psc.domain.gdims[2] > 1) {
    PIC_pez();
  }

  GET_niloc(&particles.p[0].n_part);
  // don't really reallocate, just get the new array pointer
  // if PIC_pe[xyz]() reallocated during the previous calls
  particles.p[0].particles = REALLOC_particles(particles.p[0].n_part);
  particles_fortran_put(&particles, particles_base);

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
