
#include "psc_push_particles_private.h"
#include "psc_glue.h"

#include <mrc_profile.h>

#if 0 // FIXME, convert -> vay

// ----------------------------------------------------------------------
// psc_push_particles_push_xy

static void
psc_push_particles_vay_push_xy(struct psc_push_particles *push,
				   mparticles_base_t *particles_base,
				   mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_xy", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    PIC_push_part_xy(ppsc, p, &particles.p[p], &flds.f[p]);
  }
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

// ----------------------------------------------------------------------
// psc_push_particles_push_xz

static void
psc_push_particles_vay_push_xz(struct psc_push_particles *push,
				   mparticles_base_t *particles_base,
				   mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_xz", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    PIC_push_part_xz(ppsc, p, &particles.p[p], &flds.f[p]);
  }
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

// ----------------------------------------------------------------------
// psc_push_particles_push_yz

static void
psc_push_particles_vay_push_yz(struct psc_push_particles *push,
				   mparticles_base_t *particles_base,
				   mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_yz", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    PIC_push_part_yz(ppsc, p, &particles.p[p], &flds.f[p]);
  }
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

// ----------------------------------------------------------------------
// psc_push_particles_push_xyz

static void
psc_push_particles_vay_push_xyz(struct psc_push_particles *push,
				   mparticles_base_t *particles_base,
				   mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_xyz", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    PIC_push_part_xyz(ppsc, p, &particles.p[p], &flds.f[p]);
  }
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

#endif

// ----------------------------------------------------------------------
// psc_push_particles_vay_push_a_z

static void
psc_push_particles_vay_push_a_z(struct psc_push_particles *push,
				struct psc_particles *prts_base,
				struct psc_fields *flds_base)
{
  assert(ppsc->nr_patches == 1);
  
  static int pr;
  if (!pr) {
    pr = prof_register("vay_part_z", 1., 0, 0);
  }

  struct psc_particles *prts = psc_particles_get_as(prts_base, "fortran", 0);
  struct psc_fields *flds = psc_fields_get_as(flds_base, "fortran", EX, EX + 6);
  
  prof_start(pr);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  PIC_push_part_z_vay(ppsc, prts->p, prts, flds);
  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, 0);
  psc_fields_put_as(flds, flds_base, JXI, JXI + 3);
}

// ======================================================================
// psc_push_particles: subclass "vay"

struct psc_push_particles_ops psc_push_particles_vay_ops = {
  .name                  = "vay",
  .push_a_z              = psc_push_particles_vay_push_a_z,
#if 0 // FIXME, convert to vay
  .push_xy               = psc_push_particles_vay_push_xy,
  .push_xz               = psc_push_particles_vay_push_xz,
  .push_yz               = psc_push_particles_vay_push_yz,
  .push_xyz              = psc_push_particles_vay_push_xyz,
#endif
};
