
#include "psc_push_particles_private.h"
#include "psc_glue.h"

#include <mrc_profile.h>


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
};
