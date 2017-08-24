
#include "psc_generic_c.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM DIM_YZ

#include "push_part_common.c"

void
psc_push_particles_generic_c_push_mprts_yz(struct psc_push_particles *push,
					   struct psc_mparticles *mprts,
					   struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_yz", 1., 0, 0);
  }
  
  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);

    psc_fields_zero_range(flds, JXI, JXI + 3);
    do_genc_push_part(p, flds, prts);
  }
  prof_stop(pr);
}

