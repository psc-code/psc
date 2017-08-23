#include <string.h>
#include "psc_ppu.h"
#include <mrc_profile.h>

#ifdef CELLEMU
#include "libspe2_c.h"
#else
#include <libspe2.h>
#endif

void
psc_push_particles_cbe_push_xy(struct psc_push_particles *push, 
			       struct psc_mparticles *particles_base,
			       struct psc_mfields *flds_base)
{
  
  struct psc_mfields flds;
  struct psc_mparticles particles;
  fields_get(&flds, EX, EX +6,flds_base);
  mparticles_get(&particles, particles_base);
  
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_part_2d", 1., 0, 0);
  }
  prof_start(pr);

  int job = SPU_PART;

  psc_mfields_zero(&flds, JXI);
  psc_mfields_zero(&flds, JYI);
  psc_mfields_zero(&flds, JZI);
  // I'm not really sure about this. I mean, it's the nicest way to do this 
  // in terms of reusing code, but it also means I'm basically dooming the ppu to not
  // do anything but manage the spes. I think we could do more with it, but that's
  // a task for another day...
  foreach_patch(p) {
    // So, another thing I'm not too please about. The function that's getting
    // called here will basically stall the ppu until it manages to start the patch
    // on an spe...
    cell_run_patch(p,&flds.f[p], &particles.p[p], job);
  }

  wait_all_spe();

  prof_stop(pr);

  fields_put(&flds, JXI, JXI + 3, flds_base);
  mparticles_put(&particles, particles_base);
}
