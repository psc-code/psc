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
			       mparticles_base_t *particles_base,
			       mfields_base_t *flds_base)
{
  
  mfields_t flds;
  mparticles_t particles;
  fields_get(&flds, EX, EX +6,flds_base);
  particles_get(&particles, particles_base);
  
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_part_2d", 1., 0, 0);
  }
  prof_start(pr);

  int job = SPU_PART;

  // I'm not really sure about this. I mean, it's the nicest way to do this 
  // in terms of reusing code, but it also means I'm basically dooming the ppu to not
  // do anything but manage the spes. I think we could do more with it, but that's
  // a task for another day...
  foreach_patch(p) {
    // So, another thing I'm not too please about. The function that's getting
    // called here will basically stall the ppu until it manages to start the patch
    // on an spe...
    fields_zero(&flds.f[p], JXI);
    fields_zero(&flds.f[p], JYI);
    fields_zero(&flds.f[p], JZI);
    cell_run_patch(p,&flds.f[p], &particles.p[p], job);
  }

  wait_all_spe();

  prof_stop(pr);

  fields_put(&flds, JXI, JXI + 3, flds_base);
  particles_put(&particles, particles_base);
}
