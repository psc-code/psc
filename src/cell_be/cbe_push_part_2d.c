#include <string.h>
#include "psc_ppu.h"
#include <mrc_profile.h>

#ifdef CELLEMU
#include "libspe2_c.h"
#else
#include <libspe2.h>
#endif

void
cbe_push_part_2d(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mfields_t flds;
  mparticles_t particles;
  fields_get(&flds, EX, EX +6,flds_base);
  particles_get(&particles, particles_base);

  particles_cbe_get(&pp);
  fields_c_get(&pf, EX, EX+6);

  if(spes_inited)
    psc_init_spes();

  // This may not be needed anymore...
  cbe_field_blocks_get(&pf,EX,EX+6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_part_yz", 1., 0, psc.pp.n_part * 9 * sizeof(cbe_real));
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
    cell_run_patch(&flds.f[p], &particles.p[p], job);
  }

  wait_all_spe(void);

  prof_stop(pr);

  fields_put(&pf,JXI,JXI+3);
  particles_put(&pp);
}
