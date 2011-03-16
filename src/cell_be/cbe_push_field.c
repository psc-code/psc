#include <string.h>
#include "psc_ppu.h"
#include <mrc_profile.h>

#include "psc_push_fields_private.h"
#include "psc_bnd.h"

#ifdef CELLEMU
#include "libspe2_c.h"
#else
#include <libspe2.h>
#endif

void
cbe_push_field_a_2d(mfields_base_t *flds_base)
{
  mfields_t flds;
  fields_get(&flds, JXI,JXI+9,flds_base);

  if(!spes_inited)
    psc_init_spes();
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_field_2d_a", 1., 0, 0);
  }
  prof_start(pr);

  int job = SPU_FIELD_A;

  particles_t null_parts; 
  null_parts.particles = NULL;
  null_parts.n_part = 0;


  // I'm not really sure about this. I mean, it's the nicest way to do this 
  // in terms of reusing code, but it also means I'm basically dooming the ppu to not
  // do anything but manage the spes. I think we could do more with it, but that's
  // a task for another day...
  foreach_patch(p) {
    // So, another thing I'm not too please about. The function that's getting
    // called here will basically stall the ppu until it manages to start the patch
    // on an spe...
    cell_run_patch(&flds.f[p], &null_parts, job);
  }

  wait_all_spe();

  prof_stop(pr);

  fields_put(&flds, EX, EX + 6, flds_base);

  psc_bnd_fill_ghosts(psc.bnd,flds_base, EX, EX + 6);  
}

void
cbe_push_field_b_2d(mfields_base_t *flds_base)
{
  mfields_t flds;
  fields_get(&flds, JXI,JXI+9,flds_base);

  if(!spes_inited)
    psc_init_spes();
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_field_2d_b", 1., 0, 0);
  }
  prof_start(pr);

  int job = SPU_FIELD_B;

  particles_t null_parts; 
  null_parts.particles = NULL;
  null_parts.n_part = 0;


  // I'm not really sure about this. I mean, it's the nicest way to do this 
  // in terms of reusing code, but it also means I'm basically dooming the ppu to not
  // do anything but manage the spes. I think we could do more with it, but that's
  // a task for another day...
  foreach_patch(p) {
    // So, another thing I'm not too please about. The function that's getting
    // called here will basically stall the ppu until it manages to start the patch
    // on an spe...
    cell_run_patch(&flds.f[p], &null_parts, job);
  }

  wait_all_spe();

  prof_stop(pr);

  fields_put(&flds, EX, EX + 6, flds_base);

  psc_bnd_fill_ghosts(psc.bnd,flds_base, EX, EX + 6);
  
}

static void
psc_push_fields_cbe_step_a(struct psc_push_fields *push, mfields_base_t *flds)
{
  if (psc.domain.use_pml) {
    assert(0);
  } else {
    cbe_push_field_a_2d(flds);
  }
}

static void
psc_push_fields_cbe_step_b(struct psc_push_fields *push, mfields_base_t *flds)
{
  if (psc.domain.use_pml) {
    assert(0);
  } else {
    cbe_push_field_b_2d(flds);
  }
}


// ======================================================================
// psc_push_fields: subclass "cbe"

struct psc_push_fields_ops psc_push_fields_cbe_ops = {
  .name                  = "cbe",
  .step_a                = psc_push_fields_cbe_step_a,
  .step_b                = psc_push_fields_cbe_step_b,
};
