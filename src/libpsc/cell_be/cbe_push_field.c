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



static void
psc_push_fields_cbe_push_a_2d(struct psc_push_fields *push, struct psc_mfields *flds_base)
{
  struct psc_mfields flds;
  fields_get(&flds, JXI,JXI+9,flds_base);

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
    cell_run_patch(p,&flds.f[p], &null_parts, job);
  }

  wait_all_spe();

  fields_put(&flds, EX, EX + 6, flds_base);

}

static void
psc_push_fields_cbe_push_b_2d(struct psc_push_fields *push, struct psc_mfields *flds_base)
{
  struct psc_mfields flds;
  fields_get(&flds, JXI,JXI+9,flds_base);
  
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
    cell_run_patch(p,&flds.f[p], &null_parts, job);
  }

  wait_all_spe();

  fields_put(&flds, EX, EX + 6, flds_base);

  
}

// So, the latest round of refactoring has made this implementation
// inconsistent with the rest of the rest of the code. Allow me to explain 
// the problem. For the cell implementation, the most expensive thing is 
// moving data on and off the spes. If I split each step of the field push into 
// two parts (E and H), I have to move both fields onto the spu 4 each timestep, 
// and move one of them off each timestep. But it together, and for each patch 
// that's 32*32*12~=12KB moving on and 32*32*6=6KB off each timestep. The way I 
// have set up now (which is mostly equivalent, except it calculates into the ghost
// points to avoid the ghost cell exchange) only moves the fields on twice per timestep,
// and off twice per timestep. Resulting in 32*32*9=9KB on and the same 6KB off. 
// Those three KB might be insignificant, and might be cheaper than the additional 
// calculation, but I doubt it. I need to test for sure. But, for now I'm going to 
// fool the main code, and have everything done during the first push of each step, 
// and just have the second field be empty. 

// Hence the dummy function.

static void
dummy_function(struct psc_push_fields *push, struct psc_mfields *flds_base)
{
}


// ======================================================================
// psc_push_fields: subclass "cbe"


static void 
cbe_push_setup(struct psc_push_fields *push)
{
  // Initialize the spes and create the context.
  psc_init_spes();
}

static void 
cbe_push_destroy(struct psc_push_fields *push)
{
  // The spes and free the context
  psc_kill_spes();
}


struct psc_push_fields_ops psc_push_fields_cbe_ops = {
  .name                  = "cbe",
  .push_a_E              = psc_push_fields_cbe_push_a_2d,
  .push_a_H              = dummy_function,
  .push_b_H              = psc_push_fields_cbe_push_b_2d,
  .push_b_E              = dummy_function,
  .setup                 = cbe_push_setup,
  .destroy               = cbe_push_destroy,
};
