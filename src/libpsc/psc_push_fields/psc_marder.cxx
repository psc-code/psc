
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"

#include <mrc_io.h>
#include <mrc_profile.h>

#include <math.h>

// ----------------------------------------------------------------------
// marder_calc_aid_fields

static void
marder_calc_aid_fields(struct psc_marder *marder, 
		       struct psc_mfields *flds, struct psc_mparticles *particles,
		       struct psc_mfields *div_e, struct psc_mfields *rho)
{
  psc_output_fields_item_run(marder->item_div_e, flds, particles, div_e); // FIXME, should accept NULL for particles
  
  if (marder->dump) {
    static int cnt;
    mrc_io_open(marder->io, "w", cnt, cnt);//ppsc->timestep, ppsc->timestep * ppsc->dt);
    cnt++;
    psc_mfields_write_as_mrc_fld(rho, marder->io);
    psc_mfields_write_as_mrc_fld(div_e, marder->io);
    mrc_io_close(marder->io);
  }

  div_e->axpy(0, -1., rho, 0);
  // FIXME, why is this necessary?
  psc_bnd_fill_ghosts(marder->bnd, div_e, 0, 1);
}

// ----------------------------------------------------------------------
// psc_marder_run
//
// On ghost cells:
// It is possible (variant = 1) that ghost cells are set before this is called
// and the subsequent code expects ghost cells to still be set on return.
// We're calling fill_ghosts at the end of each iteration, so that's fine.
// However, for variant = 0, ghost cells aren't set on entry, and they're not
// expected to be set on return (though we do that, anyway.)

void
psc_marder_run(struct psc_marder *marder, 
	       struct psc_mfields *flds, struct psc_mparticles *particles)
{
  static int pr, pr_A, pr_B;
  if (!pr) {
    pr   = prof_register("psc_marder_run", 1., 0, 0);
    pr_A = prof_register("psc_marder rho", 1., 0, 0);
    pr_B = prof_register("psc_marder iter", 1., 0, 0);
  }

  if (marder->every_step < 0 || ppsc->timestep % marder->every_step != 0) 
    return;

  prof_start(pr);
  struct psc_marder_ops *ops = psc_marder_ops(marder);
  assert(ops);

  if (ops->run) {
    ops->run(marder, flds, particles);
  } else {
    prof_start(pr_A);

    psc_output_fields_item_run(marder->item_rho, flds, particles, marder->rho);

    // need to fill ghost cells first (should be unnecessary with only variant 1) FIXME
    psc_bnd_fill_ghosts(ppsc->bnd, flds, EX, EX+3);

    prof_stop(pr_A);

    prof_start(pr_B);
    for (int i = 0; i < marder->loop; i++) {
      marder_calc_aid_fields(marder, flds, particles, marder->div_e, marder->rho);
      ops->correct(marder, flds, marder->div_e);
      psc_bnd_fill_ghosts(ppsc->bnd, flds, EX, EX+3);
    }
    prof_stop(pr_B);
  }
  
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_marder_init

static void
psc_marder_init()
{
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_single_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_cuda_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_vpic_ops);
#endif
}

// ----------------------------------------------------------------------
// psc_marder description

#define VAR(x) (void *)offsetof(struct psc_marder, x)
static struct param psc_marder_descr[] = {
  { "every_step"       , VAR(every_step)       , PARAM_INT(-1),     },
  { "diffusion"        , VAR(diffusion)        , PARAM_DOUBLE(0.9), },
  { "loop"             , VAR(loop)             , PARAM_INT(1),      },
  { "dump"             , VAR(dump)             , PARAM_BOOL(false), },

  { "clean_div_e_interval", VAR(clean_div_e_interval), PARAM_INT(0),     },
  { "clean_div_b_interval", VAR(clean_div_b_interval), PARAM_INT(0),     },
  { "sync_shared_interval", VAR(sync_shared_interval), PARAM_INT(0),     },
  { "num_div_e_round"     , VAR(num_div_e_round)     , PARAM_INT(0),     },
  { "num_div_b_round"     , VAR(num_div_b_round)     , PARAM_INT(0),     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_marder class description

struct mrc_class_psc_marder mrc_class_psc_marder = {
  .name             = "psc_marder",
  .size             = sizeof(struct psc_marder),
  .param_descr      = psc_marder_descr,
  .init             = psc_marder_init,
};

