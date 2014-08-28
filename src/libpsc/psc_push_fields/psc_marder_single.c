
#include "psc_fields_as_single.h"
#include "psc_particles_as_single.h"

#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"

#include <mrc_io.h>

// FIXME: checkpointing won't properly restore state
// FIXME: if they subclass creates objects, it'd be cleaner to have them
// be part of the subclass

// ----------------------------------------------------------------------
// fld_create
//
// FIXME, should be consolidated with psc_checks.c, and probably other places

static struct psc_mfields *
fld_create(struct psc *psc, const char *name)
{
  struct psc_mfields *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(fld, FIELDS_TYPE);
  psc_mfields_set_domain(fld, psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_set_param_int(fld, "nr_fields", 1);
  psc_mfields_setup(fld);
  psc_mfields_set_comp_name(fld, 0, name);

  return fld;
}

// ----------------------------------------------------------------------
// psc_marder_sub_setup

static void
psc_marder_sub_setup(struct psc_marder *marder)
{
  marder->div_e = fld_create(ppsc, "div_E");
  marder->rho = fld_create(ppsc, "rho");

  marder->bnd = psc_bnd_create(psc_marder_comm(marder));
  psc_bnd_set_name(marder->bnd, "marder_bnd");
  psc_bnd_set_type(marder->bnd, FIELDS_TYPE);
  psc_bnd_set_psc(marder->bnd, ppsc);
  psc_bnd_setup(marder->bnd);

  // FIXME, output_fields should be taking care of their own psc_bnd?
  marder->item_div_e = psc_output_fields_item_create(psc_comm(ppsc));
  psc_output_fields_item_set_type(marder->item_div_e, "dive_" FIELDS_TYPE);
  psc_output_fields_item_set_psc_bnd(marder->item_div_e, marder->bnd);
  psc_output_fields_item_setup(marder->item_div_e);

  marder->item_rho = psc_output_fields_item_create(psc_comm(ppsc));
  psc_output_fields_item_set_type(marder->item_rho, "rho_1st_nc_" PARTICLE_TYPE);
  psc_output_fields_item_set_psc_bnd(marder->item_rho, marder->bnd);
  psc_output_fields_item_setup(marder->item_rho);

  if (marder->dump) {
    struct mrc_io *io = mrc_io_create(psc_comm(ppsc));
    mrc_io_set_type(io, "xdmf_collective");
    mrc_io_set_name(io, "mrc_io_marder");
    mrc_io_set_param_string(io, "basename", "marder");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);

    marder->io = io;
  }
}

// ----------------------------------------------------------------------
// psc_marder_sub_destroy

static void
psc_marder_sub_destroy(struct psc_marder *marder)
{
  psc_mfields_destroy(marder->div_e);
  psc_mfields_destroy(marder->rho);

  psc_output_fields_item_destroy(marder->item_div_e);
  psc_output_fields_item_destroy(marder->item_rho);

  psc_bnd_destroy(marder->bnd);

  if (marder->dump) {
    mrc_io_destroy(marder->io);
  }
}

// ======================================================================
// psc_marder: subclass "single"

struct psc_marder_ops psc_marder_single_ops = {
  .name                  = FIELDS_TYPE,
  .setup                 = psc_marder_sub_setup,
  .destroy               = psc_marder_sub_destroy,
};

