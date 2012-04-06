
#include "psc_output_photons_private.h"

// ======================================================================
// forward to subclass

void
psc_output_photons_run(struct psc_output_photons *output_photons,
			 mphotons_t *photons)
{
  struct psc_output_photons_ops *ops = psc_output_photons_ops(output_photons);
  assert(ops->run);
  psc_stats_start(st_time_output);
  ops->run(output_photons, photons);
  psc_stats_stop(st_time_output);
}

// ======================================================================
// psc_output_photons_init

static void
psc_output_photons_init()
{
  mrc_class_register_subclass(&mrc_class_psc_output_photons,
			      &psc_output_photons_none_ops);
#ifdef HAVE_HDF5
  mrc_class_register_subclass(&mrc_class_psc_output_photons,
			      &psc_output_photons_xdmf_compact_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_photons,
			      &psc_output_photons_xdmf_spread_ops);
#endif
}

// ======================================================================
// psc_output_photons class

struct mrc_class_psc_output_photons mrc_class_psc_output_photons = {
  .name             = "psc_output_photons",
  .size             = sizeof(struct psc_output_photons),
  .init             = psc_output_photons_init,
};

