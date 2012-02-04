
#ifndef PSC_OUTPUT_PHOTONS_PRIVATE_H
#define PSC_OUTPUT_PHOTONS_PRIVATE_H

#include <psc_output_photons.h>

struct psc_output_photons {
  struct mrc_obj obj;
};

struct psc_output_photons_ops {
  MRC_SUBCLASS_OPS(struct psc_output_photons);
  void (*run)(struct psc_output_photons *output_photons,
	      mphotons_t *photons);
};

// ======================================================================

extern struct psc_output_photons_ops psc_output_photons_none_ops;
extern struct psc_output_photons_ops psc_output_photons_c_ops;
extern struct psc_output_photons_ops psc_output_photons_xdmf_compact_ops;
extern struct psc_output_photons_ops psc_output_photons_xdmf_spread_ops;

#define psc_output_photons_ops(output_photons) ((struct psc_output_photons_ops *)((output_photons)->obj.ops))

#endif
