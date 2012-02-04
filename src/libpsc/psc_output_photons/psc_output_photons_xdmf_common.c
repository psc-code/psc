/*
 *  XDMF/HDF5 photon output by Simon Jagoda
 *	Intended to be included in case, this is the place for setter functions
 *	concerning options available in both xdmf ops
 */
 
#include "psc_output_photons_xdmf.h"

#define to_psc_output_photons_xdmf(out) ((struct psc_output_photons_xdmf *)((out)->obj.subctx))



void photons_xdmf_set_output_step (struct psc_output_photons *out, int step, bool choice) {
	struct psc_output_photons_xdmf *out_c = to_psc_output_photons_xdmf(out);
	out_c->write_photons[step]=choice;
}
