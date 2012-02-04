/*
 *  XDMF/HDF5 particle output by Simon Jagoda
 *	Intended to be included in case, this is the place for setter functions
 *	concerning options available in both xdmf ops
 */
 
#include "psc_output_particles_xdmf.h"

#define to_psc_output_particles_xdmf(out) ((struct psc_output_particles_xdmf *)((out)->obj.subctx))



void particles_xdmf_set_output_step (struct psc_output_particles *out, int step, bool choice) {
	struct psc_output_particles_xdmf *out_c = to_psc_output_particles_xdmf(out);
	out_c->write_particles[step]=choice;
}
