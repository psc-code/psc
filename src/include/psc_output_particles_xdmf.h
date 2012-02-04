/*
 *  XDMF/HDF5 particle output by Simon Jagoda
 *	This struct is used by both xdmf output ops
 */
 
#ifndef PSC_OUTPUT_PARTICLES_XDMF_H
#define PSC_OUTPUT_PARTICLES_XDMF_H

#include "psc_output_particles_private.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_c.h"

#include "hdf5.h"

#include <stdbool.h>


struct psc_output_particles_xdmf {
  char *data_dir;	//directory for file output
  char *basename;
  struct mrc_io *io;	//no clue what this is or does, but I suppose the rest of the PSC goes on a rampage if it's missing
  int mpi_size;
  bool (*filter_func)(particle_t *part);
  int first;
  int step;
  bool *write_particles; //huge, yet lightweight array of size nmax
  hid_t hdf5;		//file handlers
  FILE *xdmf;
  char hdf5filename[80]; //set when file is created: the pointer to the file is passed around, but xdmf needs to know the actual name
  char *typenames[9];//"x,y,z,px,py,pz,q,m,w"
  bool write_variable[9];
  int times_written;

};

void particles_xdmf_set_output_step(struct psc_output_particles *out, int step, bool choice);


#endif
