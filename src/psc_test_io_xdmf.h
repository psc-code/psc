
#ifndef PSC_TEST_IO_XDMF_H
#define PSC_TEST_IO_XDMF_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <mrc.h>
#include <mrc_list.h>
#include <mrc_fld.h>

struct xdmf {
  int slab_dims[3];
  int slab_off[3];
  int nr_writers;
  MPI_Comm comm_writers; //< communicator for only the writers
  int is_writer;         //< this rank is a writer
};

BEGIN_C_DECLS

void xdmf_collective_setup(struct xdmf *xdmf);
void xdmf_collective_destroy(struct xdmf *xdmf);
void xdmf_collective_write_m3(struct xdmf *xdmf, const char *path, struct mrc_fld *m3);

END_C_DECLS

#endif
