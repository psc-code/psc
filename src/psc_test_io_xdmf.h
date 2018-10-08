
#ifndef PSC_TEST_IO_XDMF_H
#define PSC_TEST_IO_XDMF_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <mrc.h>
#include <mrc_list.h>
#include <mrc_fld.h>

struct xdmf_file {
  hid_t h5_file;
  list_t xdmf_spatial_list;
};

struct xdmf {
  int slab_dims[3];
  int slab_off[3];
  struct xdmf_file file;
  struct xdmf_temporal *xdmf_temporal;
  bool use_independent_io;
  char *romio_cb_write;
  char *romio_ds_write;
  int nr_writers;
  MPI_Comm comm_writers; //< communicator for only the writers
  int *writers;          //< rank (in mrc_io comm) for each writer
  int is_writer;         //< this rank is a writer
  char *mode;            //< open mode, "r" or "w"
};

BEGIN_C_DECLS

void xdmf_collective_setup(struct xdmf *xdmf);
void xdmf_collective_destroy(struct xdmf *xdmf);
void xdmf_collective_write_m3(struct xdmf *xdmf, const char *path, struct mrc_fld *m3);

END_C_DECLS

#endif
