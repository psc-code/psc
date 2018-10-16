
#ifndef PSC_TEST_IO_XDMF_H
#define PSC_TEST_IO_XDMF_H

#include <mrc.h>

#include <stdio.h>

BEGIN_C_DECLS

// ======================================================================
// mock_domain

struct mock_patch {
  int off[3];
  int ldims[3];
  int rank;
};

struct xdmf {
  int gdims[3];

  int nr_global_patches;
  struct mock_patch *global_patches;

  int nr_patches;
  struct mock_patch *patches;

  int slow_dim;
  int slow_indices_per_writer;
  int slow_indices_rmndr;

  int nr_writers;
  MPI_Comm comm_writers; //< communicator for only the writers
  int is_writer;         //< this rank is a writer
};

void xdmf_collective_setup(struct xdmf *xdmf, int nr_writers, int gdims[3], int np[3]);
void xdmf_collective_destroy(struct xdmf *xdmf);
void xdmf_collective_write_m3(struct xdmf* xdmf);


END_C_DECLS

#endif
