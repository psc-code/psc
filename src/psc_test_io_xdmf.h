
#ifndef PSC_TEST_IO_XDMF_H
#define PSC_TEST_IO_XDMF_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <mrc.h>
#include <mrc_list.h>
#include <mrc_fld.h>

BEGIN_C_DECLS

// ======================================================================
// mock_domain

struct mock_domain {
  struct mrc_domain *domain;
  int gdims[3];

  int nr_global_patches;
  struct mrc_patch_info *patch_info;

  int nr_patches;
  struct mrc_patch *patches;
};

void mock_domain_init(struct mock_domain *mock, struct mrc_domain *domain);

struct xdmf {
  int slab_dims[3];
  int slab_off[3];
  int nr_writers;
  MPI_Comm comm_writers; //< communicator for only the writers
  int is_writer;         //< this rank is a writer
};

void xdmf_collective_setup(struct xdmf *xdmf);
void xdmf_collective_destroy(struct xdmf *xdmf);
void xdmf_collective_write_m3(struct xdmf* xdmf, struct mock_domain *mock);

END_C_DECLS

#endif
