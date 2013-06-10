
#include "mrctest.h"

#include <mrc_mod.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

// ----------------------------------------------------------------------

static void
dump_field_2d(MPI_Comm comm, struct mrc_fld *fld, int rank_diagsrv)
{
  struct mrc_io *io;
  if (rank_diagsrv >= 0) {
    io = mrc_io_create(comm);
    mrc_io_set_type(io, "combined");
    mrc_io_set_param_int(io, "rank_diagsrv", rank_diagsrv);
  } else {
    io = mrc_io_create(comm);
  }
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_field2d(io, 1., fld, DIAG_TYPE_2D_Z, 99.);
  //  mrc_fld_write(fld, io);
  mrc_io_close(io);

  mrc_io_destroy(io);
}

static void
mod_domain(struct mrc_mod *mod, void *arg)
{
  struct mrctest_domain_params *par = arg;

  MPI_Comm comm = mrc_mod_get_comm(mod);
  struct mrc_domain *domain = mrctest_create_domain(comm, par);

  int mx = 8, my = 4;
  struct mrc_fld *fld = mrc_fld_create(comm);
  mrc_fld_set_param_int_array(fld, "dims", 3, (int[3]) { mx, my, 1 });
  mrc_fld_set_comp_name(fld, 0, "test_2d_0");
  mrc_fld_set_param_obj(fld, "domain", domain);
  mrc_fld_setup(fld);

  for (int iy = 0; iy < my; iy++) {
    for (int ix = 0; ix < mx; ix++) {
      MRC_F2(fld, 0, ix,iy) = 100 * iy + ix;
    }
  }

  int rank_diagsrv = mrc_mod_get_first_node(mod, "diagsrv");
  dump_field_2d(comm, fld, rank_diagsrv);
  mrc_fld_destroy(fld);

  mrc_domain_destroy(domain);
}

int
main(int argc, char **argv)
{
  mrctest_init(&argc, &argv);
  mrctest_domain(mod_domain);
  mrctest_finalize();
  return 0;
}

