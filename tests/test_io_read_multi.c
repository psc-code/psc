
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrc_mod.h>
#include <mrctest.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static void
test_write_read_m1(struct mrc_fld *fld)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(fld));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/fld", "fld", fld);
  mrc_io_close(io);

  mrc_io_destroy(io);

  io = mrc_io_create(mrc_fld_comm(fld));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  mrc_io_open(io, "r", 0, 0.);
  struct mrc_fld *fld2 = mrc_io_read_path(io, "/fld", "fld", mrc_fld);
  mrc_io_close(io);

  mrc_io_destroy(io);

  mrctest_m1_compare(fld, fld2, 0.);
  mrc_fld_destroy(fld2);
}

static void
mod_domain(struct mrc_mod *mod, void *arg)
{
  struct mrctest_domain_params *par = arg;

  MPI_Comm comm = mrc_mod_get_comm(mod);
  struct mrc_domain *domain = mrctest_create_domain(comm, par);

  int testcase = 0;
  mrc_params_get_option_int("case", &testcase);

  switch (testcase) {
  case 1: ;
    struct mrc_fld *m1 = mrctest_create_m1_1(domain, 1);
    test_write_read_m1(m1);
    mrc_fld_destroy(m1);
    break;
  default:
    assert(0);
  }

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
