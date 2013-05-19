
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrc_mod.h>
#include <mrctest.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// TODO:
// - test read/write non-uniform crds (broken)
// - test non trivial domain simple (broken)
// - reading crds by itself is broken (dont fix?)

static void
test_write_read_domain(struct mrc_domain *domain)
{
  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/domain", "domain", domain);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_domain *domain2 = mrc_io_read_path(io, "/domain", "domain", mrc_domain);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_domain_view(domain2);
  mrctest_crds_compare(mrc_domain_get_crds(domain),
		       mrc_domain_get_crds(domain2));

  mrc_domain_destroy(domain2);
}

static void
test_write_read(struct mrc_fld *fld)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(fld));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/fld", "fld", fld);
  mrc_io_close(io);

  mrc_io_destroy(io);

  io = mrc_io_create(mrc_fld_comm(fld));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);

  mrc_io_open(io, "r", 0, 0.);
  struct mrc_fld *fld2 = mrc_io_read_path(io, "/fld", "fld", mrc_fld);
  mrc_io_close(io);

  mrc_io_destroy(io);

  mrctest_fld_compare(fld, fld2, 0.);
  mrctest_crds_compare(mrc_domain_get_crds(fld->_domain),
		       mrc_domain_get_crds(fld2->_domain));

  mrc_fld_destroy(fld2);
}

static void
test_write_read_two_fields(struct mrc_domain *domain)
{
  struct mrc_fld *fld1 = mrctest_create_field_1(domain);
  mrc_fld_set_name(fld1, "fld1");
  struct mrc_fld *fld2 = mrctest_create_field_2(domain);
  mrc_fld_set_name(fld2, "fld2");

  struct mrc_io *io = mrc_io_create(mrc_fld_comm(fld1));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/fld", "fld1", fld1);
  mrc_io_write_path(io, "/fld", "fld2", fld2);
  mrc_io_close(io);

  mrc_io_destroy(io);

  io = mrc_io_create(mrc_fld_comm(fld1));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);

  mrc_io_open(io, "r", 0, 0.);
  struct mrc_fld *fld1r = mrc_io_read_path(io, "/fld", "fld1", mrc_fld);
  struct mrc_fld *fld2r = mrc_io_read_path(io, "/fld", "fld2", mrc_fld);
  mrc_io_close(io);

  mrc_io_destroy(io);

  mrctest_fld_compare(fld1, fld1r, 0.);
  mrctest_fld_compare(fld2, fld2r, 0.);

  mrc_fld_destroy(fld1);
  mrc_fld_destroy(fld2);
}

static void
test_write_read_domain_rectilinear(MPI_Comm comm, struct mrctest_domain_params *par)
{
  struct mrc_domain *domain = mrctest_create_domain_rectilinear(comm, par);
  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/domain", "domain", domain);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_domain *domain2 = mrc_io_read_path(io, "/domain", "domain", mrc_domain);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_domain_view(domain2);
  mrctest_crds_compare(mrc_domain_get_crds(domain),
		       mrc_domain_get_crds(domain2));

  mrc_domain_destroy(domain2);
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
  case 0:
    test_write_read_domain(domain);
    break;
  case 1: {
    struct mrc_fld *fld = mrctest_create_field_1(domain);
    test_write_read(fld);
    mrc_fld_destroy(fld);
    break;
  }
  case 2: {
    struct mrc_fld *fld = mrctest_create_field_2(domain);
    test_write_read(fld);
    mrc_fld_destroy(fld);
    break;
  }
  case 3:
    test_write_read_two_fields(domain);
    break;
  case 4:
    test_write_read_domain_rectilinear(comm, par);
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
