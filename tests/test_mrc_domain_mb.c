
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_block_factory.h>
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrctest.h>
#include <mrc_trafo.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>

static void
test_read_write(struct mrc_domain *domain)
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

  mrctest_crds_compare(mrc_domain_get_crds(domain),
		       mrc_domain_get_crds(domain2));

  mrc_domain_destroy(domain2);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "mb");
  struct mrc_block_factory *blk_fac;
  mrc_domain_get_param_obj(domain, "block_factory", &blk_fac);

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  int testcase = 1;
  mrc_params_get_option_int("case", &testcase);

  switch (testcase) {
  case 1:
    mrc_block_factory_set_type(blk_fac, "simple3d");
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    test_read_write(domain);
    break;
  case 2: ;
    mrc_block_factory_set_type(blk_fac, "simple3d");
    mrc_crds_set_param_int(crds, "sw", 2);
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    mrc_domain_view(domain);
    mrctest_set_crds_rectilinear_1(domain);
    test_read_write(domain);
    break;
  case 3:
    mrc_block_factory_set_type(blk_fac, "cylindrical");
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    test_read_write(domain);
    break;
  case 4:
    mrc_block_factory_set_type(blk_fac, "half_cylinder");
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    test_read_write(domain);
    break;

  }
  mrc_domain_destroy(domain);

  MPI_Finalize();
}
