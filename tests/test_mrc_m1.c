
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>

static void
set_m1(struct mrc_fld *m1)
{
  struct mrc_patch *patches = mrc_domain_get_patches(m1->_domain, NULL);
  struct mrc_crds *crds = mrc_domain_get_crds(m1->_domain);

  mrc_m1_foreach_patch(m1, p) {
    int *off = patches[p].off;
    mrc_m1_foreach(m1, ix, 0,0) {
      MRC_M1(m1, 0, ix, p) = ix + off[0];
      MRC_M1(m1, 1, ix, p) = MRC_MCRD(crds, 0, ix, p);
    } mrc_m1_foreach_end;
  }
}

static void
check_m1(struct mrc_fld *m1)
{
  struct mrc_patch *patches = mrc_domain_get_patches(m1->_domain, NULL);

  mrc_m1_foreach_patch(m1, p) {
    int *off = patches[p].off;
    mrc_m1_foreach(m1, ix, 0,0) {
      assert(MRC_M1(m1, 0, ix, p) == ix + off[0]);
    } mrc_m1_foreach_end;
  }
}

static void
test_write_m1(struct mrc_fld *m1)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(m1));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/m1", "m1", m1);
  mrc_io_close(io);

  mrc_io_open(io, "w", 1, 1.);
  mrc_io_write_path(io, "/m1", "m1", m1);
  mrc_io_close(io);

  mrc_io_destroy(io);
}

static void
test_write_read_m1(struct mrc_fld *m1)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(m1));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/m1", "m1", m1);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_fld_comm(m1));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_fld *m1_2 = mrc_io_read_path(io, "/m1", "m1", mrc_fld);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrctest_m1_compare(m1, m1_2, 1e-7);
  mrc_fld_destroy(m1_2);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int testcase = 1;
  mrc_params_get_option_int("case", &testcase);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "multi");
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_param_int(crds, "sw", 1);
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);
  
  if (strcmp(mrc_crds_type(crds), "rectilinear") == 0) {
    mrctest_set_crds_rectilinear_1(domain);
  }

  struct mrc_fld *m1 = mrc_domain_m1_create(domain);
  mrc_fld_set_name(m1, "test_m1");
  mrc_fld_set_param_int(m1, "dim", 0);
  mrc_fld_set_param_int(m1, "nr_comps", 2);
  mrc_fld_set_from_options(m1);
  mrc_fld_setup(m1);
  mrc_fld_set_comp_name(m1, 0, "fld0");
  mrc_fld_set_comp_name(m1, 1, "fld1");
  mrc_fld_view(m1);
  
  set_m1(m1);
  check_m1(m1);
  switch (testcase) {
  case 1:
    test_write_m1(m1);
    break;
  case 2:
    test_write_read_m1(m1);
    break;
  }

  mrc_fld_destroy(m1);
  mrc_domain_destroy(domain);

  MPI_Finalize();
}
