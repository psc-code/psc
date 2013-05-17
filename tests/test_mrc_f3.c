
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>

static void
set_f3(struct mrc_f3 *f3)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f3->_domain);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(f3->_domain, 0, &info);
  int *off = info.off;
  mrc_f3_foreach(f3, ix,iy,iz, 0,0) {
    MRC_F3(f3, 0, ix,iy,iz) =
      (iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]);
    MRC_F3(f3, 1, ix,iy,iz) = MRC_CRD(crds, 0, ix);
  } mrc_f3_foreach_end;
}

static void
check_f3(struct mrc_f3 *f3)
{
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(f3->_domain, 0, &info);
  int *off = info.off;
  mrc_f3_foreach(f3, ix,iy,iz, 0,0) {
    assert(MRC_F3(f3, 0, ix,iy,iz) ==
	   (iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]));
  } mrc_f3_foreach_end;
}

static void
test_write_f3(struct mrc_f3 *f3)
{
  struct mrc_io *io = mrc_io_create(mrc_f3_comm(f3));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/f3", "f3", f3);
  mrc_io_close(io);

  mrc_io_open(io, "w", 1, 1.);
  mrc_io_write_path(io, "/f3", "f3", f3);
  mrc_io_close(io);

  mrc_io_destroy(io);
}

static void
test_write_read_f3(struct mrc_f3 *f3)
{
  struct mrc_io *io = mrc_io_create(mrc_f3_comm(f3));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/f3", "f3", f3);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_f3_comm(f3));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_f3 *f3_2 = mrc_io_read_path(io, "/f3", "f3", mrc_f3);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrctest_f3_compare(f3, f3_2, 0.);
  mrc_f3_destroy(f3_2);
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
  mrc_crds_set_type(crds, "uniform");
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);
  
  if (strcmp(mrc_crds_type(crds), "rectilinear") == 0) {
    mrctest_set_crds_rectilinear_1(domain);
  }

  struct mrc_f3 *f3 = mrc_domain_f3_create(domain, 2);
  mrc_f3_set_name(f3, "test_f3");
  mrc_f3_set_param_int(f3, "nr_comps", 2);
  mrc_f3_set_from_options(f3);
  mrc_f3_setup(f3);
  mrc_f3_set_comp_name(f3, 0, "fld0");
  mrc_f3_set_comp_name(f3, 1, "fld1");
  mrc_f3_view(f3);
  
  set_f3(f3);
  check_f3(f3);

  switch (testcase) {
  case 1:
    test_write_f3(f3);
    break;
  case 2:
    test_write_read_f3(f3);
    break;
  }
  mrc_f3_destroy(f3);
  mrc_domain_destroy(domain);

  MPI_Finalize();
}
