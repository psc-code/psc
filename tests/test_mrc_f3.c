
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>

static void
set_fld(struct mrc_fld *fld)
{
  struct mrc_crds *crds = mrc_domain_get_crds(fld->_domain);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(fld->_domain, 0, &info);
  int *off = info.off;
  mrc_fld_foreach(fld, ix,iy,iz, 0,0) {
    MRC_F3(fld, 0, ix,iy,iz) =
      (iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]);
    MRC_F3(fld, 1, ix,iy,iz) = MRC_CRD(crds, 0, ix);
  } mrc_fld_foreach_end;
}

static void
check_fld(struct mrc_fld *fld)
{
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(fld->_domain, 0, &info);
  int *off = info.off;
  mrc_fld_foreach(fld, ix,iy,iz, 0,0) {
    assert(MRC_F3(fld, 0, ix,iy,iz) ==
	   (iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]));
  } mrc_fld_foreach_end;
}

static void
test_write_fld(struct mrc_fld *fld)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(fld));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/fld", "fld", fld);
  mrc_io_close(io);

  mrc_io_open(io, "w", 1, 1.);
  mrc_io_write_path(io, "/fld", "fld", fld);
  mrc_io_close(io);

  mrc_io_destroy(io);
}

static void
test_write_read_fld(struct mrc_fld *fld)
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
  struct mrc_fld *fld_2 = mrc_io_read_path(io, "/fld", "fld", mrc_fld);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrctest_fld_compare(fld, fld_2, 0.);
  mrc_fld_destroy(fld_2);
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

  struct mrc_fld *fld = mrc_domain_fld_create(domain, 2, "fld0:fld1");
  mrc_fld_set_name(fld, "test_fld");
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);
  
  set_fld(fld);
  check_fld(fld);

  switch (testcase) {
  case 1:
    test_write_fld(fld);
    break;
  case 2:
    test_write_read_fld(fld);
    break;
  }
  mrc_fld_destroy(fld);
  mrc_domain_destroy(domain);

  MPI_Finalize();
}
