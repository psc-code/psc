
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>

static void
set_m3(struct mrc_fld *m3)
{
  struct mrc_patch *patches = mrc_domain_get_patches(m3->_domain, NULL);
  struct mrc_crds *crds = mrc_domain_get_crds(m3->_domain);

  mrc_fld_foreach_patch(m3, p) {
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    int *off = patches[p].off;
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      MRC_M3(m3p, 0, ix,iy,iz) =
	(iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]);
      MRC_M3(m3p, 1, ix,iy,iz) = MRC_MCRD(crds, 0, ix, p);
    } mrc_m3_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

static void
check_m3(struct mrc_fld *m3)
{
  struct mrc_patch *patches = mrc_domain_get_patches(m3->_domain, NULL);

  mrc_fld_foreach_patch(m3, p) {
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    int *off = patches[p].off;
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      assert(MRC_M3(m3p, 0, ix,iy,iz) ==
	     (iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]));
    } mrc_m3_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

static void
test_write_m3(struct mrc_fld *m3)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(m3));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/m3", "m3", m3);
  mrc_io_close(io);

  mrc_io_open(io, "w", 1, 1.);
  mrc_io_write_path(io, "/m3", "m3", m3);
  mrc_io_close(io);

  mrc_io_destroy(io);
}

static void
test_write_read_m3(struct mrc_fld *m3)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(m3));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/m3", "m3", m3);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_fld_comm(m3));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_fld *m3_2 = mrc_io_read_path(io, "/m3", "m3", mrc_fld);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrctest_m3_compare(m3, m3_2);
  mrc_fld_destroy(m3_2);
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
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);
  
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  if (strcmp(mrc_crds_type(crds), "rectilinear") == 0) {
    mrctest_set_crds_rectilinear_1(domain);
  }

  struct mrc_fld *m3 = mrc_domain_m3_create(domain);
  mrc_fld_set_name(m3, "test_m3");
  mrc_fld_set_param_int(m3, "nr_comps", 2);
  mrc_fld_set_from_options(m3);
  mrc_fld_setup(m3);
  mrc_fld_set_comp_name(m3, 0, "fld0");
  mrc_fld_set_comp_name(m3, 1, "fld1");
  mrc_fld_view(m3);
  
  set_m3(m3);
  check_m3(m3);

  switch (testcase) {
  case 1:
    test_write_m3(m3);
    break;
  case 2:
    test_write_read_m3(m3);
    break;
  }
  mrc_fld_destroy(m3);
  mrc_domain_destroy(domain);

  MPI_Finalize();
}
