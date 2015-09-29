
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>
#include <math.h>

// ----------------------------------------------------------------------
// mrc_fld_print_2d

static void
mrc_fld_print_2d(struct mrc_fld *fld)
{
  const int *offs = mrc_fld_ghost_offs(fld);
  const int *dims = mrc_fld_ghost_dims(fld);
  printf("mrc_fld[%d:%d,%d:%d,%d:%d]\n",
	 offs[0], offs[0] + dims[0],
	 offs[1], offs[1] + dims[1],
	 offs[2], offs[2] + dims[2]);

  for (int iz = offs[2]; iz < offs[2] + dims[2]; iz++) {
    for (int iy = offs[1]; iy < offs[1] + dims[1]; iy++) {
      for (int ix = offs[0]; ix < offs[0] + dims[0]; ix++) {
	printf(" %8g", MRC_S3(fld, ix,iy,iz));
      }
      printf("\n");
    }
    printf("\n");
  }
}

// ----------------------------------------------------------------------
// test_1

static void
test_1(int sw)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  mrc_fld_set_name(fld, "test_fld");
  mrc_fld_set_param_int3(fld, "dims", (int [3]) { 3, 4, 1 });
  mrc_fld_set_param_int3(fld, "offs", (int [3]) { 1, 2, 0 });
  mrc_fld_set_param_int3(fld, "sw", (int [3]) { sw, sw, sw });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  mrc_fld_foreach_bnd(fld, ix,iy,iz) {
    MRC_S3(fld, ix,iy,iz) = ix * 10000 + iy * 100 + iz;
  } mrc_fld_foreach_bnd_end;

  mrc_fld_foreach_bnd(fld, ix,iy,iz) {
    assert(MRC_S3(fld, ix,iy,iz) == ix * 10000 + iy * 100 + iz);
  } mrc_fld_foreach_bnd_end;

  mrc_fld_print_2d(fld);

  printf("VIEW 1:4,2:6,0:1 (identical, though shifted)\n");
  struct mrc_fld *fld2 =
    mrc_fld_create_view(fld, 3, (int [3]) { 3, 4, 1 }, (int [3]) { 1, 2, 0 });
  mrc_fld_print_2d(fld2);
  mrc_fld_destroy(fld2);

  printf("VIEW 1:3,2:5,0:1\n");
  struct mrc_fld *fld3 =
    mrc_fld_create_view(fld, 3, (int [3]) { 2, 3, 1 }, (int [3]) { 1, 2, 0 });
  mrc_fld_print_2d(fld3);
  mrc_fld_destroy(fld3);

  printf("VIEW 2:4,3:6,0:1\n");
  struct mrc_fld *fld4 =
    mrc_fld_create_view(fld, 3, (int [3]) { 2, 3, 1 }, (int [3]) { 2, 3, 0 });
  mrc_fld_print_2d(fld4);
  mrc_fld_destroy(fld4);

  printf("VIEW 2:4,3:6,0:1, set sw = 0,0,0\n");
  struct mrc_fld *fld5 =
    mrc_fld_create_view_ext(fld, 3, (int [3]) { 2, 3, 1 }, (int [3]) { 2, 3, 0},
			    NULL, NULL);
  mrc_fld_print_2d(fld5);
  mrc_fld_destroy(fld5);

  printf("VIEW 2:4,3:6,0:1 -> set sw = 0,0,0, shift by 2,3,0\n");
  struct mrc_fld *fld6 =
    mrc_fld_create_view_ext(fld, 3, (int [3]) { 2, 3, 1 }, (int [3]) { 2, 3, 0},
			    NULL, (int [3]) { 2, 3, 0 });
  mrc_fld_print_2d(fld6);
  mrc_fld_destroy(fld6);

  mrc_fld_destroy(fld);
}

// ----------------------------------------------------------------------
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int testcase = 1;
  mrc_params_get_option_int("case", &testcase);

  switch (testcase) {
  case 1: test_1(0); break;
  case 2: test_1(1); break;
  default: assert(0);
  }

  MPI_Finalize();
}
