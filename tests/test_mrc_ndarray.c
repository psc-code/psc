
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>
#include <math.h>

#define S3(nd, i,j,k) MRC_NDARRAY(nd, float, i,j,k,0,0)

// ----------------------------------------------------------------------
// test_0
//
// tests that we can actually store and retrieve back values in the ndarray

static void
test_0()
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  mrc_fld_set_name(fld, "test_fld");
  mrc_fld_set_param_int3(fld, "offs", (int [3]) { 1, 2, 3 });
  mrc_fld_set_param_int3(fld, "dims", (int [3]) { 2, 3, 4 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  struct mrc_ndarray *nd = fld->_nd;

  for (int k = 3; k < 3 + 4; k++) {
    for (int j = 2; j < 2 + 3; j++) {
      for (int i = 1; i < 1 + 2; i++) {
	S3(nd, i,j,k) = i * 10000 + j * 100 + k;
      }
    }
  }

  for (int k = 3; k < 3 + 4; k++) {
    for (int j = 2; j < 2 + 3; j++) {
      for (int i = 1; i < 1 + 2; i++) {
	assert(S3(nd, i,j,k) == i * 10000 + j * 100 + k);
      }
    }
  }

  mrc_fld_destroy(fld);
}

// ----------------------------------------------------------------------
// main

typedef void (*test_func)(void);

static test_func tests[] = {
  [0] = test_0,
};

static int n_tests = sizeof(tests)/sizeof(tests[0]);

static void
run_test(int n)
{
  mpi_printf(MPI_COMM_WORLD, "==============================================\n");
  mpi_printf(MPI_COMM_WORLD, "TEST %d\n", n);
  mpi_printf(MPI_COMM_WORLD, "==============================================\n");

  tests[n]();

  mpi_printf(MPI_COMM_WORLD, "=== TEST %d SUCCEEDED\n", n);
}


int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int testcase = -1;
  mrc_params_get_option_int("case", &testcase);

  if (testcase == -1) {
    for (int n = 0; n < n_tests; n++) {
      run_test(n);
    }
  } else if (testcase < n_tests) {
    run_test(testcase);
  } else {
    mprintf("ERROR: invalid case %d\n", testcase);
    assert(0);
  }

  MPI_Finalize();
}
