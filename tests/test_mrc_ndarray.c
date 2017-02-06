
#include <mrc_ndarray.h>

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define S3(nd, i,j,k) MRC_NDARRAY(nd, float, i,j,k,0,0)

// ----------------------------------------------------------------------
// mrc_ndarray_print_3d

static void
mrc_ndarray_print_3d(struct mrc_ndarray *nd)
{
  assert(mrc_ndarray_n_dims(nd) == 3);
  const int *offs = mrc_ndarray_offs(nd);
  const int *dims = mrc_ndarray_dims(nd);
  printf("mrc_ndarray [%d:%d,%d:%d,%d:%d]\n",
	 offs[0], offs[0] + dims[0],
	 offs[1], offs[1] + dims[1],
	 offs[2], offs[2] + dims[2]);

  struct mrc_ndarray_it it;
  mrc_ndarray_it_all(&it, nd);
  for (; !mrc_ndarray_it_done(&it); mrc_ndarray_it_next(&it)) {
    printf(" [%d,%d,%d] %g", it.idx[0], it.idx[1], it.idx[2], IT_S(&it));
    if (it.idx[0] == offs[0] + dims[0] - 1) {
      printf("\n");
      if (it.idx[1] == offs[1] + dims[1] - 1) {
	printf("\n");
      }
    }
  }
}

// ----------------------------------------------------------------------
// test_0
//
// tests that we can actually store and retrieve back values in the ndarray

static void
test_0()
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  mrc_ndarray_set_param_int3(nd, "offs", (int [3]) { 1, 2, 3 });
  mrc_ndarray_set_param_int3(nd, "dims", (int [3]) { 2, 3, 4 });
  mrc_ndarray_set_param_int3(nd, "perm", (int [3]) { 0, 1, 2 });
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);

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

  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_1
//
// like test_0, but don't set perm (tests the default)

static void
test_1()
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  mrc_ndarray_set_param_int3(nd, "offs", (int [3]) { 1, 2, 3 });
  mrc_ndarray_set_param_int3(nd, "dims", (int [3]) { 2, 3, 4 });
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);

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

  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_2
//
// like test_1, but don't set offs (tests the default)

static void
test_2()
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  mrc_ndarray_set_param_int3(nd, "dims", (int [3]) { 2, 3, 4 });
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);

  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 2; i++) {
	S3(nd, i,j,k) = i * 10000 + j * 100 + k;
      }
    }
  }

  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 2; i++) {
	assert(S3(nd, i,j,k) == i * 10000 + j * 100 + k);
      }
    }
  }

  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_3
//
// test "all" iterator (otherwise like test_1)

static void
test_3()
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  mrc_ndarray_set_param_int3(nd, "offs", (int [3]) { 1, 2, 3 });
  mrc_ndarray_set_param_int3(nd, "dims", (int [3]) { 2, 3, 4 });
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);

  struct mrc_ndarray_it it;

  for (mrc_ndarray_it_all(&it, nd); !mrc_ndarray_it_done(&it); mrc_ndarray_it_next(&it)) {
    IT_S(&it) = it.idx[0] * 10000 + it.idx[1] * 100 + it.idx[2];
  }

  mrc_ndarray_print_3d(nd);

  for (mrc_ndarray_it_all(&it, nd); !mrc_ndarray_it_done(&it); mrc_ndarray_it_next(&it)) {
    assert(IT_S(&it) == it.idx[0] * 10000 + it.idx[1] * 100 + it.idx[2]);
  }

  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_4
//
// test subset iterator (otherwise like test_3)

static void
test_4()
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  mrc_ndarray_set_param_int3(nd, "offs", (int [3]) { 1, 2, 3 });
  mrc_ndarray_set_param_int3(nd, "dims", (int [3]) { 2, 3, 4 });
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);

  struct mrc_ndarray_it it;

  mrc_ndarray_it_beg_end(&it, nd, (int [3]) { 1, 2, 4 }, (int [3]) { 3, 5, 6 });
  for (; !mrc_ndarray_it_done(&it); mrc_ndarray_it_next(&it)) {
    IT_S(&it) = it.idx[0] * 10000 + it.idx[1] * 100 + it.idx[2];
  }

  mrc_ndarray_print_3d(nd);

  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_5
//
// test set()

static void
test_5()
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  mrc_ndarray_set_param_int3(nd, "offs", (int [3]) { 1, 2, 3 });
  mrc_ndarray_set_param_int3(nd, "dims", (int [3]) { 2, 3, 4 });
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);

  mrc_ndarray_set(nd, 3.);
  mrc_ndarray_print_3d(nd);

  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// main

typedef void (*test_func)(void);

static test_func tests[] = {
  [0] = test_0,
  [1] = test_1,
  [2] = test_2,
  [3] = test_3,
  [4] = test_4,
  [5] = test_5,
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
