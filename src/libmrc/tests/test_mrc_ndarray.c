
#include <mrc_ndarray.h>

#include <mrc_io.h>

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define S3(nd, i,j,k) MRC_NDARRAY(nd, float, i,j,k,0,0)

// ----------------------------------------------------------------------
// setup_nd

static struct mrc_ndarray *
setup_nd(const int *offs, const int *dims, const int *perm)
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  if (offs) {
    mrc_ndarray_set_param_int3(nd, "offs", offs);
  }
  if (dims) {
    mrc_ndarray_set_param_int3(nd, "dims", dims);
  }
  if (perm) {
    mrc_ndarray_set_param_int3(nd, "perm", perm);
  }
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);

  return nd;
}

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
    printf(" [%d,%d,%d] %03g", it.idx[0], it.idx[1], it.idx[2], IT_S(&it));
    if (it.idx[0] == offs[0] + dims[0] - 1) {
      printf("\n");
      if (it.idx[1] == offs[1] + dims[1] - 1) {
	printf("\n");
      }
    }
  }
}

// ----------------------------------------------------------------------
// set_and_assert_3d

static void
set_and_assert_3d(struct mrc_ndarray *nd)
{
  int *offs = mrc_ndarray_offs(nd), *dims = mrc_ndarray_dims(nd);

  for (int k = offs[2]; k < offs[2] + dims[2]; k++) {
    for (int j = offs[1]; j < offs[1] + dims[1]; j++) {
      for (int i = offs[0]; i < offs[0] + dims[0]; i++) {
	S3(nd, i,j,k) = i * 10000 + j * 100 + k;
      }
    }
  }

  for (int k = offs[2]; k < offs[2] + dims[2]; k++) {
    for (int j = offs[1]; j < offs[1] + dims[1]; j++) {
      for (int i = offs[0]; i < offs[0] + dims[0]; i++) {
	assert(S3(nd, i,j,k) == i * 10000 + j * 100 + k);
      }
    }
  }
}

// ----------------------------------------------------------------------
// setup_and_set_nd

static struct mrc_ndarray *
setup_and_set_nd(int *offs, int *dims, int *perm)
{
  struct mrc_ndarray *nd = setup_nd(offs, dims, perm);

  struct mrc_ndarray_it it;
  for (mrc_ndarray_it_all(&it, nd); !mrc_ndarray_it_done(&it); mrc_ndarray_it_next(&it)) {
    IT_S(&it) = it.idx[0] * 100 + it.idx[1] * 10 + it.idx[2];
  }

  return nd;
}

// ----------------------------------------------------------------------
// make_view

struct mrc_ndarray *
make_view(struct mrc_ndarray *nd_base, const int *offs, const int *dims,
	  const int *view_offs)
{
  struct mrc_ndarray *nd = mrc_ndarray_create(MPI_COMM_WORLD);
  if (offs) {
    mrc_ndarray_set_param_int3(nd, "offs", offs);
  }
  if (dims) {
    mrc_ndarray_set_param_int3(nd, "dims", dims);
  }
  if (view_offs) {
    mrc_ndarray_set_param_int3(nd, "view_offs", view_offs);
  }
  mrc_ndarray_set_param_obj(nd, "view_base", nd_base);
  mrc_ndarray_set_from_options(nd);
  mrc_ndarray_setup(nd);
  mrc_ndarray_view(nd);
  return nd;
}		   

// ----------------------------------------------------------------------
// test_0
//
// tests that we can actually store and retrieve back values in the ndarray

static void
test_0()
{
  struct mrc_ndarray *nd = setup_nd((int [3]) { 1, 2, 3 }, (int [3]) { 2, 3, 4 },
				    (int [3]) { 0, 1, 2 });
  set_and_assert_3d(nd);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_1
//
// like test_0, but don't set perm (tests the default)

static void
test_1()
{
  struct mrc_ndarray *nd = setup_nd((int [3]) { 1, 2, 3 }, (int [3]) { 2, 3, 4 },
				    NULL);
  set_and_assert_3d(nd);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_2
//
// like test_1, but don't set offs (tests the default)

static void
test_2()
{
  struct mrc_ndarray *nd = setup_nd(NULL, (int [3]) { 2, 3, 4 }, NULL);

  set_and_assert_3d(nd);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_3
//
// test "all" iterator (otherwise like test_1)

static void
test_3()
{
  struct mrc_ndarray *nd = setup_nd((int [3]) { 1, 2, 3 }, (int [3]) { 2, 3, 4 },
				    NULL);

  struct mrc_ndarray_it it;

  for (mrc_ndarray_it_all(&it, nd); !mrc_ndarray_it_done(&it); mrc_ndarray_it_next(&it)) {
    IT_S(&it) = it.idx[0] * 10000 + it.idx[1] * 100 + it.idx[2];
  }

  mrc_ndarray_print_3d(nd);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_4
//
// test subset iterator (otherwise like test_3)

static void
test_4()
{
  struct mrc_ndarray *nd = setup_nd((int [3]) { 1, 2, 3 }, (int [3]) { 2, 3, 4 },
				    NULL);

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
  struct mrc_ndarray *nd = setup_nd((int [3]) { 1, 2, 3 }, (int [3]) { 2, 3, 4 },
				    NULL);
  mrc_ndarray_set(nd, 3.);
  mrc_ndarray_print_3d(nd);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_6
//
// test view of same size

static void
test_6()
{
  struct mrc_ndarray *nd = setup_and_set_nd((int [3]) { 1, 2, 0 }, (int [3]) { 3, 4, 1 },
					    NULL);

  printf("VIEW 1:4,2:6,0:1 (identical)\n");
  struct mrc_ndarray *nd_view =
    make_view(nd, (int [3]) { 1, 2, 0 }, (int [3]) { 3, 4, 1 }, NULL);
  
  mrc_ndarray_print_3d(nd);
  mrc_ndarray_print_3d(nd_view);

  S3(nd_view, 1, 2, 0) = 999;

  printf("first element set to 999 (should be in both ndarrays)\n\n");
  mrc_ndarray_print_3d(nd);
  mrc_ndarray_print_3d(nd_view);

  mrc_ndarray_destroy(nd_view);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_7
//
// test view that's actually a subset

static void
test_7()
{
  struct mrc_ndarray *nd = setup_and_set_nd((int [3]) { 1, 2, 0 }, (int [3]) { 3, 4, 1 },
					    NULL);
  printf("VIEW 1:3,2:5,0:1\n");
  struct mrc_ndarray *nd_view =
    make_view(nd, (int [3]) { 1, 2, 0 }, (int [3]) { 2, 3, 1 }, NULL);

  mrc_ndarray_print_3d(nd);
  mrc_ndarray_print_3d(nd_view);

  mrc_ndarray_destroy(nd_view);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_8
//
// test view that's actually a subset, starting not at the same corner

static void
test_8()
{
  struct mrc_ndarray *nd = setup_and_set_nd((int [3]) { 1, 2, 0 }, (int [3]) { 3, 4, 1 },
					    NULL);
  printf("VIEW 2:4,3:5,0:1\n");
  struct mrc_ndarray *nd_view =
    make_view(nd, (int [3]) { 2, 3, 0 }, (int [3]) { 2, 3, 1 }, NULL);

  mrc_ndarray_print_3d(nd);
  mrc_ndarray_print_3d(nd_view);

  mrc_ndarray_destroy(nd_view);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_9
//
// test subset view with shifted indices

static void
test_9()
{
  struct mrc_ndarray *nd = setup_and_set_nd((int [3]) { 1, 2, 0 }, (int [3]) { 3, 4, 1 },
					    NULL);
  printf("VIEW 0:2,0:3,0:1\n");
  struct mrc_ndarray *nd_view =
    make_view(nd, NULL, (int [3]) { 2, 3, 1 }, (int [3]) { 1, 3, 0 });

  mrc_ndarray_print_3d(nd);
  mrc_ndarray_print_3d(nd_view);

  mrc_ndarray_destroy(nd_view);
  mrc_ndarray_destroy(nd);
}

// ----------------------------------------------------------------------
// test_10
//
// test direct copy (no view involved)

static void
test_10()
{
  struct mrc_ndarray *nd = setup_and_set_nd((int [3]) { 1, 2, 0 }, (int [3]) { 3, 4, 1 },
					    NULL);
  struct mrc_ndarray *nd2 = setup_and_set_nd(NULL, (int [3]) { 3, 4, 1 }, NULL);

  printf("nd\n");
  mrc_ndarray_print_3d(nd);
  printf("nd2\n");
  mrc_ndarray_print_3d(nd2);

  mrc_ndarray_copy(nd, nd2);

  printf("nd after mrc_ndarray_copy(nd, nd2)\n");
  mrc_ndarray_print_3d(nd);
}

// ----------------------------------------------------------------------
// test_11
//
// test copy into view

static void
test_11()
{
  struct mrc_ndarray *nd = setup_and_set_nd(NULL, (int [3]) { 3, 4, 1 }, NULL);
  printf("VIEW 1:2,1:3,0:1\n");
  struct mrc_ndarray *nd_view =
    make_view(nd, NULL, (int [3]) { 1, 2, 1 }, (int [3]) { 1, 1, 0 });

  struct mrc_ndarray *nd2 = setup_and_set_nd(NULL, (int [3]) { 1, 2, 1 }, NULL);

  mrc_ndarray_print_3d(nd);
  mrc_ndarray_print_3d(nd_view);

  printf("nd2 0:1,0:2,0:1\n");
  mrc_ndarray_print_3d(nd2);

  mrc_ndarray_copy(nd_view, nd2);

  printf("after copy\n");
  mrc_ndarray_print_3d(nd);
  mrc_ndarray_print_3d(nd_view);
  
  mrc_ndarray_destroy(nd_view);
  mrc_ndarray_destroy(nd);
  mrc_ndarray_destroy(nd2);
}

// ----------------------------------------------------------------------
// test 12
//
// test f_contiguous

static void
test_12()
{
  struct mrc_ndarray *nd = setup_and_set_nd(NULL, (int [3]) { 3, 4, 1 }, NULL);
  printf("VIEW 1:2,1:3,0:1\n");
  struct mrc_ndarray *nd_view =
    make_view(nd, NULL, (int [3]) { 1, 2, 1 }, (int [3]) { 1, 1, 0 });

  assert(mrc_ndarray_f_contiguous(nd));
  assert(!mrc_ndarray_f_contiguous(nd_view));
}

// ----------------------------------------------------------------------
// test 13
//
// test write / read

static void
test_13()
{
  struct mrc_ndarray *nd = setup_and_set_nd(NULL, (int [3]) { 3, 4, 1 }, NULL);

  // write it
  struct mrc_io *io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_type(io, "hdf5_serial");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/nd", "nd", nd);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_ndarray_print_3d(nd);
  mrc_ndarray_destroy(nd);

  io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_type(io, "hdf5_serial");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  nd = mrc_io_read_path(io, "/nd", "nd", mrc_ndarray);
  mrc_io_close(io);
  mrc_io_destroy(io);

  printf("read back:\n");
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
  [6] = test_6,
  [7] = test_7,
  [8] = test_8,
  [9] = test_9,
  [10] = test_10,
  [11] = test_11,
  [12] = test_12,
  [13] = test_13,
};

static int n_tests = sizeof(tests)/sizeof(tests[0]);

static void
run_test(int n)
{
  mpi_printf(MPI_COMM_WORLD, "\n=== TEST %d\n", n);

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
