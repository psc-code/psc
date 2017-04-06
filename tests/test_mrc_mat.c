
#include <mrc_params.h>
#include <mrc_fld.h>
#include <mrc_fld_as_double.h>
#include <mrc_mat.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

static void
mrc_fld_print(struct mrc_fld *x, const char *name)
{
  MPI_Comm comm = mrc_fld_comm(x);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  for (int r = 0; r < size; r++) {
    MPI_Barrier(comm);
    if (r == rank) {
      for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
	mrc_fld_foreach(x, i,j,k, 0, 0) {
	  mprintf("%s[%d,%d,%d, %d] = %g\n", name, i, j, k, p, M3(x, 0, i,j,k, p));
	} mrc_fld_foreach_end;
	mprintf("\n");
      }
    }
  }
}

struct entry {
  int col;
  mrc_fld_data_t val;
};

const int N_0 = 8;
static const struct entry *mat_0[8] = {
  [0] = (struct entry[]) {             { 0, -2. }, { 1, 1. }, { -1, } },
  [1] = (struct entry[]) { { 0,  1. }, { 1, -2. }, { 2, 1. }, { -1, } },
  [2] = (struct entry[]) { { 1,  1. }, { 2, -2. }, { 3, 1. }, { -1, } },
  [3] = (struct entry[]) { { 2,  1. }, { 3, -2. }, { 4, 1. }, { -1, } },
  [4] = (struct entry[]) { { 3,  1. }, { 4, -2. }, { 5, 1. }, { -1, } },
  [5] = (struct entry[]) { { 4,  1. }, { 5, -2. }, { 6, 1. }, { -1, } },
  [6] = (struct entry[]) { { 5,  1. }, { 6, -2. }, { 7, 1. }, { -1, } },
  [7] = (struct entry[]) { { 6,  1. }, { 7, -2. },            { -1, } },
};

const int N_1 = 12;
static const struct entry *mat_1[12] = {
  [ 0] = (struct entry[]) { {  0, 1. },                        { 9, 9. }, { -1, } },
  [ 1] = (struct entry[]) { {  1, 1. },             { 8, 8. },            { -1, } },
  [ 2] = (struct entry[]) { {  2, 1. }, {  6, 6. },                       { -1, } },
  [ 3] = (struct entry[]) { {  3, 1. }, { -1, } },
  [ 4] = (struct entry[]) { {  4, 1. }, { -1, } },
  [ 5] = (struct entry[]) { {  5, 1. }, { -1, } },
  [ 6] = (struct entry[]) { {  6, 1. }, { -1, } },
  [ 7] = (struct entry[]) { {  7, 1. }, { -1, } },
  [ 8] = (struct entry[]) { {  8, 1. }, { -1, } },
  [ 9] = (struct entry[]) { {  9, 1. }, { -1, } },
  [10] = (struct entry[]) { { 10, 1. }, { -1, } },
  [11] = (struct entry[]) { { 11, 1. }, { -1, } },
};

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int testcase = 0;
  mrc_params_get_option_int("testcase", &testcase);

  int N;
  const struct entry **mat;

  switch (testcase) {
  case 0: N = N_0; mat = mat_0; break;
  case 1: N = N_1; mat = mat_1; break;
  default: assert(0);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  struct mrc_domain *domain = mrc_domain_create(comm);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_param_int3(domain, "m", (int [3]) { N, 1, 1});
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);
  
  struct mrc_fld *x = mrc_domain_fld_create(domain, 0, "x0");
  mrc_fld_set_type(x, FLD_TYPE);
  mrc_fld_setup(x);
  mrc_fld_view(x);
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    mrc_fld_foreach(x, i,j,k, 0, 0) {
      M3(x, 0, i,0,0, p) = i;
    } mrc_fld_foreach_end;
  }

  struct mrc_fld *y = mrc_domain_fld_create(domain, 0, "y0");
  mrc_fld_set_type(y, FLD_TYPE);
  mrc_fld_setup(y);

  struct mrc_mat *A = mrc_mat_create(comm);
  mrc_mat_set_type(A, "csr_mpi");
  mrc_mat_set_param_int(A, "m", mrc_fld_len(y));
  mrc_mat_set_param_int(A, "n", mrc_fld_len(x));
  mrc_mat_set_from_options(A);
  mrc_mat_setup(A);
  mrc_mat_view(A);

  // FIXME, this really should be on offsets in mrc_vec,
  // not making assumptions about an underlying domain, equal
  // patches, etc.
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(domain, 0, &info);
  int row_off = mrc_fld_len(x) * info.global_patch;
  
  for (int i = 0; i < mrc_fld_len(x); i++) {
    MRC_D1(x, i) = i + row_off;
  }
  
  for (int i = 0; i < mrc_fld_len(x); i++) {
    int row_idx = i + row_off;
    
    for (const struct entry *e = mat[row_idx]; e && e->col >= 0; e++) {
      /* mprintf("row %d col %d val %g\n", row_idx, e->col, e->val); */
      mrc_mat_add_value(A, row_idx, e->col, e->val);
    }
  }
  
  //////
  int rank;
  MPI_Comm_rank(mrc_fld_comm(x), &rank);
  int _r = 2;
  // this test to see if _r is local only works if y is evenly
  // divided among procs
  if (testcase == 0 && rank * mrc_fld_len(y) <= _r && _r < (rank + 1) * mrc_fld_len(y)) {
    // test adding to rows/cols out of order
    // mprintf("> test adding values out of order\n");
    mrc_mat_add_value(A, _r, 0, 10.0);
    mrc_mat_add_value(A, _r, 1, 11.0);
    mrc_mat_add_value(A, _r, 2, 12.0);
    mrc_mat_add_value(A, _r, 3, 13.0);
    mrc_mat_add_value(A, _r, 4, 14.0);
    mrc_mat_add_value(A, _r, 5, 15.0);
    mrc_mat_add_value(A, _r, 6, 16.0);
    mrc_mat_add_value(A, _r, 7, 17.0);
    // mrc_mat_print(A);
    // test removing values
    // mprintf("> test removing values\n");
    mrc_mat_add_value(A, _r, 1, -11.0);
    mrc_mat_add_value(A, _r, 3, -13.0);
    mrc_mat_add_value(A, _r, 5, -15.0);
    mrc_mat_add_value(A, _r, 7, -17.0);
    // mrc_mat_print(A);
    mrc_mat_add_value(A, _r, 0, -10.0);
    mrc_mat_add_value(A, _r, 2, -12.0);
    mrc_mat_add_value(A, _r, 4, -14.0);
    mrc_mat_add_value(A, _r, 6, -16.0);
  }
  //////
  
  mrc_mat_assemble(A);
  mrc_mat_print(A);

  mrc_fld_print(x, "x");

  mrc_mat_apply(y->_nd->vec, A, x->_nd->vec);

  MPI_Barrier(comm);
  mrc_fld_print(y, "y");

  mrc_fld_destroy(x);
  mrc_fld_destroy(y);
  mrc_mat_destroy(A);

  MPI_Finalize();
}
