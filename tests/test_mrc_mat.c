
#include <mrc_params.h>
#include <mrc_fld.h>
#include <mrc_fld_as_double.h>
#include <mrc_mat.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

const int N = 7;

static void
mrc_fld_print(struct mrc_fld *x, const char *name)
{
  for (int i = 0; i < N; i++) {
    mprintf("%s[%d] = %g\n", name, i, F3(x, 0, i,0,0));
  }
  mprintf("\n");
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int testcase = 0;
  mrc_params_get_option_int("case", &testcase);

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
  for (int i = 0; i < N; i++) {
    F3(x, 0, i,0,0) = i;
  }

  struct mrc_fld *y = mrc_domain_fld_create(domain, 0, "y0");
  mrc_fld_set_type(y, FLD_TYPE);
  mrc_fld_setup(y);

  struct mrc_mat *A = mrc_mat_create(comm);
  mrc_mat_set_param_int(A, "m", N);
  mrc_mat_set_param_int(A, "n", N);
  mrc_mat_setup(A);
  mrc_mat_view(A);

  for (int i = 0; i < N; i++) {
    mrc_mat_add_value(A, i, i, -2.);
  }
  mrc_mat_assemble(A);

  mrc_fld_print(x, "x");

  mrc_mat_apply(y, A, x);

  mrc_fld_print(y, "y");

  mrc_fld_destroy(x);
  mrc_mat_destroy(A);

  MPI_Finalize();
}
