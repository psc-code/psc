
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_crds_gen.h>
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrctest.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>

// ----------------------------------------------------------------------
// print_x_crds

static void
print_x_crds(struct mrc_crds *crds)
{
  struct mrc_fld *crdx = crds->crd[0];
  int m = mrc_fld_dims(crdx)[0];
  
  for (int p = 0; p < mrc_fld_nr_patches(crdx); p++) {
    printf("patch %d:\n", p);
    printf("MRC_DMCRDX:");
    for (int i = 0; i < m; i++) {
      printf(" %g", MRC_DMCRDX(crds, i, p));
    }
    printf("\n");

    printf("MRC_DMCRDX_NC:");
    for (int i = 0; i <= m; i++) {
      printf(" %g", MRC_DMCRDX_NC(crds, i, p));
    }
    printf("\n");
  }
}

// ----------------------------------------------------------------------
// norm_test

static void
norm_test(double norm_length, double norm_length_scale)
{
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_param_int3(domain, "m", (int [3]) { 8, 8, 8 });
			  
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_double3(crds, "l", (double [3]) { -1., -2., -3. });
  mrc_crds_set_param_double3(crds, "h", (double [3]) {  1.,  2.,  3. });
  mrc_crds_set_param_double(crds, "norm_length", norm_length);
  mrc_crds_set_param_double(crds, "norm_length_scale", norm_length_scale);

  mrc_domain_setup(domain);
  mrc_domain_view(domain);

  print_x_crds(crds);

  mrc_domain_destroy(domain);
}

// ----------------------------------------------------------------------
// norm_test_non_uniform

static void
norm_test_non_uniform(double norm_length, double norm_length_scale)
{
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_param_int3(domain, "m", (int [3]) { 8, 8, 8 });
			  
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "rectilinear");
  mrc_crds_set_param_double3(crds, "l", (double [3]) { -1., -2., -3. });
  mrc_crds_set_param_double3(crds, "h", (double [3]) {  1.,  2.,  3. });
  mrc_crds_set_param_double(crds, "norm_length", norm_length);
  mrc_crds_set_param_double(crds, "norm_length_scale", norm_length_scale);
  struct mrc_crds_gen *gen_x = crds->crds_gen[0];
  mrc_crds_gen_set_type(gen_x, "ggcm_yz");
  mrc_crds_gen_set_param_double(gen_x, "center_spacing", .2);

  mrc_domain_setup(domain);
  mrc_domain_view(domain);

  print_x_crds(crds);

  mrc_domain_destroy(domain);
}

// ----------------------------------------------------------------------
// test_0
//
// test basic coordinates without scaling

static void
test_0()
{
  norm_test(1., 1.);
}

// ----------------------------------------------------------------------
// test_1
//
// test basic coordinates where the "l", "h" are now in "km".

static void
test_1()
{
  norm_test(1., 1000.);
}

// ----------------------------------------------------------------------
// test_2
//
// test basic coordinates where the "l", "h" are now in R_E,
// but code units are normalized to R_E, too

static void
test_2()
{
  const double R_E = 6370e3;
  norm_test(R_E, R_E);
}

// ----------------------------------------------------------------------
// test_3
//
// test non-unform coordinates

static void
test_3()
{
  norm_test_non_uniform(1., 1.);
}

// ----------------------------------------------------------------------
// test_4
//
// test non-unform coordinates with rescaled coordinates

static void
test_4()
{
  norm_test_non_uniform(1., 1000.);
}

// ----------------------------------------------------------------------
// tests

typedef void (*test_func)(void);

static test_func tests[] = {
  [0] = test_0,
  [1] = test_1,
  [2] = test_2,
  [3] = test_3,
  [4] = test_4,
};

static int n_tests = sizeof(tests)/sizeof(tests[0]);

// ----------------------------------------------------------------------
// run_test

static void
run_test(int n)
{
  mpi_printf(MPI_COMM_WORLD, "\n=== TEST %d\n", n);

  tests[n]();

  mpi_printf(MPI_COMM_WORLD, "=== TEST %d SUCCEEDED\n", n);
}

// ----------------------------------------------------------------------
// main

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
