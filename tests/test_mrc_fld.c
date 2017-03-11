
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>
#include <math.h>
// ----------------------------------------------------------------------
// test_12

static void
test_12(int sw)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  mrc_fld_set_name(fld, "test_fld");
  mrc_fld_set_param_int3(fld, "offs", (int [3]) { 1, 2, 3 });
  mrc_fld_set_param_int3(fld, "dims", (int [3]) { 2, 3, 4 });
  mrc_fld_set_param_int3(fld, "sw", (int [3]) { sw, sw, sw });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    MRC_S3(fld, ix,iy,iz) = ix * 10000 + iy * 100 + iz;
  } mrc_fld_foreach_end;

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    assert(MRC_S3(fld, ix,iy,iz) == ix * 10000 + iy * 100 + iz);
  } mrc_fld_foreach_end;

  mrc_fld_destroy(fld);
}

// ----------------------------------------------------------------------
// test_34
//
// same as test_12 as "double"

static void
test_34(int sw)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  mrc_fld_set_name(fld, "test_fld");
  mrc_fld_set_type(fld, "double");
  mrc_fld_set_param_int3(fld, "offs", (int [3]) { 1, 2, 3 });
  mrc_fld_set_param_int3(fld, "dims", (int [3]) { 2, 3, 4 });
  mrc_fld_set_param_int3(fld, "sw", (int [3]) { sw, sw, sw });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    MRC_D3(fld, ix,iy,iz) = ix * 10000 + iy * 100 + iz;
  } mrc_fld_foreach_end;

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    assert(MRC_D3(fld, ix,iy,iz) == ix * 10000 + iy * 100 + iz);
  } mrc_fld_foreach_end;

  mrc_fld_destroy(fld);
}



// ----------------------------------------------------------------------
// test_56
//
// 4-d field

static void
test_56(int sw)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  mrc_fld_set_param_int_array(fld, "offs", 4, (int []) { 1, 2, 3, 0 });
  mrc_fld_set_param_int_array(fld, "dims", 4, (int []) { 2, 3, 4, 2 });
  mrc_fld_set_param_int_array(fld, "sw", 4, (int []) { sw, sw, sw, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    MRC_S4(fld, ix,iy,iz,1) = ix * 10000 + iy * 100 + iz;
  } mrc_fld_foreach_end;

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    assert(MRC_S4(fld, ix,iy,iz,1) == ix * 10000 + iy * 100 + iz);
  } mrc_fld_foreach_end;

  mrc_fld_destroy(fld);
}

// ----------------------------------------------------------------------
// test_78
//
// field on a mrc_domain

static void
test_78(int sw)
{
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);
  
  struct mrc_fld *fld = mrc_domain_fld_create(domain, sw, "test0:test1");
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    MRC_S4(fld, ix,iy,iz,1) = ix * 10000 + iy * 100 + iz;
  } mrc_fld_foreach_end;

  mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
    assert(MRC_S4(fld, ix,iy,iz,1) == ix * 10000 + iy * 100 + iz);
  } mrc_fld_foreach_end;

  mrc_fld_destroy(fld);
}

// ----------------------------------------------------------------------
// test_as
//
// test get_as/put_as

static void
test_as(int sw, const char *type, const char *as_type)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  mrc_fld_set_type(fld, type);
  mrc_fld_set_param_int_array(fld, "offs", 4, (int []) { 1, 2, 3, 0 });
  mrc_fld_set_param_int_array(fld, "dims", 4, (int []) { 2, 3, 4, 2 });
  mrc_fld_set_param_int_array(fld, "sw", 4, (int []) { sw, sw, sw, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  if (strcmp(type, "float") == 0) {
    mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
      MRC_S4(fld, ix,iy,iz,1) = ix * 10000 + iy * 100 + iz;
    } mrc_fld_foreach_end;
  } else if (strcmp(type, "double") == 0) {
    mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
      MRC_D4(fld, ix,iy,iz,1) = ix * 10000 + iy * 100 + iz;
    } mrc_fld_foreach_end;
  } else {
    assert(0);
  }

  struct mrc_fld *f = mrc_fld_get_as(fld, as_type);
  mrc_fld_view(f);

  if (strcmp(as_type, "float") == 0) {
    mrc_fld_foreach(f, ix,iy,iz, sw, sw) {
      assert(MRC_S4(f, ix,iy,iz,1) == ix * 10000 + iy * 100 + iz);
    } mrc_fld_foreach_end;
  } else if (strcmp(as_type, "double") == 0) {
    mrc_fld_foreach(f, ix,iy,iz, sw, sw) {
      assert(MRC_D4(f, ix,iy,iz,1) == ix * 10000 + iy * 100 + iz);
    } mrc_fld_foreach_end;
  } else {
    assert(0);
  }

  mrc_fld_put_as(f, fld);

  mrc_fld_destroy(fld);
}

// ----------------------------------------------------------------------
// test_fld5d

static void
test_fld5d(int sw)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  mrc_fld_set_name(fld, "test_fld");
  mrc_fld_set_type(fld, "double");
  mrc_fld_set_param_int_array(fld, "offs", 5, (int [5]) { 1, 2, 3, 4, 5 });
  mrc_fld_set_param_int_array(fld, "dims", 5, (int [5]) { 2, 3, 4, 5, 6 });
  mrc_fld_set_param_int_array(fld, "sw"  , 5, (int [5]) { sw, sw, sw, 0, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  for (int p = 5; p < 5 + 6; p++) {
    for (int m = 4; m < 4 + 5; m++) {
      mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
	MRC_D5(fld, ix,iy,iz,m,p) = p * 1e8 + m * 1e6 + ix * 10000 + iy * 100 + iz;
      } mrc_fld_foreach_end;
    }
  }

  for (int p = 5; p < 5 + 6; p++) {
    for (int m = 4; m < 4 + 5; m++) {
      mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
	assert(MRC_D5(fld, ix,iy,iz,m,p) == p * 1e8 + m * 1e6 + ix * 10000 + iy * 100 + iz);
      } mrc_fld_foreach_end;
    }
  }

  mrc_fld_destroy(fld);
}

static void
test_set(const char *type)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  int sw = 2;
  mrc_fld_set_type(fld, type);
  mrc_fld_set_param_int_array(fld, "offs", 4, (int []) { 1, 2, 3, 0 });
  mrc_fld_set_param_int_array(fld, "dims", 4, (int []) { 2, 3, 4, 2 });
  mrc_fld_set_param_int_array(fld, "sw", 4, (int []) { sw, sw, sw, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);
  
  // FIXME: should just be using double, but the set function doesn't support
  // that for some legacy reason.
  float set_val = M_PI;
  mrc_fld_set(fld, set_val);

  if (strcmp(type, "float") == 0) {
    float check_val = M_PI;
    mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
      assert(MRC_S4(fld, ix,iy,iz,1) == check_val);
    } mrc_fld_foreach_end;
  } else if (strcmp(type, "double") == 0) {
    // FIXME: should check against double, but set_val only takes float
    // so a double check is useless.
    float check_val = M_PI;
    mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
      assert(MRC_D4(fld, ix,iy,iz,1) == check_val);
    } mrc_fld_foreach_end;
  } else if (strcmp(type, "int") == 0){
    int check_val = M_PI;
    mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
      assert(MRC_FLD(fld, int, ix,iy,iz,1,0) == check_val);
    } mrc_fld_foreach_end;
  } else {
    assert(0);
  }
  mrc_fld_destroy(fld);  
}

static void
test_axpy(const char *type)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  int sw = 2;
  mrc_fld_set_type(fld, type);
  mrc_fld_set_param_int_array(fld, "offs", 4, (int []) { 1, 2, 3, 0 });
  mrc_fld_set_param_int_array(fld, "dims", 4, (int []) { 2, 3, 4, 2 });
  mrc_fld_set_param_int_array(fld, "sw", 4, (int []) { sw, sw, sw, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  struct mrc_fld *fld2 = mrc_fld_duplicate(fld);
  
  mrc_fld_set(fld, 1.0);
  mrc_fld_set(fld2, 1.0);
  
  void (*fld_axpy)(struct mrc_fld *, float, struct mrc_fld *);
  fld_axpy = (void (*)(struct mrc_fld *, float, struct mrc_fld *)) mrc_fld_get_method(fld, "axpy");
  assert(fld_axpy);
  fld_axpy(fld2, 2.0, fld);
  
  if (strcmp(type, "float") == 0) {
    float check_val = 1.0f + 2.0f * 1.0f;
    mrc_fld_foreach(fld2, ix,iy,iz, sw, sw) {
      assert(MRC_S4(fld2, ix,iy,iz,1) == check_val);
    } mrc_fld_foreach_end;
  } else if (strcmp(type, "double") == 0) {
    // Doesn't really matter in this case, but we really should be accounting
    // for casting behavior.
    double check_val = (double) 1.0f + (double) 2.0f * (double) 1.0f;
    mrc_fld_foreach(fld2, ix,iy,iz, sw, sw) {
      assert(MRC_D4(fld2, ix,iy,iz,1) == check_val);
    } mrc_fld_foreach_end;
  } else if (strcmp(type, "int") == 0){
    int check_val = 1 + 2 * 1;
    mrc_fld_foreach(fld2, ix,iy,iz, sw, sw) {
      assert(MRC_FLD(fld2, int, ix,iy,iz,1,0) == check_val);
    } mrc_fld_foreach_end;
  }
  mrc_fld_destroy(fld);  
  mrc_fld_destroy(fld2);
}

static void
test_waxpy(const char *type)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  int sw = 2;
  mrc_fld_set_type(fld, type);
  mrc_fld_set_param_int_array(fld, "offs", 4, (int []) { 1, 2, 3, 0 });
  mrc_fld_set_param_int_array(fld, "dims", 4, (int []) { 2, 3, 4, 2 });
  mrc_fld_set_param_int_array(fld, "sw", 4, (int []) { sw, sw, sw, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  struct mrc_fld *fld2 = mrc_fld_duplicate(fld);
  struct mrc_fld *fld3 = mrc_fld_duplicate(fld);
  
  mrc_fld_set(fld, 1.0);
  mrc_fld_set(fld2, 1.0);
  
  void (*fld_waxpy)(struct mrc_fld *, float, struct mrc_fld *, struct mrc_fld *);
  fld_waxpy = (void (*)(struct mrc_fld *, float, struct mrc_fld *, struct mrc_fld *)) 
    mrc_fld_get_method(fld, "waxpy");
  assert(fld_waxpy);
  fld_waxpy(fld3, 2.0, fld, fld2);
  
  if (strcmp(type, "float") == 0) {
    float check_val = 1.0f + 2.0f * 1.0f;
    mrc_fld_foreach(fld3, ix,iy,iz, sw, sw) {
      assert(MRC_S4(fld3, ix,iy,iz,1) == check_val);
    } mrc_fld_foreach_end;
  } else if (strcmp(type, "double") == 0) {
    // Doesn't really matter in this case, but we really should be accounting
    // for casting behavior.
    double check_val = (double) 1.0f + (double) 2.0f * (double) 1.0f;
    mrc_fld_foreach(fld3, ix,iy,iz, sw, sw) {
      assert(MRC_D4(fld3, ix,iy,iz,1) == check_val);
    } mrc_fld_foreach_end;
  } else if (strcmp(type, "int") == 0){
    int check_val = 1 + 2 * 1;
    mrc_fld_foreach(fld3, ix,iy,iz, sw, sw) {
      assert(MRC_FLD(fld3, int, ix,iy,iz,1,0) == check_val);
    } mrc_fld_foreach_end;
  }
  mrc_fld_destroy(fld);  
  mrc_fld_destroy(fld2);
  mrc_fld_destroy(fld3);
}

// ----------------------------------------------------------------------
// test_make_view
//
// view of an mrc_fld on a mrc_domain

static void
test_make_view(int sw)
{
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);
  
  struct mrc_fld *fld = mrc_domain_fld_create(domain, sw, "test0:test1");
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  for (int m = 0; m < 2; m++) {
    mrc_fld_foreach(fld, ix,iy,iz, sw, sw) {
      MRC_S4(fld, ix,iy,iz, m) = m * 1000 + ix * 100 + iy * 10 + iz;
    } mrc_fld_foreach_end;
  }

  struct mrc_fld *fld2 = mrc_fld_make_view(fld, 1, 2);
  
  mrc_fld_foreach(fld2, ix,iy,iz, sw, sw) {
    //mprintf("[%d,%d,%d,%d] = %g\n", 0, ix, iy, iz, MRC_S4(fld2, ix,iy,iz, 0));
    assert(MRC_S4(fld2, ix,iy,iz, 0) == 1 * 1000 + ix * 100 + iy * 10 + iz);
  } mrc_fld_foreach_end;
  
  mrc_fld_destroy(fld2);
  mrc_fld_destroy(fld);
}
  
static void
test_norm(const char *type)
{
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_WORLD);
  int sw = 2;
  mrc_fld_set_type(fld, type);
  mrc_fld_set_param_int_array(fld, "offs", 4, (int []) { 1, 2, 3, 0 });
  mrc_fld_set_param_int_array(fld, "dims", 4, (int []) { 2, 3, 4, 2 });
  mrc_fld_set_param_int_array(fld, "sw", 4, (int []) { sw, sw, sw, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_view(fld);

  mrc_fld_set(fld, 3.);
  
  float norm = mrc_fld_norm(fld);
  assert(norm == 3.);
  
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
  case 1: test_12(0); break;
  case 2: test_12(1); break;
  case 3: test_34(0); break;
  case 4: test_34(1); break;
  case 5: test_56(0); break;
  case 6: test_56(1); break;
  case 7: test_78(0); break;
  case 8: test_78(1); break;
  case 9: test_as(1, "float", "float"); break;
  case 10: test_as(1, "float", "double"); break;
  case 11: test_as(1, "double", "float"); break;
  case 12: test_fld5d(0); break;
  case 13: test_fld5d(1); break;
  case 14: test_set("float"); break;
  case 15: test_set("double"); break;
  case 16: test_set("int"); break;
  case 17: test_axpy("float"); break;
  case 18: test_axpy("double"); break;
  case 19: test_axpy("int"); break;
  case 20: test_waxpy("float"); break;
  case 21: test_waxpy("double"); break;
  case 22: test_waxpy("int"); break;
  case 23: test_make_view(0); break;
  case 24: test_norm("float"); break;
  case 25: test_norm("double"); break;
  case 26: test_norm("int"); break;
  default: assert(0);
  }

  MPI_Finalize();
}
