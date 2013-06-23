
#include <mrctest.h>
#include <mrc_obj.h>
#include <mrc_io.h>

#include <stdio.h>
#include <assert.h>

// ----------------------------------------------------------------------
// test_1

struct mrc_test {
  struct mrc_obj obj;
};

static void
_mrc_test_destroy(struct mrc_test *test)
{
  mprintf("destroy called\n");
}

MRC_CLASS_DECLARE(mrc_test, struct mrc_test);

struct mrc_class_mrc_test mrc_class_mrc_test = {
  .name         = "mrc_test",
  .size         = sizeof(struct mrc_test),
  .destroy      = _mrc_test_destroy,
};

static void
test_1()
{
  struct mrc_test *test = mrc_test_create(MPI_COMM_SELF);
  mrc_test_set_from_options(test);
  struct mrc_test *test2 = mrc_test_get(test);
  mrc_test_destroy(test);
  mprintf("before destroy\n");
  mrc_test_put(test2);
}

// ----------------------------------------------------------------------
// test_2

static void
test_2()
{
  struct mrc_obj *test = mrc_obj_create(MPI_COMM_WORLD);

  mrc_obj_dict_add_int(test, "test_int", 1);
  mrc_obj_dict_add_bool(test, "test_bool", true);
  mrc_obj_dict_add_float(test, "test_float", 1.5);
  mrc_obj_dict_add_double(test, "test_double", 2.5);
  mrc_obj_dict_add_string(test, "test_string", "string");

  mrc_obj_setup(test);
  mrc_obj_view(test);

  struct mrc_io *io = mrc_io_create(mrc_obj_comm(test));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/test", "test", test);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_obj_comm(test));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_obj *test2 = mrc_io_read_path(io, "/test", "test", mrc_obj);
  mrc_io_close(io);
  mrc_io_destroy(io);

  union param_u *pv, *pv2;
  mrc_obj_get_var(test, "test_float", &pv);
  mrc_obj_get_var(test2, "test_float", &pv2);
  assert(pv->u_float == pv2->u_float);

  mrc_obj_destroy(test);
  mrc_obj_destroy(test2);
}

// ----------------------------------------------------------------------
// test_3

struct mrc_test3 {
  struct mrc_obj obj;
  struct mrc_domain *domain;
};

MRC_CLASS_DECLARE(mrc_test3, struct mrc_test3);

#define VAR(x) (void *)offsetof(struct mrc_test3, x)
static struct param mrc_test3_descr[] = {
  { "domain"          , VAR(domain)          , MRC_VAR_OBJ(mrc_domain)     },
  {},
};
#undef VAR

struct mrc_class_mrc_test3 mrc_class_mrc_test3 = {
  .name         = "mrc_test3",
  .size         = sizeof(struct mrc_test3),
  .param_descr  = mrc_test3_descr,
};

static void
test_3()
{
  struct mrc_test3 *test = mrc_test3_create(MPI_COMM_SELF);
  mrc_test3_set_from_options(test);
  mrc_test3_setup(test);
  mrc_test3_view(test);

  struct mrc_io *io = mrc_io_create(mrc_test3_comm(test));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/test", "test", test);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_test3_destroy(test);
}

// ----------------------------------------------------------------------
// main

int
main(int argc, char **argv)
{
  mrctest_init(&argc, &argv);

  int testcase = 1;
  mrc_params_get_option_int("case", &testcase);

  switch (testcase) {
  case 1: test_1(); break;
  case 2: test_2(); break;
  case 3: test_3(); break;
  default: assert(0);
  }

  mrctest_finalize();
  return 0;
}
