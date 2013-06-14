
#include <mrctest.h>
#include <mrc_obj.h>
#include <mrc_io.h>

#include <stdio.h>
#include <assert.h>

// ----------------------------------------------------------------------

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
test_dict()
{
  struct mrc_test *test = mrc_test_create(MPI_COMM_WORLD);

  mrc_test_dict_add_int(test, "test_int", 1);
  mrc_test_dict_add_bool(test, "test_bool", true);
  mrc_test_dict_add_float(test, "test_float", 1.5);
  mrc_test_dict_add_double(test, "test_double", 2.5);
  mrc_test_dict_add_string(test, "test_string", "string");

  mrc_test_setup(test);
  mrc_test_view(test);

  struct mrc_io *io = mrc_io_create(mrc_test_comm(test));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/test", "test", test);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_test_comm(test));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_test *test2 = mrc_io_read_path(io, "/test", "test", mrc_test);
  mrc_io_close(io);
  mrc_io_destroy(io);

  union param_u *pv, *pv2;
  mrc_test_get_var(test, "test_float", &pv);
  mrc_test_get_var(test2, "test_float", &pv2);
  assert(pv->u_float == pv2->u_float);

  mrc_test_destroy(test);
  mrc_test_destroy(test2);
}

// ----------------------------------------------------------------------

int
main(int argc, char **argv)
{
  mrctest_init(&argc, &argv);

  struct mrc_test *test = mrc_test_create(MPI_COMM_SELF);
  mrc_test_set_from_options(test);
  struct mrc_test *test2 = mrc_test_get(test);
  mrc_test_destroy(test);
  mprintf("before destroy\n");
  mrc_test_put(test2);

  test_dict();

  mrctest_finalize();
  return 0;
}
