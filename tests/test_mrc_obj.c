
#include <mrctest.h>
#include <mrc_obj.h>

#include <stdio.h>

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

// ----------------------------------------------------------------------

int
main(int argc, char **argv)
{
  mrctest_init(&argc, &argv);

  struct mrc_test *test = mrc_test_create(MPI_COMM_SELF);
  struct mrc_test *test2 = mrc_test_get(test);
  mrc_test_destroy(test);
  mprintf("before destroy\n");
  mrc_test_put(test2);

  mrctest_finalize();
  return 0;
}
