
#include <mrc_dict.h>
#include <mrctest.h>

int
main(int argc, char **argv)
{
  mrctest_init(&argc, &argv);

  struct mrc_dict *dict = mrc_dict_create(MPI_COMM_WORLD);
  mrc_dict_add_int(dict, "test_int", 1);
  mrc_dict_add_bool(dict, "test_bool", true);
  mrc_dict_add_float(dict, "test_float", 1.5);
  mrc_dict_add_double(dict, "test_double", 2.5);
  mrc_dict_add_string(dict, "test_string", "string");

  mrc_dict_view(dict);
  mrc_dict_destroy(dict);

  mrctest_finalize();
  return 0;
}

