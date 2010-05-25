
#include "params.h"

struct test_param {
  int i1, i2;
  bool b1, b2, b3;
  double x1, x2;
  const char *s1, *s2;
};

#define VAR(x) (void *)offsetof(struct test_param, x)

static struct param test_param_descr[] = {
  { "i1"         , VAR(i1)          , PARAM_INT(1)               },
  { "i2"         , VAR(i2)          , PARAM_INT(2)               },
  { "b1"         , VAR(b1)          , PARAM_BOOL(1)              },
  { "b2"         , VAR(b2)          , PARAM_BOOL(0)              },
  { "b3"         , VAR(b3)          , PARAM_BOOL(0)              },
  { "x1"         , VAR(x1)          , PARAM_DOUBLE(1.)           },
  { "x2"         , VAR(x2)          , PARAM_DOUBLE(2.)           },
  { "s1"         , VAR(s1)          , PARAM_STRING(NULL)         },
  { "s2"         , VAR(s2)          , PARAM_STRING("test")       },
  {},
};

#undef VAR

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  params_init(argc, argv);
  params_print_all();

  struct test_param par;
  params_parse_cmdline(&par, test_param_descr, "Test", MPI_COMM_WORLD);
  params_print(&par, test_param_descr, "Test", MPI_COMM_WORLD);

  MPI_Finalize();
  return 0;
}
