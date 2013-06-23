
#include "mrctest.h"

#include <mrc_mod.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

// ----------------------------------------------------------------------

struct test_io_attrs {
  struct mrc_obj obj;
  char *run;
};

#define VAR(x) (void *)offsetof(struct test_io_attrs, x)
static struct param test_io_attrs_descr[] = {
  { "run"            , VAR(run)           , PARAM_STRING(NULL)      },
  {},
};
#undef VAR

MRC_CLASS_DECLARE(test_io_attrs, struct test_io_attrs);

struct mrc_class_test_io_attrs mrc_class_test_io_attrs = {
  .name         = "test_io_attrs",
  .size         = sizeof(struct test_io_attrs),
  .param_descr  = test_io_attrs_descr,
};

// ----------------------------------------------------------------------

static void
dump_field(struct mrc_fld *fld, int rank_diagsrv)
{
  struct mrc_domain *domain = fld->_domain;
  assert(domain);

  MPI_Comm comm = mrc_domain_comm(domain);

  struct test_io_attrs *attrs = test_io_attrs_create(comm);
  test_io_attrs_set_param_string(attrs, "run", "testrun");

  struct mrc_obj *dict = mrc_obj_create(comm);
  mrc_obj_dict_add_string(dict, "dict_test", "string");

  struct mrc_io *io;
  if (rank_diagsrv >= 0) {
    io = mrc_io_create(comm);
    mrc_io_set_type(io, "combined");
    mrc_io_set_param_int(io, "rank_diagsrv", rank_diagsrv);
  } else {
    io = mrc_io_create(comm);
  }
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);

  mrc_io_open(io, "w", 0, 0.);
  test_io_attrs_write(attrs, io);
  mrc_obj_write(dict, io);
  mrc_fld_write(fld, io);
  mrc_io_close(io);

  mrc_io_open(io, "w", 1, 1.);
  mrc_obj_write(dict, io);
  mrc_fld_write(fld, io);
  mrc_io_close(io);

  mrc_io_destroy(io);

  test_io_attrs_destroy(attrs);
  mrc_obj_destroy(dict);
}

static void
mod_domain(struct mrc_mod *mod, void *arg)
{
  struct mrctest_domain_params *par = arg;

  MPI_Comm comm = mrc_mod_get_comm(mod);
  struct mrc_domain *domain = mrctest_create_domain(comm, par);

  struct mrc_fld *fld = mrctest_create_field_1(domain);
  int rank_diagsrv = mrc_mod_get_first_node(mod, "diagsrv");
  dump_field(fld, rank_diagsrv);
  mrc_fld_destroy(fld);

  mrc_domain_destroy(domain);
}

int
main(int argc, char **argv)
{
  mrctest_init(&argc, &argv);
  mrctest_domain(mod_domain);
  mrctest_finalize();
  return 0;
}

