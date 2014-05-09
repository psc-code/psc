
#include "mrctest.h"

#include <mrc_diag.h>
#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

struct test_diag_params {
  int gdims[3];
  int nproc[3];
  bool use_diagsrv;

  int ldims[3];
  int nproc_domain;
};

#define VAR(x) (void *)offsetof(struct test_diag_params, x)

static struct param test_diag_params_descr[] = {
  { "mx"              , VAR(gdims[0])        , PARAM_INT(128)        },
  { "my"              , VAR(gdims[1])        , PARAM_INT(64)         },
  { "mz"              , VAR(gdims[2])        , PARAM_INT(32)         },
  { "npx"             , VAR(nproc[0])        , PARAM_INT(1)          },
  { "npy"             , VAR(nproc[1])        , PARAM_INT(1)          },
  { "npz"             , VAR(nproc[2])        , PARAM_INT(1)          },
  { "use_diagsrv"     , VAR(use_diagsrv)     , PARAM_BOOL(false)     },
  {},
};

#undef VAR

static void
dump_field(struct mrc_fld *fld, int m, int rank_diagsrv, const char *fmt)
{
  struct mrc_domain *domain = fld->domain;
  assert(domain);

  struct diag_format *format =
    diag_format_create(((struct mrc_obj *)domain)->comm, fmt);
  if (rank_diagsrv) {
    diag_format_set_param_int(format, "rank_diagsrv", rank_diagsrv);
  }
  diag_format_set_from_options(format);
  diag_format_view(format);
  diag_format_setup(format);
  diag_format_open(format, -1, DIAG_TYPE_3D, 0, 0., "timestr");
  diag_format_write_field(format, 1., fld, m);
  diag_format_close(format);
  diag_format_destroy(format);
}

static void
init_values(struct mrc_fld *f)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->domain);

  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float xx = crds->crd[0][ix];

    MRC_FLD(f, 0, ix,iy,iz) = 2.f + .2f * sin(xx);
  } mrc_fld_foreach_end;
}

static void
do_test_diag(MPI_Comm comm, struct test_diag_params *par, int rank_diagsrv,
	     const char *fmt)
{
  struct mrc_domain_simple_params simple_par = {
    .ldims    = { par->ldims[0], par->ldims[1], par->ldims[2] },
    .nr_procs = { par->nproc[0], par->nproc[1], par->nproc[2] },
  };
  struct mrc_domain *domain = mrc_domain_create(comm, "simple");
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_domain_simple_set_params(domain, &simple_par);
  mrc_crds_set_param_int(crds, "sw", SW_2);
  mrc_crds_set_param_double(crds, "xl", -30.);
  mrc_crds_set_param_double(crds, "yl", -20.);
  mrc_crds_set_param_double(crds, "zl", -20.);
  mrc_crds_set_param_double(crds, "xh",  50.);
  mrc_crds_set_param_double(crds, "yh",  20.);
  mrc_crds_set_param_double(crds, "zh",  20.);
  mrc_domain_set_from_options(domain);
  mrc_domain_view(domain);
  mrc_domain_setup(domain);

  struct mrc_fld fld;
  mrc_domain_fld_alloc(domain, &fld, 1, SW_2);
  fld.name[0] = strdup("test");
  init_values(&fld);
  dump_field(&fld, 0, rank_diagsrv, fmt);
  mrc_fld_free(&fld);

  mrc_domain_destroy(domain);
}

static void
test_diag(struct test_diag_params *par)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  assert(par->nproc_domain <= size);
  
  MPI_Comm comm;
  int color = MPI_UNDEFINED;
  if (rank < par->nproc_domain) {
    color = 0;
  }
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
  
  if (color == 0) {
    do_test_diag(comm, par, 0, "xdmf");
  }
  if (comm != MPI_COMM_NULL) {
    MPI_Comm_free(&comm);
  }
}

static void
test_diag_diagsrv(struct test_diag_params *par)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int rank_diagsrv = par->nproc_domain;
  assert(rank_diagsrv < size);

  MPI_Comm comm;
  int color = MPI_UNDEFINED;
  if (rank < par->nproc_domain) {
    color = 0;
  } else if (rank == rank_diagsrv) {
    color = 1;
  }
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);

  if (color == 0) {
    do_test_diag(comm, par, rank_diagsrv, "combined");
  } else if (color == 1) {
    diagsrv_one("xdmf_serial", "cache", par->nproc_domain);
  }
  if (comm != MPI_COMM_NULL) {
    MPI_Comm_free(&comm);
  }
}

int
main(int argc, char **argv)
{
  mrctest_init(&argc, &argv);
  struct test_diag_params par;
  mrc_params_parse(&par, test_diag_params_descr, "test_diag", MPI_COMM_WORLD);
  mrc_params_print(&par, test_diag_params_descr, "test_diag", MPI_COMM_WORLD);

  for (int d = 0; d < 3; d++) {
    assert(par.gdims[d] % par.nproc[d] == 0);
    par.ldims[d] = par.gdims[d] / par.nproc[d];
  }
  par.nproc_domain = par.nproc[0] * par.nproc[1] * par.nproc[2];

  if (par.use_diagsrv) {
    test_diag_diagsrv(&par);
  } else {
    test_diag(&par);
  }

  mrctest_finalize();
  return 0;
}

