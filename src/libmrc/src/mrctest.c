
#include "mrctest.h"

#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_params.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <mrc_mod.h>

#include <math.h>
#include <assert.h>
#include <string.h>

void
mrc_f3_init_values(struct mrc_f3 *f3, struct mrc_f3_init_values_info *iv_info)
{
  struct mrc_domain *domain = f3->domain;
  assert(domain);
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  assert(crds);

  mrc_f3_set(f3, 0.);

  for (int i = 0; iv_info->ini_flds[i].ini; i++) {
    int m = iv_info->ini_flds[i].m;
    mrc_f3_foreach(f3, ix,iy,iz, 0,0) {
      float xx = MRC_CRDX(crds, ix), yy = MRC_CRDY(crds, iy), zz = MRC_CRDZ(crds, iz);
      MRC_F3(f3,m, ix,iy,iz) = iv_info->ini_flds[i].ini(xx, yy, zz);
    } mrc_f3_foreach_end;
  }
}

void
mrctest_init(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
  libmrc_params_init(*argc, *argv);
}

void
mrctest_finalize()
{
  prof_print();
  MPI_Finalize();
}

// ----------------------------------------------------------------------
// mrctest_domain

#define VAR(x) (void *)offsetof(struct mrctest_domain_params, x)
static struct param mrctest_domain_params_descr[] = {
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

void
mrctest_domain_init(struct mrctest_domain_params *par)
{
  mrc_params_parse(par, mrctest_domain_params_descr, "mrctest_domain", MPI_COMM_WORLD);
  mrc_params_print(par, mrctest_domain_params_descr, "mrctest_domain", MPI_COMM_WORLD);
}

struct mrc_domain *
mrctest_create_domain(MPI_Comm comm, struct mrctest_domain_params *par)
{
  struct mrc_domain_simple_params simple_par = {
    .ldims    = { par->gdims[0] / par->nproc[0],
		  par->gdims[1] / par->nproc[1],
		  par->gdims[2] / par->nproc[2] },
    .nr_procs = { par->nproc[0], par->nproc[1], par->nproc[2] },
  };
  struct mrc_domain *domain = mrc_domain_create(comm);
  mrc_domain_set_type(domain, "simple");
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_domain_simple_set_params(domain, &simple_par);
  mrc_crds_set_param_int(crds, "sw", SW_2);
  mrc_crds_set_param_float(crds, "xl", -30.);
  mrc_crds_set_param_float(crds, "yl", -20.);
  mrc_crds_set_param_float(crds, "zl", -20.);
  mrc_crds_set_param_float(crds, "xh",  50.);
  mrc_crds_set_param_float(crds, "yh",  20.);
  mrc_crds_set_param_float(crds, "zh",  20.);
  mrc_domain_set_from_options(domain);
  mrc_domain_view(domain);
  mrc_domain_setup(domain);

  return domain;
}

struct mrc_domain *
mrctest_create_domain_rectilinear(MPI_Comm comm, struct mrctest_domain_params *par)
{
  struct mrc_domain_simple_params simple_par = {
    .ldims    = { par->gdims[0] / par->nproc[0],
		  par->gdims[1] / par->nproc[1],
		  par->gdims[2] / par->nproc[2] },
    .nr_procs = { par->nproc[0], par->nproc[1], par->nproc[2] },
  };
  struct mrc_domain *domain = mrc_domain_create(comm);
  mrc_domain_simple_set_params(domain, &simple_par);
  mrc_domain_set_type(domain, "simple");
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "rectilinear");
  mrc_crds_set_param_int(crds, "sw", SW_2);
  mrc_crds_set_param_float(crds, "xl", -30.);
  mrc_crds_set_param_float(crds, "yl", -20.);
  mrc_crds_set_param_float(crds, "zl", -20.);
  mrc_crds_set_param_float(crds, "xh",  50.);
  mrc_crds_set_param_float(crds, "yh",  20.);
  mrc_crds_set_param_float(crds, "zh",  20.);
  mrc_domain_set_from_options(domain);
  mrc_domain_view(domain);
  mrc_domain_setup(domain);
  int sw, ldims[3];
  mrc_crds_get_param_int(crds, "sw", &sw);
  mrc_domain_get_local_offset_dims(domain, NULL, ldims);
  for (int ix = 0; ix < ldims[0] + 2 * sw; ix++) {
    MRC_CRDX(crds, ix) = ix*ix;
  }

  return domain;
}

void
mrctest_domain_init_values_0(struct mrc_f3 *f)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->domain);

  mrc_f3_foreach(f, ix,iy,iz, 2, 2) {
    float xx = MRC_CRDX(crds, ix);

    MRC_F3(f,0, ix,iy,iz) = 2.f + .2f * sin(xx);
  } mrc_f3_foreach_end;
}

static void
mrctest_domain_init_values_1(struct mrc_f3 *f)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->domain);

  mrc_f3_foreach(f, ix,iy,iz, 2, 2) {
    float yy = MRC_CRDY(crds, iy);

    MRC_F3(f,1, ix,iy,iz) = 2.f + .2f * sin(yy);
  } mrc_f3_foreach_end;
}

struct mrc_f3 *
mrctest_create_field_1(struct mrc_domain *domain)
{
  struct mrc_f3 *f3 = mrc_domain_f3_create(domain, SW_2);
  f3->name[0] = strdup("test");
  mrc_f3_setup(f3);
  mrctest_domain_init_values_0(f3);
  return f3;
}

struct mrc_f3 *
mrctest_create_field_2(struct mrc_domain *domain)
{
  struct mrc_f3 *f3 = mrc_domain_f3_create(domain, SW_2);
  mrc_f3_set_nr_comps(f3, 2);
  f3->name[0] = strdup("test0");
  f3->name[1] = strdup("test1");
  mrc_f3_setup(f3);
  mrctest_domain_init_values_0(f3);
  mrctest_domain_init_values_1(f3);
  return f3;
}

static void
mod_diagsrv(struct mrc_mod *mod, void *arg)
{
  int nr_procs_domain = mrc_mod_get_nr_procs(mod, "domain");
  mrc_io_server("xdmf_serial", "cache", nr_procs_domain);
}

void
mrctest_domain(void (*mod_domain)(struct mrc_mod *mod, void *arg))
{
  struct mrctest_domain_params par;
  mrctest_domain_init(&par);

  int nproc_domain = par.nproc[0] * par.nproc[1] * par.nproc[2];

  struct mrc_mod *mod = mrc_mod_create(MPI_COMM_WORLD);
  mrc_mod_register(mod, "domain", nproc_domain, mod_domain, &par);
  if (par.use_diagsrv) {
    mrc_mod_register(mod, "diagsrv", 1, mod_diagsrv, &par);
  }
  mrc_mod_view(mod);
  mrc_mod_setup(mod);
  mrc_mod_run(mod);
  mrc_mod_destroy(mod);
}

// ----------------------------------------------------------------------
// mrctest_f3_compare

void
mrctest_f3_compare(struct mrc_f3 *f1, struct mrc_f3 *f2, float eps)
{
  int sw = f1->sw;

  assert(f1->nr_comp == f2->nr_comp);
  for (int m = 0; m < f1->nr_comp; m++) {
    float diff = 0.;
    mrc_f3_foreach(f1, ix,iy,iz, sw, sw) {
      diff = fmaxf(diff, fabsf(MRC_F3(f1,m, ix,iy,iz) - MRC_F3(f2,m, ix,iy,iz)));
    } mrc_f3_foreach_end;
    if (diff > eps) {
      mprintf("mrctest_f3_compare: m = %d diff = %g\n", m, diff);
      assert(0);
    }
  }
}

// ----------------------------------------------------------------------
// mrctest_crds_compare

void
mrctest_crds_compare(struct mrc_crds *crds1, struct mrc_crds *crds2)
{
  int sw = crds1->par.sw;

  assert(crds1->par.sw == crds2->par.sw);
  for (int d = 0; d < 3; d++) {
    assert(crds1->par.xl[d] == crds2->par.xl[d]);
    assert(crds1->par.xh[d] == crds2->par.xh[d]);
  }

  int ldims[3];
  mrc_domain_get_local_offset_dims(crds1->domain, NULL, ldims);
  float diff = 0.;
  for (int ix = 0; ix < ldims[0] + 2 * sw; ix++) {
    diff = fmaxf(diff, fabsf(MRC_CRDX(crds1, ix) - MRC_CRDX(crds2, ix)));
    if (diff > 0.) {
      mprintf("mrctest_crds_compare: ix = %d diff = %g\n", ix, diff);
      assert(0);
    }
  }
}

