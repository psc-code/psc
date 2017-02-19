
#include "mrctest.h"

#include <mrc.h>
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
mrc_fld_init_values(struct mrc_fld *fld, struct mrc_fld_init_values_info *iv_info)
{
  struct mrc_domain *domain = fld->_domain;
  assert(domain);
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  assert(crds);

  mrc_fld_set(fld, 0.);

  for (int i = 0; iv_info->ini_flds[i].ini; i++) {
    int m = iv_info->ini_flds[i].m;
    mrc_fld_foreach(fld, ix,iy,iz, 2, 2) {
      float xx = MRC_CRDX(crds, ix), yy = MRC_CRDY(crds, iy), zz = MRC_CRDZ(crds, iz);
      MRC_F3(fld,m, ix,iy,iz) = iv_info->ini_flds[i].ini(xx, yy, zz);
    } mrc_fld_foreach_end;
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
  libmrc_finalize();
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
  struct mrc_domain *domain = mrc_domain_create(comm);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_param_int3(domain, "m", par->gdims);
  mrc_domain_set_param_int3(domain, "np", par->nproc);
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_param_int(crds, "sw", SW_2);
  mrc_crds_set_param_double3(crds, "l", (double[3]) { -30., -20., -20. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  50.,  20.,  20. });
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);

  return domain;
}

struct mrc_domain *
mrctest_create_domain_rectilinear(MPI_Comm comm, struct mrctest_domain_params *par)
{
  struct mrc_domain *domain = mrc_domain_create(comm);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_param_int3(domain, "m", par->gdims);
  mrc_domain_set_param_int3(domain, "np", par->nproc);
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "rectilinear");
  mrc_crds_set_param_int(crds, "sw", SW_2);
  mrc_crds_set_param_double3(crds, "l", (double[3]) { -30., -20., -20. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  50.,  20.,  20. });
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  int sw;
  mrc_crds_get_param_int(crds, "sw", &sw);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &nr_patches);
  assert(nr_patches == 1);
  for (int ix = -sw; ix < patches[0].ldims[0] + sw; ix++) {
    MRC_CRDX(crds, ix) = ix*ix;
  }

  return domain;
}

void
mrctest_set_crds_rectilinear_1(struct mrc_domain *domain)
{
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  int sw;
  mrc_crds_get_param_int(crds, "sw", &sw);
  struct mrc_patch *patches = mrc_domain_get_patches(domain, NULL);
  for (int d = 0; d < 3; d++) {
    mrc_m1_foreach_patch(crds->crd[d], p) {
      mrc_m1_foreach(crds->crd[d], ix, sw, sw) {
	int jx = ix + patches[p].off[d];
	MRC_M1(crds->crd[d], 0, ix, p) = jx*jx;
      } mrc_m1_foreach_end;
    }
  }
}

void
mrctest_domain_init_values_0(struct mrc_fld *f)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->_domain);

  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float xx = MRC_CRDX(crds, ix);

    MRC_F3(f,0, ix,iy,iz) = 2.f + .2f * sin(xx);
  } mrc_fld_foreach_end;
}

static void
mrctest_domain_init_values_1(struct mrc_fld *f)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->_domain);

  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float yy = MRC_CRDY(crds, iy);

    MRC_F3(f,1, ix,iy,iz) = 2.f + .2f * sin(yy);
  } mrc_fld_foreach_end;
}

struct mrc_fld *
mrctest_create_field_1(struct mrc_domain *domain)
{
  struct mrc_fld *fld = mrc_domain_fld_create(domain, SW_2, "test");
  mrc_fld_setup(fld);
  mrctest_domain_init_values_0(fld);
  return fld;
}

struct mrc_fld *
mrctest_create_field_2(struct mrc_domain *domain)
{
  struct mrc_fld *fld = mrc_domain_fld_create(domain, SW_2, "test0:test1");
  mrc_fld_setup(fld);
  mrctest_domain_init_values_0(fld);
  mrctest_domain_init_values_1(fld);
  return fld;
}

struct mrc_fld *
mrctest_create_m1_1(struct mrc_domain *domain, int dim)
{
  struct mrc_fld *m1 = mrc_domain_m1_create(domain);
  mrc_fld_set_param_int(m1, "nr_ghosts", 2);
  mrc_fld_set_param_int(m1, "dim", dim);
  mrc_fld_setup(m1);
  mrc_fld_set_comp_name(m1, 0, "test");

  mrc_m1_foreach_patch(m1, p) {
    mrc_m1_foreach(m1, ix, 2, 2) {
      MRC_M1(m1, 0, ix, p) = 1.f + ix * ix;
    } mrc_m1_foreach_end;
  }
  return m1;
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
// mrctest_fld_compare

void
mrctest_fld_compare(struct mrc_fld *fld1, struct mrc_fld *fld2, float eps)
{
  assert(mrc_fld_nr_comps(fld1) == mrc_fld_nr_comps(fld2));
  int nr_comps = mrc_fld_nr_comps(fld1);
  for (int m = 0; m < nr_comps; m++) {
    float diff = 0.;
    mrc_fld_foreach(fld1, ix,iy,iz, 0, 0) {
      diff = fmaxf(diff, fabsf(MRC_F3(fld1,m, ix,iy,iz) - MRC_F3(fld2,m, ix,iy,iz)));
    } mrc_fld_foreach_end;
    if (diff > eps) {
      mprintf("mrctest_fld_compare: m = %d diff = %g\n", m, diff);
      assert(0);
    }
  }
}

// ----------------------------------------------------------------------
// mrctest_m1_compare

void
mrctest_m1_compare(struct mrc_fld *m1_1, struct mrc_fld *m1_2, float eps)
{
  assert(mrc_fld_same_shape(m1_1, m1_2));
  int sw = m1_1->_sw.vals[0];
  for (int m = 0; m < mrc_fld_nr_comps(m1_2); m++) {
    float diff = 0.;
    mrc_m1_foreach_patch(m1_1, p) {
      mrc_m1_foreach(m1_1, ix, sw, sw) {
	diff = fmaxf(diff, fabsf(MRC_M1(m1_1, m, ix, p) - MRC_M1(m1_2, m, ix, p)));
      } mrc_m1_foreach_end;
    }
    if (diff > eps) {
      mprintf("mrctest_m1_compare: m = %d diff = %g\n", m, diff);
      assert(0);
    }
  }
}

// ----------------------------------------------------------------------
// mrctest_m3_compare

void
mrctest_m3_compare(struct mrc_fld *m3_1, struct mrc_fld *m3_2)
{
  assert(mrc_fld_same_shape(m3_1, m3_2));
  mrc_fld_foreach_patch(m3_1, p) {
    struct mrc_fld_patch *m3p_1 = mrc_fld_patch_get(m3_1, p);
    struct mrc_fld_patch *m3p_2 = mrc_fld_patch_get(m3_2, p);
    float diff = 0.;
    for (int m = 0; m < mrc_fld_nr_comps(m3_1); m++) {
      mrc_m3_foreach_bnd(m3p_1, ix,iy,iz) {
	diff = fmaxf(diff, fabsf(MRC_M3(m3p_1, 0, ix,iy,iz) - MRC_M3(m3p_2, 0, ix,iy,iz)));
	if (diff > 0.) {
	  mprintf("mrc_fld_compare: ix = %d,%d,%d m = %d diff = %g %g/%g\n", ix, iy,iz, m,
		  diff, MRC_M3(m3p_1, 0, ix,iy,iz), MRC_M3(m3p_2, 0, ix,iy,iz));
	  assert(0);
	}
      } mrc_m3_foreach_end;
    }
    mrc_fld_patch_put(m3_1);
    mrc_fld_patch_put(m3_2);
  }
}

// ----------------------------------------------------------------------
// mrctest_crds_compare

void
mrctest_crds_compare(struct mrc_crds *crds1, struct mrc_crds *crds2)
{
  int sw = crds1->sw;

  assert(crds1->sw == crds2->sw);
  const double *lo1 = mrc_crds_lo(crds1), *lo2 = mrc_crds_lo(crds2);
  const double *hi1 = mrc_crds_hi(crds1), *hi2 = mrc_crds_hi(crds2);
  for (int d = 0; d < 3; d++) {
    assert(lo1[d] == lo2[d]);
    assert(hi1[d] == hi2[d]);
  }

  assert(strcmp(mrc_crds_type(crds1), mrc_crds_type(crds2)) == 0);
  for (int d = 0; d < 3; d++) {
    mrc_m1_foreach_patch(crds1->crd[d], p) {
      float diff = 0.;
      mrc_m1_foreach(crds1->crd[d], ix, sw, sw) {
	diff = fmaxf(diff, fabsf(MRC_M1(crds1->crd[d], 0, ix, p) - MRC_M1(crds2->crd[d], 0, ix, p)));
	if (diff > 0.) {
	  mprintf("mrctest_crds_compare: ix = %d diff = %g %g/%g\n", ix, diff,
		  MRC_M1(crds1->crd[d], 0, ix, p), MRC_M1(crds2->crd[d], 0, ix, p));
	  assert(0);
	}
      } mrc_m1_foreach_end;
    }

    struct mrc_ndarray *global_crd1 = crds1->global_crd[d];
    struct mrc_ndarray *global_crd2 = crds2->global_crd[d];
    mrc_ndarray_view(global_crd1);
    mrc_ndarray_view(global_crd2);
    assert(mrc_ndarray_same_shape(global_crd1, global_crd2));

    struct mrc_ndarray_it it1, it2;
    mrc_ndarray_it_all(&it1, global_crd1);
    mrc_ndarray_it_all(&it2, global_crd2);
#define TYPE double
    double diff = 0.;
    for (; !mrc_ndarray_it_done(&it1); mrc_ndarray_it_next(&it1), mrc_ndarray_it_next(&it2)) {
      diff = fmax(diff, fabs(IT_TYPE(&it1, TYPE) - IT_TYPE(&it2, TYPE)));
      if (diff > 0.) {
	mprintf("mrctest_crds_compare: idx = %d:%d diff = %g %g/%g\n", it1.idx[0], it1.idx[1], diff,
		IT_TYPE(&it1, TYPE), IT_TYPE(&it2, TYPE));
	assert(0);
      }
    }
  }
}

// ======================================================================
// mrc_domain "amr" sample setups

void
mrctest_set_amr_domain_0(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 1, 0 });
}

void
mrctest_set_amr_domain_1(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
}

void
mrctest_set_amr_domain_2(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
}

void
mrctest_set_amr_domain_3(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 5, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 5, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
}

void
mrctest_set_amr_domain_4(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 1, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 3, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 5, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 5, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
}

