
#include <mrc_ts.h>
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrc_nr.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================

// FIXME
#define sqr(a) ((a) * (a))

// ======================================================================

#define BND 1

enum {
  OM_I,
  PSI_R,
  BZ_I,
  VZ_R,
  NR_FLDS,
};

struct rmhd {
  struct mrc_obj obj;

  // parameters
  float Lx;
  float lambda;
  float S;
  float d_i;
  float ky;
  float cfl;
  float dx0;
  float xn;
  float xm;
  int mx;

  struct mrc_domain *domain;
  struct mrc_f1 *By0;
};

MRC_CLASS_DECLARE(rmhd, struct rmhd);

// ======================================================================

// FIXME BND
#define CRDX(ix) (MRC_CRDX(crds, (ix)+BND))

static void
_rmhd_create(struct rmhd *rmhd)
{
  rmhd->domain = mrc_domain_create(rmhd_comm(rmhd));
  mrc_domain_set_param_int3(rmhd->domain, "m", (int [3]) { 100, 1, 1 });
}

static void
_rmhd_set_from_options(struct rmhd *rmhd)
{
  mrc_domain_set_from_options(rmhd->domain);
}

static struct mrc_f1 *
rmhd_get_fld(struct rmhd *rmhd, int nr_comps, const char *name)
{
  struct mrc_f1 *x = mrc_domain_f1_create(rmhd->domain, BND);
  mrc_f1_set_name(x, name);
  mrc_f1_set_param_int(x, "nr_comps", nr_comps);
  mrc_f1_setup(x);
  return x;
}

static void
rmhd_setup_crd(struct rmhd *rmhd)
{
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);

  float xm = rmhd->xm, xn = rmhd->xn, dx0 = rmhd->dx0;
  int mx = rmhd->mx;
  float xc = mx / 2 - .5;
  float a = (pow((rmhd->Lx / 2.) / (dx0 * xc), 1./xm) - 1.) / pow(xc, 2.*xn);
  mrc_f1_foreach(crds->crd[0], ix, 1, 1) {
    float xi = ix - xc;
    float s = 1 + a*(pow(xi, (2. * xn)));
    float sm = pow(s, xm);
    float g = dx0 * xi * sm;
    //    float dg = rmhd->dx0 * (sm + rmhd->xm*xi*2.*rmhd->xn*a*(pow(xi, (2.*rmhd->xn-1.))) * sm / s);
    CRDX(ix) = g;
  } mrc_f1_foreach_end;
}

static void
_rmhd_setup(struct rmhd *rmhd)
{
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  mrc_crds_set_type(crds, "rectilinear");
  mrc_crds_set_param_float3(crds, "l", (float [3]) { -rmhd->Lx / 2. });
  mrc_crds_set_param_float3(crds, "h", (float [3]) {  rmhd->Lx / 2. });
  mrc_crds_set_param_int(crds, "sw", 1);
  mrc_domain_setup(rmhd->domain);

  rmhd_setup_crd(rmhd);

  rmhd->By0 = rmhd_get_fld(rmhd, 1, "By0");
}

static void
_rmhd_destroy(struct rmhd *rmhd)
{
  mrc_f1_destroy(rmhd->By0);
}

#define Dxx(x, m_x, ix)							\
  (((MRC_F1(x, m_x, ix+1) - MRC_F1(x, m_x, ix  )) / (CRDX(ix+1) - CRDX(ix  )) - \
    (MRC_F1(x, m_x, ix  ) - MRC_F1(x, m_x, ix-1)) / (CRDX(ix  ) - CRDX(ix-1))) \
   / (.5 * (CRDX(ix+1) - CRDX(ix-1))))

#define Lapl(x, m_x, ix)					\
  (Dxx(x, m_x, ix) - sqr(rmhd->ky) * MRC_F1(x, m_x, ix))

static void
solve_poisson(struct rmhd *rmhd, struct mrc_f1 *x, int m_x,
	      struct mrc_f1 *b, int m_b)
{
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  static float *aa, *bb, *cc;

  if (!aa) {
    aa = calloc(rmhd->mx, sizeof(*aa));
    bb = calloc(rmhd->mx, sizeof(*bb));
    cc = calloc(rmhd->mx, sizeof(*cc));

    for (int ix = 0; ix < rmhd->mx; ix++) {
      aa[ix] =  2. / (CRDX(ix+1) - CRDX(ix-1)) * 1. / (CRDX(ix) - CRDX(ix-1));
      bb[ix] = -2. / ((CRDX(ix+1) - CRDX(ix)) * (CRDX(ix) - CRDX(ix-1))) - sqr(rmhd->ky);
      cc[ix] =  2. / (CRDX(ix+1) - CRDX(ix-1)) * 1. / (CRDX(ix+1) - CRDX(ix));
    }
  }
  
  mrc_nr_tridag(aa, bb, cc, &MRC_F1(b, m_b, 0), &MRC_F1(x, m_x, 0), rmhd->mx);
}

static void
rmhd_diag(void *ctx, float time, struct mrc_f1 *x, FILE *file)
{
  float absmax[NR_FLDS] = {};
  mrc_f1_foreach(x, ix, 0, 0) {
    for (int m = 0; m < NR_FLDS; m++) {
      absmax[m] = fmaxf(absmax[m], fabsf(MRC_F1(x,m, ix)));
    }
  } mrc_f1_foreach_end;
  
  fprintf(file, "%g", time);
  for (int m = 0; m < NR_FLDS; m++) {
    fprintf(file, " %g", absmax[m]);
  }
  fprintf(file, "\n");
  fflush(file);
}

static void
set_bnd_zero(struct rmhd *rmhd, struct mrc_f1 *x, int m_x)
{
  int mx = rmhd->mx;
  MRC_F1(x, m_x , -1) = 0.;
  MRC_F1(x, m_x , mx) = 0.;
}

#define By0(ix) MRC_F1(By0, 0, ix)

enum {
  J_R,
};

enum {
  PHI_I,
};

#define J_R(ix) MRC_F1(j, J_R, ix)
#define PHI_I(ix) MRC_F1(phi, PHI_I, ix)

#define OM_I(ix) MRC_F1(x, OM_I, ix)
#define PSI_R(ix) MRC_F1(x, PSI_R, ix)
#define BZ_I(ix) MRC_F1(x, BZ_I, ix)
#define VZ_R(ix) MRC_F1(x, VZ_R, ix)

static void
calc_rhs(void *ctx, struct mrc_f1 *rhs, struct mrc_f1 *x)
{
  struct rmhd *rmhd = ctx;
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  struct mrc_f1 *By0 = rmhd->By0;
  struct mrc_f1 *phi = rmhd_get_fld(rmhd, 1, "phi");

  set_bnd_zero(rmhd, x, OM_I);
  set_bnd_zero(rmhd, x, PSI_R);

  solve_poisson(rmhd, phi, PHI_I, x, OM_I);

  mrc_f1_foreach(x, ix, 0, 0) {
    float By0pp = 
      ((MRC_F1(By0,0, ix+1) - MRC_F1(By0,0, ix  )) / (CRDX(ix+1) - CRDX(ix  )) -
       (MRC_F1(By0,0, ix  ) - MRC_F1(By0,0, ix-1)) / (CRDX(ix  ) - CRDX(ix-1)))
      / (.5 * (CRDX(ix+1) - CRDX(ix-1)));
    float J_r = Lapl(x, PSI_R, ix);

    MRC_F1(rhs, VZ_R, ix) =
      - rmhd->ky * By0(ix) * BZ_I(ix);

    MRC_F1(rhs, OM_I, ix) =
      rmhd->ky * (By0(ix) * J_r - By0pp * PSI_R(ix));

    MRC_F1(rhs, PSI_R, ix) =
      1. / rmhd->S * J_r -
      rmhd->ky * By0(ix) * PHI_I(ix) +
      rmhd->d_i * By0(ix) * BZ_I(ix);

    MRC_F1(rhs, BZ_I, ix) =
      1. / rmhd->S * Lapl(x, BZ_I, ix) +
      rmhd->ky * By0(ix) * VZ_R(ix) +
      rmhd->d_i * rmhd->ky * (By0(ix) * J_r - By0pp * PSI_R(ix));
  }

  mrc_f1_destroy(phi);
}

// ======================================================================

#define VAR(x) (void *)offsetof(struct rmhd, x)
static struct param rmhd_param_descr[] = {
  { "Lx"              , VAR(Lx)             , PARAM_FLOAT(40.)      },
  { "lambda"          , VAR(lambda)         , PARAM_FLOAT(1.)       },
  { "S"               , VAR(S)              , PARAM_FLOAT(1000.)    },
  { "d_i"             , VAR(d_i)            , PARAM_FLOAT(0.)       },
  { "ky"              , VAR(ky)             , PARAM_FLOAT(.5)       },
  { "cfl"             , VAR(cfl)            , PARAM_FLOAT(.5)       },
  { "dx0"             , VAR(dx0)            , PARAM_FLOAT(.1)       },
  { "xn"              , VAR(xn)             , PARAM_FLOAT(2.)       },
  { "xm"              , VAR(xm)             , PARAM_FLOAT(.5)       },
  { "mx"              , VAR(mx)             , PARAM_INT(100)        },
  {},
};
#undef VAR

struct mrc_class_rmhd mrc_class_rmhd = {
  .name             = "rmhd",
  .size             = sizeof(struct rmhd),
  .param_descr      = rmhd_param_descr,
  .create           = _rmhd_create,
  .set_from_options = _rmhd_set_from_options,
  .setup            = _rmhd_setup,
  .destroy          = _rmhd_destroy,
};

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct rmhd *rmhd = rmhd_create(MPI_COMM_WORLD);
  rmhd_set_from_options(rmhd);
  rmhd_setup(rmhd);

  // i.c.
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  struct mrc_f1 *By0 = rmhd->By0;
  struct mrc_f1 *x = rmhd_get_fld(rmhd, NR_FLDS, "x");

  // setup initial equilibrium and perturbation
  mrc_f1_foreach(x, ix, 1, 1) {
    MRC_F1(By0, 0, ix) = tanh(rmhd->lambda * CRDX(ix));
    MRC_F1(x, PSI_R, ix) = exp(-sqr(CRDX(ix)));
  } mrc_f1_foreach_end;

  // write out equilibrium, etc.
  struct mrc_io *io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_param_string(io, "basename", "rmhd-ini");
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_f1_write(rmhd->By0, io);
  mrc_f1_write(x, io);
  mrc_io_close(io);
  mrc_io_destroy(io);

  float dx = rmhd->Lx / rmhd->mx; // FIXME
  float dt = rmhd->cfl * fminf(dx, rmhd->S * sqr(dx));

  // run time integration
  struct mrc_ts *ts = mrc_ts_create(MPI_COMM_WORLD);
  mrc_ts_set_context(ts, rmhd);
  mrc_ts_set_rhs_function(ts, calc_rhs);
  mrc_ts_set_diag_function(ts, rmhd_diag);
  mrc_ts_set_dt(ts, dt);
  mrc_ts_set_state(ts, x);
  mrc_ts_set_from_options(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_destroy(ts);

  mrc_f1_destroy(x);

  rmhd_destroy(rmhd);

  MPI_Finalize();
  return 0;
}
