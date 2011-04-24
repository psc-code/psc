
#include <mrc_ts.h>
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================

#define sqr(a) ((a) * (a))

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

  struct mrc_f1 *crd;
  struct mrc_f1 *By0;
};

MRC_CLASS_DECLARE(rmhd, struct rmhd);

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
  .name         = "rmhd",
  .size         = sizeof(struct rmhd),
  .param_descr  = rmhd_param_descr,
};

// ======================================================================

#define CRDX(ix) (MRC_F1(crd,0, (ix)))

// FIXME
/* Taken from NR */
void
tridag(float *a, float *b, float *c, float *r, float *u, int n)
{
  int j;
  float *gam, bet;
  
  gam = calloc(sizeof(*gam), n);
  assert(fabs(b[0]) >= 1e-7);
  bet = b[0];
  u[0] = r[0] / bet;
  for (j = 1; j < n; j++) {
    gam[j] = c[j-1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (j < n-1 || fabs(bet) >= 1e-7) {
      assert(fabs(bet) >= 1e-7);
      u[j] = (r[j] - a[j] * u[j-1]) / bet;
    } else {
      assert(fabs(r[j] - a[j] * u[j-1]) < 1e-7);
      u[j] = 0.;
    }
  }
    
  for (j = n - 2; j >= 0; j--) {
    u[j] = u[j] - gam[j+1] * u[j+1];
  }
  free(gam);
}

// ======================================================================

static struct mrc_f1 *
rmhd_get_fld(struct rmhd *par, int nr_comps, const char *name)
{
  struct mrc_f1 *x = mrc_f1_create(MPI_COMM_WORLD);
  mrc_f1_set_name(x, name);
  mrc_f1_set_param_int(x, "ibx", -1);
  mrc_f1_set_param_int(x, "imx", par->mx + 2);
  mrc_f1_set_param_int(x, "nr_comps", nr_comps);
  mrc_f1_setup(x);
  return x;
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
  static float *aa, *bb, *cc;
  struct mrc_f1 *crd = rmhd->crd;

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
  
  tridag(aa, bb, cc, &MRC_F1(b, m_b, 0), &MRC_F1(x, m_x, 0), rmhd->mx);

#if 0
  // lapl
  { int ix = 0;
    MRC_F1(b, m_b, ix) =
      bb[ix] * MRC_F1(x, m_x, ix) + 
      cc[ix] * MRC_F1(x, m_x, ix+1);
  }
  for (int ix = 1; ix < rmhd->mx - 1; ix++) {
    MRC_F1(b, m_b, ix) =
      aa[ix] * MRC_F1(x, m_x, ix-1) +
      bb[ix] * MRC_F1(x, m_x, ix) + 
      cc[ix] * MRC_F1(x, m_x, ix+1);
  }
  { int ix = rmhd->mx - 1;
    MRC_F1(b, m_b, ix) =
      aa[ix] * MRC_F1(x, m_x, ix-1) +
      bb[ix] * MRC_F1(x, m_x, ix);
  }
#endif
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
  struct mrc_f1 *crd = rmhd->crd;
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

static float
acoff(int n, float y, float xm, float xn, float d0)
{
  float x = n - .5;
  float yy = y;
  yy /= d0 * x;
  yy = pow(yy, 1./xm);
  yy -= 1.;
  yy /= pow(x, 2.*xn);
  return yy;
}

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
  struct mrc_f1 *x = rmhd_get_fld(rmhd, NR_FLDS, "x");
  struct mrc_f1 *By0 = rmhd_get_fld(rmhd, 1, "By0");
  rmhd->By0 = By0;
  rmhd->crd = rmhd_get_fld(rmhd, 1, "crd");
  struct mrc_f1 *crd = rmhd->crd;

  float a = acoff(rmhd->mx / 2, rmhd->Lx / 2, rmhd->xm, rmhd->xn, rmhd->dx0);
  float dx = rmhd->Lx / rmhd->mx;
  mrc_f1_foreach(x, ix, 1, 1) {
    float xi = ix - (rmhd->mx / 2 - .5);
    float s = 1 + a*(pow(xi, (2. * rmhd->xn)));
    float sm = pow(s, rmhd->xm);
    float g = rmhd->dx0 * xi * sm;
    //    float dg = rmhd->dx0 * (sm + rmhd->xm*xi*2.*rmhd->xn*a*(pow(xi, (2.*rmhd->xn-1.))) * sm / s);
    CRDX(ix) = g;
  } mrc_f1_foreach_end;

  mrc_f1_foreach(x, ix, 1, 1) {
    MRC_F1(By0, 0, ix) = tanh(rmhd->lambda * CRDX(ix));
    MRC_F1(x, PSI_R, ix) = exp(-sqr(CRDX(ix)));
  } mrc_f1_foreach_end;

  // write out equilibrium, etc.
  struct mrc_io *io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_param_string(io, "basename", "rmhd-ini");
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_f1_write(By0, io);
  mrc_f1_write(x, io);
  mrc_io_close(io);
  mrc_io_destroy(io);

  float dt = rmhd->cfl * fminf(dx, rmhd->S * sqr(dx));

  // r.h.s
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
  mrc_f1_destroy(By0);
  mrc_f1_destroy(crd);

  rmhd_destroy(rmhd);

  MPI_Finalize();
  return 0;
}
