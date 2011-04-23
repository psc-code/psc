
#include <mrc_fld.h>
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define sqr(a) ((a) * (a))

enum {
  OM_I,
  PSI_R,
  BZ_I,
  VZ_R,
  NR_FLDS,
};

struct rmhd_params {
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
  int n_end;
  int out_every;
  int diag_every;

  struct mrc_f1 *By0;
};

#define VAR(x) (void *)offsetof(struct rmhd_params, x)
static struct param rmhd_params_descr[] = {
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
  { "n_end"           , VAR(n_end)          , PARAM_INT(10000)      },
  { "out_every"       , VAR(out_every)      , PARAM_INT(1000)       },
  { "diag_every"      , VAR(diag_every)     , PARAM_INT(10)         },
  {},
};
#undef VAR

// FIXME
static struct mrc_f1 *crd;
#define CRDX(ix) (MRC_F1(crd,0, (ix)))

// FIXME
static void
mrc_f1_dump(struct mrc_f1 *f1, const char *pfx)
{
  char fname[strlen(pfx) + 5];
  sprintf(fname, "%s.asc", pfx);

  FILE *f = fopen(fname, "w");
  mrc_f1_foreach(f1, ix, 0, 0) {
    fprintf(f, "%g", CRDX(ix));
    for (int m = 0; m < f1->nr_comp; m++) {
      fprintf(f, " %g", MRC_F1(f1, m, ix));
    }
    fprintf(f, "\n");
  } mrc_f1_foreach_end;
  fclose(f);
}

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
get_fld(struct rmhd_params *par, int nr_comps)
{
  struct mrc_f1 *x = mrc_f1_create(MPI_COMM_WORLD);
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

#define Lapl(x, m_x, ix) \
  (Dxx(x, m_x, ix) - sqr(par->ky) * MRC_F1(x, m_x, ix))

static void
solve_poisson(struct rmhd_params *par, struct mrc_f1 *x, int m_x,
	      struct mrc_f1 *b, int m_b)
{
  static float *aa, *bb, *cc;

  if (!aa) {
    aa = calloc(par->mx, sizeof(*aa));
    bb = calloc(par->mx, sizeof(*bb));
    cc = calloc(par->mx, sizeof(*cc));

    for (int ix = 0; ix < par->mx; ix++) {
      aa[ix] =  2. / (CRDX(ix+1) - CRDX(ix-1)) * 1. / (CRDX(ix) - CRDX(ix-1));
      bb[ix] = -2. / ((CRDX(ix+1) - CRDX(ix)) * (CRDX(ix) - CRDX(ix-1))) - sqr(par->ky);
      cc[ix] =  2. / (CRDX(ix+1) - CRDX(ix-1)) * 1. / (CRDX(ix+1) - CRDX(ix));
    }
  }
  
  tridag(aa, bb, cc, &MRC_F1(b, m_b, 0), &MRC_F1(x, m_x, 0), par->mx);

#if 0
  // lapl
  { int ix = 0;
    MRC_F1(b, m_b, ix) =
      bb[ix] * MRC_F1(x, m_x, ix) + 
      cc[ix] * MRC_F1(x, m_x, ix+1);
  }
  for (int ix = 1; ix < par->mx - 1; ix++) {
    MRC_F1(b, m_b, ix) =
      aa[ix] * MRC_F1(x, m_x, ix-1) +
      bb[ix] * MRC_F1(x, m_x, ix) + 
      cc[ix] * MRC_F1(x, m_x, ix+1);
  }
  { int ix = par->mx - 1;
    MRC_F1(b, m_b, ix) =
      aa[ix] * MRC_F1(x, m_x, ix-1) +
      bb[ix] * MRC_F1(x, m_x, ix);
  }
#endif
}

static void
set_bnd_zero(struct rmhd_params *par, struct mrc_f1 *x, int m_x)
{
  int mx = par->mx;
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
calc_rhs(struct rmhd_params *par, struct mrc_f1 *rhs, struct mrc_f1 *x)
{
  struct mrc_f1 *By0 = par->By0;
  struct mrc_f1 *phi = get_fld(par, 1);

  set_bnd_zero(par, x, OM_I);
  set_bnd_zero(par, x, PSI_R);

  solve_poisson(par, phi, PHI_I, x, OM_I);
  //  mrc_f1_dump(phi, "phi");

  mrc_f1_foreach(x, ix, 0, 0) {
    float By0pp = 
      ((MRC_F1(By0,0, ix+1) - MRC_F1(By0,0, ix  )) / (CRDX(ix+1) - CRDX(ix  )) -
       (MRC_F1(By0,0, ix  ) - MRC_F1(By0,0, ix-1)) / (CRDX(ix  ) - CRDX(ix-1)))
      / (.5 * (CRDX(ix+1) - CRDX(ix-1)));
    float J_r = Lapl(x, PSI_R, ix);

    MRC_F1(rhs, VZ_R, ix) =
      - par->ky * By0(ix) * BZ_I(ix);

    MRC_F1(rhs, OM_I, ix) =
      par->ky * (By0(ix) * J_r - By0pp * PSI_R(ix));

    MRC_F1(rhs, PSI_R, ix) =
      1. / par->S * J_r -
      par->ky * By0(ix) * PHI_I(ix) +
      par->d_i * By0(ix) * BZ_I(ix);

    MRC_F1(rhs, BZ_I, ix) =
      1. / par->S * Lapl(x, BZ_I, ix) +
      par->ky * By0(ix) * VZ_R(ix) +
      par->d_i * par->ky * (By0(ix) * J_r - By0pp * PSI_R(ix));
  }

  //  mrc_f1_dump(rhs, "rhs");

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

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct rmhd_params par;
  mrc_params_parse(&par, rmhd_params_descr, "rmhd", MPI_COMM_WORLD);
  mrc_params_print(&par, rmhd_params_descr, "rmhd", MPI_COMM_WORLD);

  // i.c.
  struct mrc_f1 *x = get_fld(&par, NR_FLDS);
  struct mrc_f1 *By0 = get_fld(&par, 1);
  struct mrc_f1 *spacing = get_fld(&par, 1);
  par.By0 = By0;
  crd = get_fld(&par, 1);

  float a = acoff(par.mx / 2, par.Lx / 2, par.xm, par.xn, par.dx0);
  float dx = par.Lx / par.mx;
  mrc_f1_foreach(x, ix, 1, 1) {
    float xi = ix - (par.mx / 2 - .5);
    float s = 1 + a*(pow(xi, (2.*par.xn)));
    float sm = pow(s, par.xm);
    float g = par.dx0 * xi * sm;
    float dg = par.dx0 * (sm + par.xm*xi*2.*par.xn*a*(pow(xi, (2.*par.xn-1.))) * sm / s);
    CRDX(ix) = g;
    MRC_F1(spacing,0, ix) = dg;
    //    dx = fminf(dx, dg);
    //    CRDX(ix) = (xb + ((ix) + .5 ) * dx);

    MRC_F1(By0, 0, ix) = tanh(par.lambda * CRDX(ix));
    MRC_F1(x, PSI_R, ix) = exp(-sqr(CRDX(ix)));
    //    MRC_F1(x, X_PSI, ix) = sin(2.*M_PI * (xx - par.Lx / 2.) / par.Lx);
  } mrc_f1_foreach_end;

  mrc_f1_dump(By0, "By0");
  mrc_f1_dump(x, "x");
  mrc_f1_dump(spacing, "dx");
  mrc_f1_destroy(spacing);

  float dt = par.cfl * fminf(dx, par.S * sqr(dx));

  // r.h.s
  struct mrc_f1 *rhs = get_fld(&par, NR_FLDS);
  struct mrc_f1 *xm = get_fld(&par, NR_FLDS);

  FILE *f_diag = fopen("diag.asc", "w");
  for (int n = 0; n <= par.n_end; n++) {
    if (n % par.out_every == 0) {
      char pfx[10];
      sprintf(pfx, "x-%d", n / par.out_every);
      mrc_f1_dump(x, pfx);
    }

    if (n % par.diag_every == 0) {
      float absmax[NR_FLDS] = {};
      mrc_f1_foreach(x, ix, 0, 0) {
	for (int m = 0; m < NR_FLDS; m++) {
	  absmax[m] = fmaxf(absmax[m], fabsf(MRC_F1(x,m, ix)));
	}
      } mrc_f1_foreach_end;

      fprintf(f_diag, "%g", n*dt);
      for (int m = 0; m < NR_FLDS; m++) {
	fprintf(f_diag, " %g", absmax[m]);
      }
      fprintf(f_diag, "\n");
      fflush(f_diag);
    }

    calc_rhs(&par, rhs, x);
    //  mrc_f1_dump(rhs, "rhs");
    mrc_f1_foreach(x, ix, 0, 0) {
      for (int m = 0; m < NR_FLDS; m++) {
	MRC_F1(xm,m, ix) = MRC_F1(x,m, ix) + .5 * dt * MRC_F1(rhs,m, ix);
      }
    } mrc_f1_foreach_end;

    calc_rhs(&par, rhs, xm);
    //  mrc_f1_dump(rhs, "rhs");
    mrc_f1_foreach(x, ix, 0, 0) {
      for (int m = 0; m < NR_FLDS; m++) {
	MRC_F1(x,m, ix) += dt * MRC_F1(rhs,m, ix);
      }
    } mrc_f1_foreach_end;
  }
  fclose(f_diag);

  mrc_f1_destroy(x);
  mrc_f1_destroy(rhs);

  MPI_Finalize();
  return 0;
}
