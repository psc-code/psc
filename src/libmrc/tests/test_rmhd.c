
#include <mrc_fld.h>
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define sqr(a) ((a) * (a))

enum {
  X_OM,
  X_PSI,
  X_NR,
};

struct rmhd_params {
  float Lx;
  float lambda;
  float S;
  float ky;

  int mx;
  int n_end;
  int out_every;
  int diag_every;

  struct mrc_f1 *By0;
};

#define VAR(x) (void *)offsetof(struct rmhd_params, x)
static struct param rmhd_params_descr[] = {
  { "Lx"              , VAR(Lx)             , PARAM_FLOAT(10.)      },
  { "lambda"          , VAR(lambda)         , PARAM_FLOAT(1.)       },
  { "S"               , VAR(S)              , PARAM_FLOAT(100.)     },
  { "ky"              , VAR(ky)             , PARAM_FLOAT(1.)       },
  { "mx"              , VAR(mx)             , PARAM_INT(100)        },
  { "n_end"           , VAR(n_end)          , PARAM_INT(100)        },
  { "out_every"       , VAR(out_every)      , PARAM_INT(10)         },
  { "diag_every"      , VAR(diag_every)     , PARAM_INT(10)         },
  {},
};
#undef VAR

// FIXME
static float xb, dx;
#define CRDX(ix) (xb + ((ix) + .5 ) * dx)

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

static void
calc_laplace(struct rmhd_params *par, struct mrc_f1 *b, int m_b,
	     struct mrc_f1 *x, int m_x)
{
  mrc_f1_foreach(x, ix, 0, 0) {
    MRC_F1(b, m_b, ix) =
      (MRC_F1(x, m_x, ix+1) - 2.*MRC_F1(x, m_x, ix) + MRC_F1(x, m_x, ix-1)) / sqr(dx) 
      - sqr(par->ky) * MRC_F1(x, m_x, ix);
  }
}

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
      aa[ix] = 1. / sqr(dx);
      bb[ix] = -2. / sqr(dx) - sqr(par->ky);
      cc[ix] = 1. / sqr(dx);
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

static void
calc_rhs(struct rmhd_params *par, struct mrc_f1 *rhs, struct mrc_f1 *x)
{
  struct mrc_f1 *By0 = par->By0;
  struct mrc_f1 *j = get_fld(par, 1);
  struct mrc_f1 *phi = get_fld(par, 1);

  set_bnd_zero(par, x, X_OM);
  set_bnd_zero(par, x, X_PSI);

  calc_laplace(par, j, 0, x, X_PSI);
  //  set_bnd_zero(par, j, 0);
  //  mrc_f1_dump(j, "j");
  solve_poisson(par, phi, 0, x, X_OM);
  //  set_bnd_zero(par, phi, 0);
  //  mrc_f1_dump(phi, "phi");

  mrc_f1_foreach(x, ix, 0, 0) {
    MRC_F1(rhs, X_OM, ix) =
      par->ky * (MRC_F1(By0,0, ix) * MRC_F1(j,0, ix) -
		 (MRC_F1(By0,0, ix+1) - 2.*MRC_F1(By0,0, ix) + MRC_F1(By0,0, ix-1)) / sqr(dx) *
		 MRC_F1(x, X_PSI, ix));
    MRC_F1(rhs, X_PSI, ix) =
      1. / par->S * MRC_F1(j,0, ix) -
      par->ky * MRC_F1(By0,0, ix) * MRC_F1(phi,0, ix);
  }

  //  mrc_f1_dump(rhs, "rhs");

  mrc_f1_destroy(j);
  mrc_f1_destroy(phi);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct rmhd_params par;
  mrc_params_parse(&par, rmhd_params_descr, "rmhd", MPI_COMM_WORLD);
  mrc_params_print(&par, rmhd_params_descr, "rmhd", MPI_COMM_WORLD);

  dx = par.Lx / par.mx;
  xb = - par.Lx / 2.;
  double dt = .5 * fminf(dx, par.S * sqr(dx));

  // i.c.
  struct mrc_f1 *x = get_fld(&par, X_NR);
  struct mrc_f1 *By0 = get_fld(&par, 1);
  par.By0 = By0;

  mrc_f1_foreach(x, ix, 1, 1) {
    float xx = CRDX(ix);
    MRC_F1(By0, 0, ix) = tanh(par.lambda * xx);
    //    MRC_F1(x, X_OM, ix) = xx * exp(-sqr(xx));
    MRC_F1(x, X_PSI, ix) = exp(-sqr(xx));
    //    MRC_F1(x, X_PSI, ix) = sin(2.*M_PI * (xx - par.Lx / 2.) / par.Lx);
  } mrc_f1_foreach_end;

  mrc_f1_dump(By0, "By0");
  mrc_f1_dump(x, "x");

  // r.h.s
  struct mrc_f1 *rhs = get_fld(&par, X_NR);
  struct mrc_f1 *xm = get_fld(&par, X_NR);

  FILE *f_diag = fopen("diag.asc", "w");
  for (int n = 0; n <= par.n_end; n++) {
    calc_rhs(&par, rhs, x);
    //  mrc_f1_dump(rhs, "rhs");
    mrc_f1_foreach(x, ix, 0, 0) {
      for (int m = 0; m < X_NR; m++) {
	MRC_F1(xm,m, ix) = MRC_F1(x,m, ix) + .5 * dt * MRC_F1(rhs,m, ix);
      }
    } mrc_f1_foreach_end;

    calc_rhs(&par, rhs, xm);
    //  mrc_f1_dump(rhs, "rhs");
    mrc_f1_foreach(x, ix, 0, 0) {
      for (int m = 0; m < X_NR; m++) {
	MRC_F1(x,m, ix) += dt * MRC_F1(rhs,m, ix);
      }
    } mrc_f1_foreach_end;

    if (n % par.out_every == 0) {
      char pfx[10];
      sprintf(pfx, "x-%d", n / par.out_every);
      mrc_f1_dump(x, pfx);
    }

    if (n % par.diag_every == 0) {
      float absmax[X_NR] = {};
      mrc_f1_foreach(x, ix, 0, 0) {
	for (int m = 0; m < X_NR; m++) {
	  absmax[m] = fmaxf(absmax[m], fabsf(MRC_F1(x,m, ix)));
	}
      } mrc_f1_foreach_end;

      fprintf(f_diag, "%g %g %g\n", n*dt, absmax[0], absmax[1]);
      fflush(f_diag);
    }
  }
  fclose(f_diag);

  mrc_f1_destroy(x);
  mrc_f1_destroy(rhs);

  MPI_Finalize();
  return 0;
}
