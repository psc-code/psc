
#include <mrc_ts.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <mrc_bits.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================

struct rlc {
  struct mrc_obj obj;

  float R;
  float L;
  float C;
  float omega;
  float V0;
};

MRC_CLASS_DECLARE(rlc, struct rlc);

#define Q(x) MRC_F1(x, 0, 0)
#define I(x) MRC_F1(x, 0, 1)

// ======================================================================

static struct mrc_fld *
rlc_get_fld(struct rlc *rlc, const char *name)
{
  struct mrc_fld *x = mrc_fld_create(rlc_comm(rlc));
  mrc_fld_set_param_int_array(x, "dims", 2, (int [2]) { 2, 1 });
  mrc_fld_set_name(x, name);
  mrc_fld_setup(x);
  mrc_fld_view(x);
  return x;
}

static void
rlc_diag(void *ctx, float time, struct mrc_obj *_x, FILE *file)
{
  struct mrc_fld *x = (struct mrc_fld *) _x;

  fprintf(file, "%g %g %g\n", time, Q(x), I(x));
}

static void
rlc_calc_rhs(void *ctx, struct mrc_obj *_rhs, float t, struct mrc_obj *_x)
{
  struct rlc *rlc = ctx;
  struct mrc_fld *rhs = (struct mrc_fld *) _rhs, *x = (struct mrc_fld *) _x;
  float R = rlc->R, L = rlc->L, C = rlc->C;
  float V0 = rlc->V0, omega = rlc->omega;
  
  float Q = Q(x);
  float I = I(x);

  float Q_dot = I;
  float I_dot = -(R/L) * I - 1./(L*C) * Q - (V0/L) * sin(omega * t);

  Q(rhs) = Q_dot;
  I(rhs) = I_dot;
}

// ======================================================================

#define VAR(x) (void *)offsetof(struct rlc, x)
static struct param rlc_descr[] = {
  { "R"             , VAR(R)            , PARAM_FLOAT(.5)       },
  { "L"             , VAR(L)            , PARAM_FLOAT(1.)       },
  { "C"             , VAR(C)            , PARAM_FLOAT(1.)       },
  { "omega"         , VAR(omega)        , PARAM_FLOAT(2.)       },
  { "V0"            , VAR(V0)           , PARAM_FLOAT(0.)       },
  {},
};
#undef VAR

struct mrc_class_rlc mrc_class_rlc = {
  .name             = "rlc",
  .size             = sizeof(struct rlc),
  .param_descr      = rlc_descr,
};

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct rlc *rlc = rlc_create(MPI_COMM_WORLD);
  rlc_set_from_options(rlc);
  rlc_setup(rlc);
  rlc_view(rlc);

  struct mrc_fld *x = rlc_get_fld(rlc, "x");

  // set initial condition
  Q(x) = 1.;
  I(x) = 0.;

  // run time integration
  struct mrc_ts *ts = mrc_ts_create_std(MPI_COMM_WORLD, rlc_diag, rlc);
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(x));
  mrc_ts_set_rhs_function(ts, rlc_calc_rhs, rlc);
  mrc_ts_set_from_options(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  //  mrc_ts_view(ts);
  mrc_ts_destroy(ts);

  mrc_fld_destroy(x);

  rlc_destroy(rlc);

  MPI_Finalize();
  return 0;
}
