
#include <mrc_ts.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <mrc_nr.h>
#include <mrc_bits.h>
#include <mrc_crds_gen.h>

#include <mrc_ts_monitor_private.h>
#include <mrc_ts_private.h>
#include <mrc_io.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

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

  struct mrc_domain *domain;
  struct mrc_fld *By0;
};

MRC_CLASS_DECLARE(rmhd, struct rmhd);

// ======================================================================

#define CRDX(ix) (MRC_CRDX(crds, ix))

static void
_rmhd_create(struct rmhd *rmhd)
{
  mrc_domain_set_param_int3(rmhd->domain, "m", (int [3]) { 100, 1, 1 });
}

static struct mrc_fld *
rmhd_get_fld(struct rmhd *rmhd, int nr_comps, const char *name)
{
  struct mrc_fld *x = mrc_domain_f1_create(rmhd->domain);
  mrc_fld_set_name(x, name);
  mrc_fld_set_param_int(x, "nr_ghosts", BND);
  mrc_fld_set_param_int(x, "nr_comps", nr_comps);
  mrc_fld_setup(x);
  return x;
}

static void
_rmhd_setup(struct rmhd *rmhd)
{
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  mrc_crds_set_param_double3(crds, "l", (double[3]) { -rmhd->Lx / 2. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  rmhd->Lx / 2. });
  mrc_crds_set_param_int(crds, "sw", BND);

  rmhd_setup_member_objs(rmhd);

  rmhd->By0 = rmhd_get_fld(rmhd, 1, "By0");
}

static void
_rmhd_destroy(struct rmhd *rmhd)
{
  mrc_fld_destroy(rmhd->By0);
}

#define Dxx(x, m_x, ix)							\
  (((MRC_F1(x, m_x, ix+1) - MRC_F1(x, m_x, ix  )) / (CRDX(ix+1) - CRDX(ix  )) - \
    (MRC_F1(x, m_x, ix  ) - MRC_F1(x, m_x, ix-1)) / (CRDX(ix  ) - CRDX(ix-1))) \
   / (.5 * (CRDX(ix+1) - CRDX(ix-1))))

#define Lapl(x, m_x, ix)					\
  (Dxx(x, m_x, ix) - sqr(rmhd->ky) * MRC_F1(x, m_x, ix))

static void
rmhd_solve_poisson(struct rmhd *rmhd, struct mrc_fld *x, int m_x,
		   struct mrc_fld *b, int m_b)
{
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  static float *aa, *bb, *cc;

  int gdims[3];
  mrc_domain_get_global_dims(rmhd->domain, gdims);
  if (!aa) {
    aa = calloc(gdims[0], sizeof(*aa));
    bb = calloc(gdims[0], sizeof(*bb));
    cc = calloc(gdims[0], sizeof(*cc));

    for (int ix = 0; ix < gdims[0]; ix++) {
      aa[ix] =  2. / (CRDX(ix+1) - CRDX(ix-1)) * 1. / (CRDX(ix) - CRDX(ix-1));
      bb[ix] = -2. / ((CRDX(ix+1) - CRDX(ix)) * (CRDX(ix) - CRDX(ix-1))) - sqr(rmhd->ky);
      cc[ix] =  2. / (CRDX(ix+1) - CRDX(ix-1)) * 1. / (CRDX(ix+1) - CRDX(ix));
    }
  }
  
  mrc_nr_tridag(aa, bb, cc, &MRC_F1(b, m_b, 0), &MRC_F1(x, m_x, 0), gdims[0]);
}

static void
rmhd_diag(void *ctx, float time, struct mrc_obj *_x, FILE *file)
{
  struct mrc_fld *x = (struct mrc_fld *) _x;
  float absmax[NR_FLDS];
  for (int m = 0; m < NR_FLDS; m++) {
    absmax[m] = mrc_fld_norm_comp(x, m);
  }
  
  fprintf(file, "%g", time);
  for (int m = 0; m < NR_FLDS; m++) {
    fprintf(file, " %g", absmax[m]);
  }

  static float absmax_last[NR_FLDS], time_last;
  for (int m = 0; m < NR_FLDS; m++) {
    fprintf(file, " %g", (log(absmax[m]) - log(absmax_last[m])) / (time - time_last));
    absmax_last[m] = absmax[m];
  }
  time_last = time;

  fprintf(file, "\n");
  fflush(file);
}

static void
rmhd_set_bnd_zero(struct rmhd *rmhd, struct mrc_fld *x, int m_x)
{
  int mx = mrc_fld_dims(x)[0];
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
rmhd_calc_rhs(void *ctx, struct mrc_obj *_rhs, float time, struct mrc_obj *_x)
{
  struct rmhd *rmhd = ctx;
  struct mrc_fld *rhs = (struct mrc_fld *) _rhs, *x = (struct mrc_fld *) _x;
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  struct mrc_fld *By0 = rmhd->By0;
  struct mrc_fld *phi = rmhd_get_fld(rmhd, 1, "phi");

  rmhd_set_bnd_zero(rmhd, x, OM_I);
  rmhd_set_bnd_zero(rmhd, x, PSI_R);

  rmhd_solve_poisson(rmhd, phi, PHI_I, x, OM_I);
  
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
      rmhd->d_i * rmhd->ky * By0(ix) * BZ_I(ix);

    MRC_F1(rhs, BZ_I, ix) =
      1. / rmhd->S * Lapl(x, BZ_I, ix) +
      rmhd->ky * By0(ix) * VZ_R(ix) +
      rmhd->d_i * rmhd->ky * (By0(ix) * J_r - By0pp * PSI_R(ix));
  }

  mrc_fld_destroy(phi);
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

  { "domain"          , VAR(domain)         , MRC_VAR_OBJ(mrc_domain)    },

  {},
};
#undef VAR

struct mrc_class_rmhd mrc_class_rmhd = {
  .name             = "rmhd",
  .size             = sizeof(struct rmhd),
  .param_descr      = rmhd_param_descr,
  .create           = _rmhd_create,
  .setup            = _rmhd_setup,
  .destroy          = _rmhd_destroy,
};

// ======================================================================

struct mrc_ts_monitor_output_phi {
  struct mrc_io *io;
  int nr;
};

static void
mrc_ts_monitor_output_phi_create(struct mrc_ts_monitor *mon)
{
  struct mrc_ts_monitor_output_phi *out =
    mrc_to_subobj(mon, struct mrc_ts_monitor_output_phi);

  out->io = mrc_io_create(mrc_ts_monitor_comm(mon));
  mrc_io_set_param_string(out->io, "basename", "run_phi");
  mrc_ts_monitor_add_child(mon, (struct mrc_obj *) out->io);
}

static void
mrc_ts_monitor_output_phi_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts)
{
  struct mrc_ts_monitor_output_phi *out =
    mrc_to_subobj(mon, struct mrc_ts_monitor_output_phi);

  mpi_printf(mrc_ts_monitor_comm(mon), "Writing output_phi %d (time = %g)\n",
	     out->nr, ts->time);
  mrc_io_open(out->io, "w", out->nr, ts->time);

  struct rmhd *rmhd = (struct rmhd *) ts->ctx_obj;
  struct mrc_fld *x = (struct mrc_fld *) ts->x;

  struct mrc_fld *phi = rmhd_get_fld(rmhd, 1, "phi");
  mrc_fld_set_comp_name(phi, 0, "phi_i");
  rmhd_solve_poisson(rmhd, phi, PHI_I, x, OM_I);
  mrc_fld_write(phi, out->io);
  mrc_fld_destroy(phi);

  mrc_io_close(out->io);
  out->nr++;
}

static struct mrc_ts_monitor_ops mrc_ts_monitor_output_phi_ops = {
  .name             = "output_phi",
  .size             = sizeof(struct mrc_ts_monitor_output_phi),
  .create           = mrc_ts_monitor_output_phi_create,
  .run              = mrc_ts_monitor_output_phi_run,
};

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor,
			      &mrc_ts_monitor_output_phi_ops);

  struct rmhd *rmhd = rmhd_create(MPI_COMM_WORLD);
  rmhd_set_from_options(rmhd);
  rmhd_setup(rmhd);
  rmhd_view(rmhd);

  // i.c.
  struct mrc_crds *crds = mrc_domain_get_crds(rmhd->domain);
  struct mrc_fld *By0 = rmhd->By0;
  struct mrc_fld *x = rmhd_get_fld(rmhd, NR_FLDS, "x");
  mrc_fld_set_comp_name(x, OM_I, "om_i");
  mrc_fld_set_comp_name(x, PSI_R, "psi_r");
  mrc_fld_set_comp_name(x, BZ_I, "bz_i");
  mrc_fld_set_comp_name(x, VZ_R, "vz_r");

  // setup initial equilibrium and perturbation
  mrc_f1_foreach(x, ix, 1, 1) {
    float By;
    float xx = CRDX(ix);
    float lambda = rmhd->lambda;
#if 0
    By = tanh(lambda * xx);
#else
    const float x0 = .92 / lambda, alpha = 1.85;
    if (xx < -x0) {
      By = -1.;
    } else if (xx > x0) {
      By = 1.;
    } else {
      By = alpha * exp(-sqr(xx*lambda)) * sqrt(M_PI) / 2. * mrc_erfi(xx*lambda);
    }
#endif
    MRC_F1(By0, 0, ix) = By;
    MRC_F1(x, PSI_R, ix) = exp(-sqr(xx));
  } mrc_f1_foreach_end;

  // write out equilibrium
  mrc_fld_dump(rmhd->By0, "By0", 0);

  // calculate dt
  int gdims[3];
  mrc_domain_get_global_dims(rmhd->domain, gdims);
  float dx = rmhd->Lx / gdims[0]; // FIXME
  float dt = rmhd->cfl * fminf(dx, rmhd->S * sqr(dx));

  // run time integration
  struct mrc_ts *ts = mrc_ts_create_std(MPI_COMM_WORLD, rmhd_diag, rmhd);

  struct mrc_ts_monitor *mon_output_phi =
    mrc_ts_monitor_create(mrc_ts_comm(ts));
  mrc_ts_monitor_set_type(mon_output_phi, "output_phi");
  mrc_ts_monitor_set_name(mon_output_phi, "mrc_ts_output_phi");
  mrc_ts_add_monitor(ts, mon_output_phi);

  mrc_ts_set_context(ts, rmhd_to_mrc_obj(rmhd));

  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(x));
  mrc_ts_set_rhs_function(ts, rmhd_calc_rhs, rmhd);
  mrc_ts_set_from_options(ts);
  mrc_ts_set_dt(ts, dt);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_view(ts);
  mrc_ts_destroy(ts);

  mrc_fld_destroy(x);

  rmhd_destroy(rmhd);

  MPI_Finalize();
  return 0;
}
