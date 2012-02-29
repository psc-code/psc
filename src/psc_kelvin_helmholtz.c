
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_moments.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>

#include <mrc_params.h>

#include <math.h>

struct psc_kh {
  // parameters
  double beta;
  double theta_B, theta_V;
  double delta;
  double mi_over_me;
  double wpe_over_wce;

  // calculated from the above
  double B0;
  double v0z;
  double T;
};

#define to_psc_kh(psc) mrc_to_subobj(psc, struct psc_kh)

#define VAR(x) (void *)offsetof(struct psc_kh, x)
static struct param psc_kh_descr[] = {
  { "theta_B"       , VAR(theta_B)         , PARAM_DOUBLE(M_PI/2. - .05) },
  { "theta_V"       , VAR(theta_V)         , PARAM_DOUBLE(M_PI/2. - .05) },
  { "delta"         , VAR(delta)           , PARAM_DOUBLE(0.8944)        }, // 2/sqrt(5)
  { "beta"          , VAR(beta)            , PARAM_DOUBLE(.5)            },
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(5.)            },
  { "wpe_over_wce"  , VAR(wpe_over_wce)    , PARAM_DOUBLE(2.)            },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_kh_create

static void
psc_kh_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.nmax = 16000;
  psc->prm.cpum = 5*24.0*60*60;
  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc->prm.nicell = 50;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 30.;
  psc->domain.length[2] = 30.;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 160;
  psc->domain.gdims[2] = 160;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_CONDUCTING_WALL;
  psc->domain.bnd_fld_hi[1] = BND_FLD_CONDUCTING_WALL;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part[2] = BND_PART_PERIODIC;

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc->push_particles, "1vb");
  psc_moments_set_type(psc->moments, "1st_cc");
}

// ----------------------------------------------------------------------
// psc_kh_setup
//
// the parameters are now set, calculate quantities to initialize fields,
// particles

static void
psc_kh_setup(struct psc *psc)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double me = 1. / kh->mi_over_me;
  double B0 = sqrt(me) / (kh->wpe_over_wce);
  double vAe = B0 / sqrt(me);
  double vAe_plane = vAe * cos(kh->theta_V);
  double v0 = 2 * vAe_plane;
  double v0z = v0 * cos(kh->theta_V - kh->theta_B);
  double T = kh->beta * sqr(B0) / 2.;

  kh->B0 = B0;
  kh->v0z = v0z;
  kh->T = T;

  psc_setup_default(psc);
}

// ----------------------------------------------------------------------
// psc_kh_init_field

static double
psc_kh_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double yl = psc->domain.length[1];
  double vz = kh->v0z * tanh((x[1] - .5 * yl) / kh->delta);

  switch (m) {
  case HX: return kh->B0 * sin(kh->theta_B);
  case HZ: return kh->B0 * cos(kh->theta_B);
  case EY: return -vz * kh->B0 * sin(kh->theta_B);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_kh_init_npt

static void
psc_kh_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double yl = psc->domain.length[1];
  double vz = kh->v0z * tanh((x[1] - .5 * yl) / kh->delta);

  npt->n = 1;
  npt->p[2] = vz;
  npt->T[0] = kh->T;
  npt->T[1] = kh->T;
  npt->T[2] = kh->T;

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1. / kh->mi_over_me;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 1.;
    break;
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_kh_output
//
// keep track of energy conservation

static void
psc_kh_output(struct psc *psc)
{
  static FILE *file;
  int rank;

  MPI_Comm_rank(psc_comm(psc), &rank);

  if (!file && rank == 0) {
    file = fopen("diag.asc", "w");
    fprintf(file, "# time E2 H2\n");
  }

  
  mfields_c_t *flds = psc_mfields_get_c(psc->flds, EX, HX + 3);
  
  double EH2[2] = {};
  psc_foreach_patch(ppsc, p) {
    fields_c_t *pf = psc_mfields_get_patch_c(flds, p);
    // FIXME, this doesn't handle non-periodic b.c. right
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      EH2[0] +=
	(sqr(F3_C(pf, EX, ix,iy,iz)) +
	 sqr(F3_C(pf, EY, ix,iy,iz)) +
	 sqr(F3_C(pf, EZ, ix,iy,iz))) * psc->dx[0] * psc->dx[1] * psc->dx[2];

      EH2[1] +=
	(sqr(F3_C(pf, HX, ix,iy,iz)) +
	 sqr(F3_C(pf, HY, ix,iy,iz)) +
	 sqr(F3_C(pf, HZ, ix,iy,iz))) * psc->dx[0] * psc->dx[1] * psc->dx[2];
    } foreach_3d_end;
  }

  MPI_Reduce(MPI_IN_PLACE, &EH2, 2, MPI_DOUBLE, MPI_SUM, 0, psc_comm(psc));
  for (int i = 0; i < 2; i++) {
    EH2[i] /= psc->domain.length[0] * psc->domain.length[1] * psc->domain.length[2];
  }
  if (rank == 0) {
    fprintf(file, "%g %g %g\n", psc->timestep * psc->dt, EH2[0], EH2[1]);
    fflush(file);
  }

  psc_mfields_put_c(flds, psc->flds, 0, 0);

  psc_output_default(psc);
}

// ======================================================================
// psc_kh_ops

struct psc_ops psc_kh_ops = {
  .name             = "kh",
  .size             = sizeof(struct psc_kh),
  .param_descr      = psc_kh_descr,
  .create           = psc_kh_create,
  .setup            = psc_kh_setup,
  .init_field       = psc_kh_init_field,
  .init_npt         = psc_kh_init_npt,
  .output           = psc_kh_output,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_kh_ops);
}
