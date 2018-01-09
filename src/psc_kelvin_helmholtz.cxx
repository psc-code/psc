
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_single.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

// ======================================================================

enum {
  KH_ELECTRON1,
  KH_ELECTRON2,
  KH_ION1,
  KH_ION2,
  NR_KH_KINDS,
};

// ======================================================================

struct psc_kh {
  // parameters
  double beta;
  double theta_B, theta_V;
  double delta;
  double mi_over_me;
  double wpe_over_wce;
  double Ti_over_Te;
  double pert;
  double pert_vpic;
  double k_vpic;

  // calculated from the above
  double B0;
  double v0z;
};

#define to_psc_kh(psc) mrc_to_subobj(psc, struct psc_kh)

#define VAR(x) (void *)offsetof(struct psc_kh, x)
static struct param psc_kh_descr[] = {
  { "theta_B"       , VAR(theta_B)         , PARAM_DOUBLE(M_PI/2. - .05) },
  { "theta_V"       , VAR(theta_V)         , PARAM_DOUBLE(M_PI/2. - .05) },
  { "delta"         , VAR(delta)           , PARAM_DOUBLE(2.)            },
  { "beta"          , VAR(beta)            , PARAM_DOUBLE(.5)            },
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(5.)            },
  { "wpe_over_wce"  , VAR(wpe_over_wce)    , PARAM_DOUBLE(2.)            },
  { "Ti_over_Te"    , VAR(Ti_over_Te)      , PARAM_DOUBLE(1.)            },
  { "pert"          , VAR(pert)            , PARAM_DOUBLE(.0)            },
  { "pert_vpic"     , VAR(pert_vpic)       , PARAM_DOUBLE(.0)            },
  { "k_vpic"        , VAR(k_vpic)          , PARAM_DOUBLE(.5)            },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_kh_create

static void
psc_kh_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
  psc->prm.nicell = 50;
  psc->prm.cfl = 0.98;
  psc->prm.initial_momentum_gamma_correction = true;

  struct psc_kind kinds[NR_KH_KINDS] = {
    [KH_ELECTRON1] = { .name = strdup("e1"), .q = -1., .m = 1, },
    [KH_ELECTRON2] = { .name = strdup("e2"), .q = -1., .m = 1, },
    [KH_ION1]      = { .name = strdup("i1"), .q =  1.,         },
    [KH_ION2]      = { .name = strdup("i2"), .q =  1.,         },
  };
  psc_set_kinds(psc, NR_KH_KINDS, kinds);

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

  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part_hi[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc->push_particles, "1vb");
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

  double me = 1.;
  double mi = me * kh->mi_over_me;
  double B0 = sqrt(me) / (kh->wpe_over_wce);
  double vAe = B0 / sqrt(me);
  double vAi = B0 / sqrt(mi);
  double vAe_plane = vAe * cos(kh->theta_V);
  double vAi_plane = vAi * cos(kh->theta_V);
  double v0 = .5 * vAi;
  double v0z = v0 * cos(kh->theta_V - kh->theta_B);
  double Te = kh->beta * (1. / (1. + kh->Ti_over_Te)) * sqr(B0) / 2.;
  double Ti = kh->beta * (1. / (1. + 1./kh->Ti_over_Te)) * sqr(B0) / 2.;
  mpi_printf(MPI_COMM_WORLD, "psc/kh: v0=%g v0z=%g\n", v0, v0z);
  mpi_printf(MPI_COMM_WORLD, "psc/kh: vAe=%g vAe_plane=%g\n", vAe, vAe_plane);
  mpi_printf(MPI_COMM_WORLD, "psc/kh: vAi=%g vAi_plane=%g\n", vAi, vAi_plane);
  mpi_printf(MPI_COMM_WORLD, "psc/kh: Te %g Ti %g\n", Te, Ti);
  mpi_printf(MPI_COMM_WORLD, "psc/kh: lambda_De %g\n", sqrt(Te));

  kh->B0 = B0;
  kh->v0z = v0z;

  // set particle kind parameters
  assert(psc->nr_kinds == NR_KH_KINDS);
  psc->kinds[KH_ELECTRON1].T = Te;
  psc->kinds[KH_ELECTRON2].T = Te;
  psc->kinds[KH_ION1].m = mi;
  psc->kinds[KH_ION2].m = mi;
  psc->kinds[KH_ION1].T = Ti;
  psc->kinds[KH_ION2].T = Ti;

  psc_setup_super(psc);
}

// ----------------------------------------------------------------------
// vz_profile

static inline double
vz_profile(struct psc *psc, double y, double z)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double yl = psc->domain.length[1], zl = psc->domain.length[2];
  double vz = kh->v0z * tanh((y - .5 * yl * (1. + kh->pert * sin(2*M_PI * z / zl))) / kh->delta);
  vz += kh->pert_vpic * kh->v0z * sin(kh->k_vpic * z / kh->delta) * exp(-sqr(y - .5 * yl)/sqr(kh->delta));
  return vz;
}

// ----------------------------------------------------------------------
// psc_kh_init_field

static double
psc_kh_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double vz = vz_profile(psc, x[1], x[2]);

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

  double vz = vz_profile(psc, x[1], x[2]);
  npt->p[2] = vz;

  double yl = psc->domain.length[1];
  double B0x = kh->B0 * sin(kh->theta_B);
  double n = 1.;
  if (kind == KH_ELECTRON1 || kind == KH_ELECTRON2) {
    n += B0x * kh->v0z / (kh->delta * sqr(cosh((x[1] - .5 * yl) / kh->delta)));
  }
  switch (kind) {
  case KH_ELECTRON1:
  case KH_ION1:
    if (vz < 0.) {
      npt->n = 0.;
    } else {
      npt->n = n;
    }
    break;
  case KH_ELECTRON2:
  case KH_ION2:
    if (vz < 0.) {
      npt->n = n;
    } else {
      npt->n = 0.;
    }
    break;
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_kh_read

static void
psc_kh_read(struct psc *psc, struct mrc_io *io)
{
  // do nothing -- but having this function is important so that
  // psc_kh_create() doesn't get called instead
  psc_read_super(psc, io);
}

// ======================================================================
// psc_kh_ops

struct psc_ops psc_kh_ops = {
  .name             = "kh",
  .size             = sizeof(struct psc_kh),
  .param_descr      = psc_kh_descr,
  .create           = psc_kh_create,
  .read             = psc_kh_read,
  .setup            = psc_kh_setup,
  .init_field       = psc_kh_init_field,
  .init_npt         = psc_kh_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_kh_ops);
}
