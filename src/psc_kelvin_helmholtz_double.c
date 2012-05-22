
#include <psc.h>
#include <psc_push_particles.h>
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
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
  psc->prm.nicell = 50;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 30.;
  psc->domain.length[2] = 30.;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 160;
  psc->domain.gdims[2] = 160;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
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

  psc_setup_super(psc);
}

// ----------------------------------------------------------------------
// vz_profile

static inline double
vz_profile(struct psc *psc, double y)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double yl = psc->domain.length[1];
  double vz = kh->v0z * (-1 + tanh((y - .25 * yl) / kh->delta) - tanh((y - .75 * yl) / kh->delta));
  return vz;
}

// ----------------------------------------------------------------------
// psc_kh_init_field

static double
psc_kh_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double vz = vz_profile(psc, x[1]);

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

  double vz = vz_profile(psc, x[1]);

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
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_kh_ops);
}
