
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

struct psc_island_coalescence {
  // parameters
  double mi_over_me;
  double wpe_over_wce;
  double eps;
  double lambda; // in terms of d_i
  double pert;
  double nb;
  double Ti_over_Te;

  // calculated from the above
  double B0;
  double di;
  double Te;
  double Ti;
};

#define psc_island_coalescence(psc) mrc_to_subobj(psc, struct psc_island_coalescence)

#define VAR(x) (void *)offsetof(struct psc_island_coalescence, x)
static struct param psc_island_coalescence_descr[] = {
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(25.)           },
  { "wpe_over_wce"  , VAR(wpe_over_wce)    , PARAM_DOUBLE(2.)            },
  { "eps"           , VAR(eps)             , PARAM_DOUBLE(.4)            },
  { "lambda"        , VAR(lambda)          , PARAM_DOUBLE(10.)           },
  { "pert"          , VAR(pert)            , PARAM_DOUBLE(.1)            },
  { "nb"            , VAR(nb)              , PARAM_DOUBLE(.2)            },
  { "Ti_over_Te"    , VAR(Ti_over_Te)      , PARAM_DOUBLE(1.)            },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_island_coalescence_create

static void
psc_island_coalescence_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
  psc->prm.nr_populations = 4;
  psc->prm.nicell = 100;
  psc->prm.cfl = 0.98;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 448;
  psc->domain.gdims[2] = 896;

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
// psc_island_coalescence_setup
//
// the parameters are now set, calculate quantities to initialize fields,
// particles

static void
psc_island_coalescence_setup(struct psc *psc)
{
  struct psc_island_coalescence *sub = psc_island_coalescence(psc);

  double me = 1.;
  double mi = me * sub->mi_over_me;
  double B0 = sqrt(me) / (sub->wpe_over_wce);
  double Te = sub->Te;
  double Ti = sub->Ti;
  sub->di = sqrt(mi);
  double lambda_de = sub->lambda * sub->di;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = lambda_de * 2. * M_PI;
  psc->domain.length[2] = lambda_de * 4. * M_PI;

  psc->domain.corner[0] = 0.; // no x-dependence
  psc->domain.corner[1] = -.5 * psc->domain.length[1];
  psc->domain.corner[2] = -.5 * psc->domain.length[2];

  MPI_Comm comm = psc_comm(psc);
  sub->B0 = B0;
  mpi_printf(comm, "island coalescence: B0 = %g\n", B0);
  mpi_printf(comm, "island coalescence: lambda = %g d_i = %g d_e\n", sub->lambda, lambda_de);
  mpi_printf(comm, "island coalescence: Ti / Te = %g\n", sub->Ti_over_Te);
  double T_sum = sqr(B0) / 2.;
  Te = T_sum / (1. + sub->Ti_over_Te);
  Ti = T_sum / (1. + (1./sub->Ti_over_Te));
  sub->Te = Te;
  sub->Ti = Ti;
  double vth_e = sqrt(2 * Te / me);
  mpi_printf(comm, "island coalescence: Te %g vth_e %g\n", Te, vth_e);
  double vth_i = sqrt(2 * Ti / mi);
  mpi_printf(comm, "island coalescence: Ti %g vth_i %g\n", Ti, vth_i);
  mpi_printf(comm, "island coalescence: lambda_De %g\n", sqrt(Te));
  double v_A = B0 / sqrt(mi), tau_A = psc->domain.length[2] / v_A;
  mpi_printf(comm, "island_coalescence: v_A %g t_A %g\n", v_A, tau_A);

  psc->kinds[KIND_ELECTRON].m = me;
  psc->kinds[KIND_ELECTRON].T = Te;
  psc->kinds[KIND_ION].m = mi;
  psc->kinds[KIND_ION].T = Ti;

  psc_setup_super(psc);
}

// ----------------------------------------------------------------------
// psi

static double
psi(struct psc *psc, double x[3])
{
  struct psc_island_coalescence *sub = psc_island_coalescence(psc);
  double lambda_de = sub->lambda * sub->di, eps = sub->eps, B0 = sub->B0, pert = sub->pert;
  double *length = psc->domain.length;

  double val = -lambda_de * B0 * log(cosh(x[1] / lambda_de) + eps * cos(x[2] / lambda_de));
  val += pert * B0 * length[2] / (2. * M_PI) * 
    cos(M_PI * x[1] / length[1]) * cos(2. * M_PI * x[2] / length[2]);
  return val;
}

// ----------------------------------------------------------------------
// current

static double
current(struct psc *psc, double x[3])
{
  struct psc_island_coalescence *sub = psc_island_coalescence(psc);
  double lambda_de = sub->lambda * sub->di, eps = sub->eps, B0 = sub->B0, pert = sub->pert;
  double *length = psc->domain.length;

  double val = B0 / lambda_de * 
    (1. - sqr(eps)) * pow(cosh(x[1] / lambda_de) + eps * cos(x[2] / lambda_de), -2.);
  val -= pert * B0 * length[2] / (2. * M_PI) * (sqr(M_PI / length[1]) + sqr(2. * M_PI / length[2])) 
      * cos(M_PI * x[1] / length[1]) * cos(2. * M_PI * x[2] / length[2]);
  return val;
}

// ----------------------------------------------------------------------
// psc_island_coalescence_init_field

static double
psc_island_coalescence_init_field(struct psc *psc, double x[3], int m)
{
  double dy = psc->patch[0].dx[1], dz = psc->patch[0].dx[2];

  switch (m) {
  case HY: return   (psi(psc, (double[3]) { x[0], x[1], x[2] + .5 * dz }) - 
		     psi(psc, (double[3]) { x[0], x[1], x[2] - .5 * dz })) / dz;
  case HZ: return - (psi(psc, (double[3]) { x[0], x[1] + .5 * dy, x[2] }) - 
		     psi(psc, (double[3]) { x[0], x[1] - .5 * dy, x[2] })) / dy;

  // not needed anymore, just for initial j output
  case JXI: return current(psc, (double [3]) { x[0], x[1], x[2] });

  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_island_coalescence_init_npt

static void
psc_island_coalescence_init_npt(struct psc *psc, int pop, double x[3],
				struct psc_particle_npt *npt)
{
  struct psc_island_coalescence *sub = psc_island_coalescence(psc);
  double lambda_de = sub->lambda * sub->di, eps = sub->eps;
  double nb = sub->nb, Te = sub->Te, Ti = sub->Ti;

  double n =
    (1. - sqr(eps)) * pow(cosh(x[1] / lambda_de) + eps * cos(x[2] / lambda_de), -2.);
  double j = current(psc, x);

  switch (pop) {
  case 0: // ion drifting
    npt->n = n;
    npt->q = 1.;
    npt->m = sub->mi_over_me;
    npt->p[0] =   j / n * Ti / (Te + Ti);
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    npt->kind = KIND_ION;
    break;
  case 1: // ion bg
    npt->n = nb;
    npt->q = 1.;
    npt->m = sub->mi_over_me;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    npt->kind = KIND_ION;
    break;
  case 2: // electron drifting
    npt->n = n;
    npt->q = -1.;
    npt->m = 1.;
    npt->p[0] = - j / n * Te / (Te + Ti);
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    npt->kind = KIND_ELECTRON;
    break;
  case 3: // electron bg
    npt->n = nb;
    npt->q = -1.;
    npt->m = 1.;
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    npt->kind = KIND_ELECTRON;
    break;
  default:
    assert(0);
  };
}

// ----------------------------------------------------------------------
// psc_island_coalescence_read

static void
psc_island_coalescence_read(struct psc *psc, struct mrc_io *io)
{
  // do nothing -- but having this function is important so that
  // psc_island_coalescence_create() doesn't get called instead FIXME?
  psc_read_super(psc, io);
}

// ======================================================================
// psc_island_coalescence_ops

struct psc_ops psc_island_coalescence_ops = {
  .name             = "kh",
  .size             = sizeof(struct psc_island_coalescence),
  .param_descr      = psc_island_coalescence_descr,
  .create           = psc_island_coalescence_create,
  .read             = psc_island_coalescence_read,
  .setup            = psc_island_coalescence_setup,
  .init_field       = psc_island_coalescence_init_field,
  .init_npt         = psc_island_coalescence_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_island_coalescence_ops);
}
