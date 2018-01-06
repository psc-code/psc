
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

struct psc_shear {
  // parameters
  double mi_over_me;
  double wpe_over_wce;
  double Te;
  double Ti;

  // calculated from the above
  double B0;
  double di;
};

#define psc_shear(psc) mrc_to_subobj(psc, struct psc_shear)

#define VAR(x) (void *)offsetof(struct psc_shear, x)
static struct param psc_shear_descr[] = {
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(25.)           },
  { "wpe_over_wce"  , VAR(wpe_over_wce)    , PARAM_DOUBLE(2.)            },
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(.000625)       },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(.000625)       },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_shear_create

static void
psc_shear_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 100000;
  psc->prm.nr_populations = 2;
  psc->prm.nicell = 100;
  psc->prm.cfl = 0.98;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 400;
  psc->domain.gdims[2] = 100;

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
// psc_shear_setup
//
// the parameters are now set, calculate quantities to initialize fields,
// particles

static void
psc_shear_setup(struct psc *psc)
{
  struct psc_shear *sub = psc_shear(psc);

  double me = 1.;
  double mi = me * sub->mi_over_me;
  double B0 = sqrt(me) / (sub->wpe_over_wce);
  double Te = sub->Te;
  double Ti = sub->Ti;
  sub->di = sqrt(mi);

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 40. * sub->di;
  psc->domain.length[2] = 10. * sub->di;

  psc->domain.corner[0] = 0.; // no x-dependence
  psc->domain.corner[1] = 0.;
  psc->domain.corner[2] = 0.;

  MPI_Comm comm = psc_comm(psc);
  sub->B0 = B0;
  mpi_printf(comm, "shear: B0 = %g\n", B0);
  mpi_printf(comm, "shear: d_i = %g d_e\n", sub->di);
  double vth_e = sqrt(2 * sub->Te / me);
  mpi_printf(comm, "shear: Te %g vth_e %g\n", Te, vth_e);
  double vth_i = sqrt(2 * sub->Ti / mi);
  mpi_printf(comm, "shear: Ti %g vth_i %g\n", Ti, vth_i);
  mpi_printf(comm, "shear: lambda_De %g\n", sqrt(Te));
  double v_A = B0 / sqrt(mi), tau_A = psc->domain.length[2] / v_A;
  mpi_printf(comm, "shear: v_A %g t_A %g\n", v_A, tau_A);

  psc->kinds[KIND_ELECTRON].m = me;
  psc->kinds[KIND_ELECTRON].T = Te;
  psc->kinds[KIND_ION].m = mi;
  psc->kinds[KIND_ION].T = Ti;

  psc_setup_super(psc);
}

// ----------------------------------------------------------------------
// psc_shear_init_field

static double
psc_shear_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_shear *sub = psc_shear(psc);

  switch (m) {
  case HZ: return sub->B0;

  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_shear_init_npt

static void
psc_shear_init_npt(struct psc *psc, int kind, double x[3],
				struct psc_particle_npt *npt)
{
  npt->n = 1.;
}

// ----------------------------------------------------------------------
// psc_shear_read

static void
psc_shear_read(struct psc *psc, struct mrc_io *io)
{
  // do nothing -- but having this function is important so that
  // psc_shear_create() doesn't get called instead FIXME?
  psc_read_super(psc, io);
}

// ======================================================================
// psc_shear_ops

struct psc_ops psc_shear_ops = {
  .name             = "kh",
  .size             = sizeof(struct psc_shear),
  .param_descr      = psc_shear_descr,
  .create           = psc_shear_create,
  .read             = psc_shear_read,
  .setup            = psc_shear_setup,
  .init_field       = psc_shear_init_field,
  .init_npt         = psc_shear_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_shear_ops);
}
