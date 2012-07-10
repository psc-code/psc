
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_collision.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

// ----------------------------------------------------------------------
// psc_test_fdtd_create

static void
psc_test_fdtd_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 100;
  psc->prm.cfl = 1.;

  psc->domain.length[0] = 1.;
  psc->domain.length[1] = 1.;
  psc->domain.length[2] = 1.;

  psc->domain.gdims[0] = 8;
  psc->domain.gdims[1] = 8;
  psc->domain.gdims[2] = 1;

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
}

// ----------------------------------------------------------------------
// psc_test_fdtd_init_field

static double
psc_test_fdtd_init_field(struct psc *psc, double x[3], int m)
{
  double kx = 2. * M_PI, ky = 2. * M_PI;

  switch (m) {
  case EX: return   1./sqrtf(2.) * sin(kx * x[0] + ky * x[1]);
  case EY: return - 1./sqrtf(2.) * sin(kx * x[0] + ky * x[1]);
  case HZ: return sin(kx * x[0] + ky * x[1]);
  default: return 0.;
  }
}

// ======================================================================
// psc_test_fdtd_ops

struct psc_ops psc_test_fdtd_ops = {
  .name             = "test_fdtd",
  .create           = psc_test_fdtd_create,
  .init_field       = psc_test_fdtd_init_field,
};

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_fdtd_ops);
}
