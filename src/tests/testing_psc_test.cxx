
#include "psc_testing.h"
#include "psc_push_fields.h"
#include "psc_push_particles.h"
#include "psc_bnd.h"
#include "psc_sort.h"

#include <math.h>

// ======================================================================
// psc_test

// ----------------------------------------------------------------------
// psc_test_create

void
psc_test_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc->prm.nicell = 100;

  psc->domain.gdims[0] = 16;
  psc->domain.gdims[1] = 16;
  psc->domain.gdims[2] = 16;

  psc->domain.length[0] = 1.;
  psc->domain.length[1] = 1.;
  psc->domain.length[2] = 1.;

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
// psc_test_step

void
psc_test_step(struct psc *psc)
{
  psc_output(psc);

  // field propagation n*dt -> (n+0.5)*dt
  psc_push_fields_step_a(psc->push_fields, psc->flds);

  // particle propagation n*dt -> (n+1.0)*dt
  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);

  psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
  psc_bnd_fill_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);

  // field propagation (n+0.5)*dt -> (n+1.0)*dt
  psc_push_fields_push_H(psc->push_fields, psc->flds, .5f);
  psc_push_fields_step_b2(psc->push_fields, psc->flds);
}

// ----------------------------------------------------------------------
// psc_test_init_field_linear

double
psc_test_init_field_linear(struct psc *psc, double x[3], int m)
{
  switch (m) {
  case EX: return x[0] + x[1] + x[2];
  case EY: return x[0] + x[1] + x[2];
  case EZ: return x[0] + x[1] + x[2];
  case HX: return x[0] + x[1] + x[2];
  case HY: return x[0] + x[1] + x[2];
  case HZ: return x[0] + x[1] + x[2];
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_init_npt

void
psc_test_init_npt_rest(struct psc *psc, int kind, double x[3],
		       struct psc_particle_npt *npt)
{
  npt->n = 1.;
  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 100;
    break;
  default:
    assert(0);
  }
}

// ======================================================================
// psc_test_ops_1

struct psc_ops psc_test_ops_1 = {
  .name             = "test",
  .size             = sizeof(struct psc_test),
  .create           = psc_test_create,
  .init_field       = psc_test_init_field_linear,
  .init_npt         = psc_test_init_npt_rest,
  .step             = psc_test_step,
};

