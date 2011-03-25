
#include "psc.h"
#include "psc_case_private.h"
#include "psc_pulse.h"
#include "psc_push_fields.h"
#include "psc_bnd_fields.h"

#include <mrc_params.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ----------------------------------------------------------------------
// case microsphere
//

struct psc_case_microsphere {
  // parameters
};

#define VAR(x) (void *)offsetof(struct psc_case_microsphere, x)

static struct param psc_case_microsphere_descr[] = {
  {},
};

#undef VAR

static double
microsphere_dens(struct psc_case *_case, double x[3])
{
  // center
  const double xc[3] = { 5. * 1e-6, 5. * 1e-6, 5. * 1e-6 };
  const double radius = 2. * 1e-6;
  const double thickness = .4 * 1e-6;

  double xr[3];
  for (int d = 0; d < 3; d++) {
    xr[d] = x[d] * psc.coeff.ld - xc[0];
  };

  double r = sqrt(sqr(xr[0]) + sqr(xr[1]) + sqr(xr[2]));

  return exp(-sqr((r - radius) / thickness));
}

static void
psc_case_microsphere_set_from_options(struct psc_case *_case)
{
  //  struct psc_case_microsphere *msphere = mrc_to_subobj(_case, struct psc_case_microsphere);

  psc.prm.nicell = 10;
  psc.prm.nmax = 201;

  psc.domain.length[0] = 10 * 1e-6;
  psc.domain.length[1] = 10 * 1e-6;
  psc.domain.length[2] = 10 * 1e-6;

  psc.domain.gdims[0] = 128;
  psc.domain.gdims[1] = 128;
  psc.domain.gdims[2] = 128;

  psc.domain.bnd_fld_lo[0] = BND_FLD_OPEN;
  psc.domain.bnd_fld_hi[0] = BND_FLD_OPEN;
  psc.domain.bnd_fld_lo[1] = BND_FLD_OPEN;
  psc.domain.bnd_fld_hi[1] = BND_FLD_OPEN;
  psc.domain.bnd_fld_lo[2] = BND_FLD_OPEN;
  psc.domain.bnd_fld_hi[2] = BND_FLD_OPEN;
  psc.domain.bnd_part[0] = BND_PART_PERIODIC; // FIXME
  psc.domain.bnd_part[1] = BND_PART_PERIODIC;
  psc.domain.bnd_part[2] = BND_PART_PERIODIC;

  double *length = psc.domain.length;
  double width_normal = 2. * 1e-6;
  double width_par    = 3. * 1e-6;

  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);

  struct psc_pulse *pulse_x1 = psc_bnd_fields_get_pulse_x1(bnd_fields);
  psc_pulse_set_type(pulse_x1, "gauss");
  psc_pulse_set_param_double3(pulse_x1, "m",  (double[3]) { 0., .5 * length[1], .5 * length[2] });
  psc_pulse_set_param_double3(pulse_x1, "dm", (double[3]) { width_normal, width_par, width_par });
  psc_pulse_set_param_double(pulse_x1, "amplitude_s", 1.);

  struct psc_pulse *pulse_x2 = psc_bnd_fields_get_pulse_x2(bnd_fields);
  psc_pulse_set_type(pulse_x2, "gauss");
  psc_pulse_set_param_double3(pulse_x2, "m",  (double[3]) { length[0], .5 * length[1], .5 * length[2] });
  psc_pulse_set_param_double3(pulse_x2, "dm", (double[3]) { width_normal, width_par, width_par });
  psc_pulse_set_param_double(pulse_x2, "amplitude_s", 1.);


  struct psc_pulse *pulse_y1 = psc_bnd_fields_get_pulse_y1(bnd_fields);
  psc_pulse_set_type(pulse_y1, "gauss");
  psc_pulse_set_param_double3(pulse_y1, "m",  (double[3]) { .5 * length[0], 0., .5 * length[2] });
  psc_pulse_set_param_double3(pulse_y1, "dm", (double[3]) { width_par, width_normal, width_par });
  psc_pulse_set_param_double(pulse_y1, "amplitude_s", 1.);

  struct psc_pulse *pulse_y2 = psc_bnd_fields_get_pulse_y2(bnd_fields);
  psc_pulse_set_type(pulse_y2, "gauss");
  psc_pulse_set_param_double3(pulse_y2, "m",  (double[3]) { .5 * length[0], length[1], .5 * length[2] });
  psc_pulse_set_param_double3(pulse_y2, "dm", (double[3]) { width_par, width_normal, width_par });
  psc_pulse_set_param_double(pulse_y2, "amplitude_s", 1.);


  struct psc_pulse *pulse_z1 = psc_bnd_fields_get_pulse_z1(bnd_fields);
  psc_pulse_set_type(pulse_z1, "gauss");
  psc_pulse_set_param_double3(pulse_z1, "m",  (double[3]) { .5 * length[0], .5 * length[1], 0.});
  psc_pulse_set_param_double3(pulse_z1, "dm", (double[3]) { width_par, width_par, width_normal });
  psc_pulse_set_param_double(pulse_z1, "amplitude_s", 1.);

  struct psc_pulse *pulse_z2 = psc_bnd_fields_get_pulse_z2(bnd_fields);
  psc_pulse_set_type(pulse_z2, "gauss");
  psc_pulse_set_param_double3(pulse_z2, "m",  (double[3]) { .5 * length[0], .5 * length[1], length[2]});
  psc_pulse_set_param_double3(pulse_z2, "dm", (double[3]) { width_par, width_par, width_normal });
  psc_pulse_set_param_double(pulse_z2, "amplitude_s", 1.);
}

static void
psc_case_microsphere_init_npt(struct psc_case *_case, int kind, double x[3],
			      struct psc_particle_npt *npt)
{
  //  struct psc_case_microsphere *msphere = mrc_to_subobj(_case, struct psc_case_microsphere);
  double dens = microsphere_dens(_case, x);

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1;
    npt->n = dens;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 100;
    npt->n = dens;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_microsphere_ops = {
  .name             = "microsphere",
  .size             = sizeof(struct psc_case_microsphere),
  .param_descr      = psc_case_microsphere_descr,
  .set_from_options = psc_case_microsphere_set_from_options,
  .init_npt         = psc_case_microsphere_init_npt,
};

