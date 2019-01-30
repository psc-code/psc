
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_push_particles.h>
#include <psc_bnd_fields.h>

#include <mrc_params.h>

#include <math.h>

struct psc_test_microsphere {
  // parameters
  double radius;         //< radius of spherical shell
  double thickness;      //< thickness of shell
  double xc[3];          //< center of sphere, physical coords
  double width_normal;   //< pulse width normal to propagation direction
  double width_parallel; //< pulse width along propagation direction
};

#define VAR(x) (void *)offsetof(struct psc_test_microsphere, x)
static struct param psc_test_microsphere_descr[] = {
  { "radius"        , VAR(radius)          , PARAM_DOUBLE(2.e-6)              },
  { "thickness"     , VAR(thickness)       , PARAM_DOUBLE(.4e-6)              },
  { "center"        , VAR(xc)              , PARAM_DOUBLE3(5e-6, 5e-6, 5e-6)  },
  { "width_normal"  , VAR(width_normal)    , PARAM_DOUBLE(3e-6)               },
  { "width_parallel", VAR(width_parallel)  , PARAM_DOUBLE(2e-6)               },
  {},
};
#undef VAR

static void
psc_test_microsphere_create(struct psc *psc)
{
  psc->prm.nicell = 10;
  psc->prm.nmax = 201;

  psc->domain.length[0] = 10 * 1e-6;
  psc->domain.length[1] = 10 * 1e-6;
  psc->domain.length[2] = 10 * 1e-6;

  psc->domain.gdims[0] = 128;
  psc->domain.gdims[1] = 128;
  psc->domain.gdims[2] = 128;

  psc->domain.bnd_fld_lo[0] = BND_FLD_OPEN;
  psc->domain.bnd_fld_hi[0] = BND_FLD_OPEN;
  psc->domain.bnd_fld_lo[1] = BND_FLD_OPEN;
  psc->domain.bnd_fld_hi[1] = BND_FLD_OPEN;
  psc->domain.bnd_fld_lo[2] = BND_FLD_OPEN;
  psc->domain.bnd_fld_hi[2] = BND_FLD_OPEN;

  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;
}

static double
microsphere_dens(struct psc *psc, double x[3])
{
  struct psc_test_microsphere *msphere = mrc_to_subobj(psc, struct psc_test_microsphere);

  double xr[3];
  for (int d = 0; d < 3; d++) {
    xr[d] = x[d] * ppsc->coeff.ld - msphere->xc[d];
  };

  double r = sqrt(sqr(xr[0]) + sqr(xr[1]) + sqr(xr[2]));

  return exp(-sqr((r - msphere->radius) / msphere->thickness));
}

static void
psc_test_microsphere_set_from_options(struct psc *psc)
{
#if 0
  struct psc_test_microsphere *msphere = mrc_to_subobj(psc, struct psc_test_microsphere);

  double *length  = ppsc->domain.length;
  double w_normal = msphere->width_normal;
  double w_par    = msphere->width_parallel;

  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);

  struct psc_pulse *pulse_x1 = psc_bnd_fields_get_pulse_x1(bnd_fields);
  psc_pulse_set_type(pulse_x1, "gauss");
  psc_pulse_set_param_double3(pulse_x1, "m",  (double[3]) { 0., .5 * length[1], .5 * length[2] });
  psc_pulse_set_param_double3(pulse_x1, "dm", (double[3]) { w_par, w_normal, w_normal });
  psc_pulse_set_param_double(pulse_x1, "amplitude_s", 1.);

  struct psc_pulse *pulse_x2 = psc_bnd_fields_get_pulse_x2(bnd_fields);
  psc_pulse_set_type(pulse_x2, "gauss");
  psc_pulse_set_param_double3(pulse_x2, "m",  (double[3]) { length[0], .5 * length[1], .5 * length[2] });
  psc_pulse_set_param_double3(pulse_x2, "dm", (double[3]) { w_par, w_normal, w_normal });
  psc_pulse_set_param_double(pulse_x2, "amplitude_s", 1.);


  struct psc_pulse *pulse_y1 = psc_bnd_fields_get_pulse_y1(bnd_fields);
  psc_pulse_set_type(pulse_y1, "gauss");
  psc_pulse_set_param_double3(pulse_y1, "m",  (double[3]) { .5 * length[0], 0., .5 * length[2] });
  psc_pulse_set_param_double3(pulse_y1, "dm", (double[3]) { w_normal, w_par, w_normal });
  psc_pulse_set_param_double(pulse_y1, "amplitude_s", 1.);

  struct psc_pulse *pulse_y2 = psc_bnd_fields_get_pulse_y2(bnd_fields);
  psc_pulse_set_type(pulse_y2, "gauss");
  psc_pulse_set_param_double3(pulse_y2, "m",  (double[3]) { .5 * length[0], length[1], .5 * length[2] });
  psc_pulse_set_param_double3(pulse_y2, "dm", (double[3]) { w_normal, w_par, w_normal });
  psc_pulse_set_param_double(pulse_y2, "amplitude_s", 1.);


  struct psc_pulse *pulse_z1 = psc_bnd_fields_get_pulse_z1(bnd_fields);
  psc_pulse_set_type(pulse_z1, "gauss");
  psc_pulse_set_param_double3(pulse_z1, "m",  (double[3]) { .5 * length[0], .5 * length[1], 0.});
  psc_pulse_set_param_double3(pulse_z1, "dm", (double[3]) { w_normal, w_normal, w_par });
  psc_pulse_set_param_double(pulse_z1, "amplitude_s", 1.);

  struct psc_pulse *pulse_z2 = psc_bnd_fields_get_pulse_z2(bnd_fields);
  psc_pulse_set_type(pulse_z2, "gauss");
  psc_pulse_set_param_double3(pulse_z2, "m",  (double[3]) { .5 * length[0], .5 * length[1], length[2]});
  psc_pulse_set_param_double3(pulse_z2, "dm", (double[3]) { w_normal, w_normal, w_par });
  psc_pulse_set_param_double(pulse_z2, "amplitude_s", 1.);
#endif
}

static void
psc_test_microsphere_init_npt(struct psc *psc, int kind, double x[3],
			      struct psc_particle_npt *npt)
{
  //  struct psc_test_microsphere *msphere = mrc_to_subobj(_case, struct psc_test_microsphere);
  double dens = microsphere_dens(psc, x);

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

// ======================================================================
// psc_test_microsphere_ops

struct psc_ops psc_test_microsphere_ops = {
  .name             = "microsphere",
  .size             = sizeof(struct psc_test_microsphere),
  .param_descr      = psc_test_microsphere_descr,
  .create           = psc_test_microsphere_create,
  .set_from_options = psc_test_microsphere_set_from_options,
  .init_npt         = psc_test_microsphere_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_microsphere_ops);
}
