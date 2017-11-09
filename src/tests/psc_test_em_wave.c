
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_collision.h>
#include <psc_balance.h>

#include <psc_fields_as_c.h>

#include <mrc_params.h>

#include <math.h>

struct psc_test_em_wave {
  // params
  double k[3];
  double B[3];
  double tol;

  // state
  double om;
  double E[3];
};

#define psc_test_em_wave(psc) mrc_to_subobj(psc, struct psc_test_em_wave)

#define VAR(x) (void *)offsetof(struct psc_test_em_wave, x)
static struct param psc_test_em_wave_descr[] = {
  { "k"               , VAR(k)                 , PARAM_DOUBLE3(2., 0., 1.)  },
  { "B"               , VAR(B)                 , PARAM_DOUBLE3(0., 1., 0.)  },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_test_em_wave_create

static void
psc_test_em_wave_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nicell = 1;
  psc->prm.cfl = 1.;

  psc->domain.gdims[0] = 16;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 16;

  psc->domain.length[0] = 2. * M_PI;
  psc->domain.length[1] = 2. * M_PI;
  psc->domain.length[2] = 2. * M_PI;

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

  psc_push_particles_set_type(psc->push_particles, "1vbec_double");
}

// ----------------------------------------------------------------------
// psc_test_em_wave_setup

static void
psc_test_em_wave_setup(struct psc *psc)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);

  double *k = sub->k, *B = sub->B, *E = sub->E;

  assert(k[0]*B[0] + k[1]*B[1] + k[2]*B[2] == 0);
  double k_norm = sqrt(sqr(k[0]) + sqr(k[1]) + sqr(k[2]));
  sub->om = /* c * */ k_norm;
  E[0] = 1./sub->om * (k[1] * B[2] - k[2] * B[1]);
  E[1] = 1./sub->om * (k[2] * B[0] - k[0] * B[2]);
  E[2] = 1./sub->om * (k[0] * B[1] - k[1] * B[0]);

  psc_setup_super(psc);

}

// ----------------------------------------------------------------------
// psc_test_em_wave_init_field

static double
psc_test_em_wave_init_field(struct psc *psc, double crd[3], int m)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);
  double *k = sub->k, *B = sub->B, *E = sub->E;
  double x = crd[0], y = crd[1], z = crd[2];

  switch (m) {
  case EX: return E[0] * sin(k[0]*x + k[1]*y + k[2]*z);
  case EY: return E[1] * sin(k[0]*x + k[1]*y + k[2]*z);
  case EZ: return E[2] * sin(k[0]*x + k[1]*y + k[2]*z);
  case HX: return B[0] * sin(k[0]*x + k[1]*y + k[2]*z);
  case HY: return B[1] * sin(k[0]*x + k[1]*y + k[2]*z);
  case HZ: return B[2] * sin(k[0]*x + k[1]*y + k[2]*z);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_em_wave_step

static void
psc_test_em_wave_step(struct psc *psc)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);

  psc_output(psc);

  struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, FIELDS_TYPE, EX, EX + 6);

  for (int p = 0; p < mflds->nr_patches; p++) {
  }

  psc_mfields_put_as(mflds, psc->flds, 0, 0);
  
  psc_push_fields_push_H(psc->push_fields, psc->flds, .5);
  psc_push_fields_step_b2(psc->push_fields, psc->flds);
  psc_push_fields_step_a(psc->push_fields, psc->flds);
}

// ======================================================================
// psc_test_em_wave_ops

struct psc_ops psc_test_em_wave_ops = {
  .name             = "test_em_wave",
  .size             = sizeof(struct psc_test_em_wave),
  .param_descr      = psc_test_em_wave_descr,
  .create           = psc_test_em_wave_create,
  .setup            = psc_test_em_wave_setup,
  .init_field       = psc_test_em_wave_init_field,
  .step             = psc_test_em_wave_step,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_em_wave_ops);
}
