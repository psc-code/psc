
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

struct psc_spitzer {
  // parameters
  double nu_ei;

  double ez0;
};

#define to_psc_spitzer(psc) mrc_to_subobj(psc, struct psc_spitzer)

#define VAR(x) (void *)offsetof(struct psc_spitzer, x)
static struct param psc_spitzer_descr[] = {
  { "nu_ei"           , VAR(nu_ei)           , PARAM_DOUBLE(1.)          },

  {},
};
#undef VAR


// ----------------------------------------------------------------------
// psc_spitzer_create

static void
psc_spitzer_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
  psc->prm.nicell = 10000;
  psc->prm.cfl = 0.98;

  psc->domain.length[0] = 4.; // no x-dependence
  psc->domain.length[1] = 4.; // no y-dependence
  psc->domain.length[2] = 1.;

  psc->domain.gdims[0] = 4;
  psc->domain.gdims[1] = 4;
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

  // psc_moments_set_type(psc->moments, "c_1st_cc");

  psc->kinds[KIND_ELECTRON].T = .0001;
  psc->kinds[KIND_ION].T = .0001;
  psc->kinds[KIND_ION].m = 1836;

  psc_push_fields_set_type(psc->push_fields, "none");
  psc_sort_set_type(psc->sort, "countsort2");
  psc_collision_set_type(psc->collision, "c");
  psc_sort_set_param_int(psc->sort, "every", 1);
  psc_collision_set_param_int(psc->collision, "every", 1);
}

// ----------------------------------------------------------------------
// psc_spitzer_setup

static void
psc_spitzer_setup(struct psc *psc)
{
  struct psc_spitzer *spitzer = to_psc_spitzer(psc);

  double Z = 1.;
  double nu = 1. / Z * spitzer->nu_ei * pow(psc->kinds[KIND_ELECTRON].T, 1.5) * 3 * sqrt(M_PI/2);
  spitzer->ez0 = .03 * spitzer->nu_ei * sqrt(psc->kinds[KIND_ELECTRON].T);
  psc_collision_set_param_double(psc->collision, "nu", nu);

  psc_setup_super(psc);

  psc->dt = .01;
}

// ----------------------------------------------------------------------
// psc_spitzer_init_field

static double
psc_spitzer_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_spitzer *spitzer = to_psc_spitzer(psc);

  switch (m) {
  case EZ: return spitzer->ez0;
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_spitzer_init_npt

static void
psc_spitzer_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt)
{
  npt->n = 1.;
}

// ======================================================================
// psc_spitzer_ops

struct psc_ops psc_spitzer_ops = {
  .name             = "spitzer",
  .size             = sizeof(struct psc_spitzer),
  .param_descr      = psc_spitzer_descr,
  .create           = psc_spitzer_create,
  .setup            = psc_spitzer_setup,
  .init_field       = psc_spitzer_init_field,
  .init_npt         = psc_spitzer_init_npt,
};

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_spitzer_ops);
}
