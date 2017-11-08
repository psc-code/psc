
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_collision.h>
#include <psc_balance.h>

#include <psc_particles_as_single.h>

#include <mrc_params.h>

#include <math.h>

struct psc_test_em_wave {
  // params
  double k[3];
  double tol;
};

#define psc_test_em_wave(psc) mrc_to_subobj(psc, struct psc_test_em_wave)

#define VAR(x) (void *)offsetof(struct psc_test_em_wave, x)
static struct param psc_test_em_wave_descr[] = {
  { "k"               , VAR(k)                 , PARAM_DOUBLE3(2., 0., 1.)  },
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
// psc_test_em_wave_init_field

static double
psc_test_em_wave_init_field(struct psc *psc, double crd[3], int m)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);
  double *k = sub->k;
  double k_norm_inv = -1. / sqrt(sqr(k[0]) + sqr(k[1]) + sqr(k[2]));
  double b[3] = { 0., 1., 0. };
  double e[3] = { k_norm_inv * (k[1] * b[2] - k[2] * b[1]),
		  k_norm_inv * (k[2] * b[0] - k[0] * b[2]),
		  k_norm_inv * (k[0] * b[1] - k[1] * b[0]) };
  double x = crd[0], y = crd[1], z = crd[2];

  switch (m) {
  case EX: return e[0] * sin(k[0]*x + k[1]*y + k[2]*z);
  case EY: return e[1] * sin(k[0]*x + k[1]*y + k[2]*z);
  case EZ: return e[2] * sin(k[0]*x + k[1]*y + k[2]*z);
  case HX: return b[0] * sin(k[0]*x + k[1]*y + k[2]*z);
  case HY: return b[1] * sin(k[0]*x + k[1]*y + k[2]*z);
  case HZ: return b[2] * sin(k[0]*x + k[1]*y + k[2]*z);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_em_wave_step

static void
psc_test_em_wave_step(struct psc *psc)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);

  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);

  int n = psc->timestep;
  double ux = 1. * (n+1);
  double uy = 2. * (n+1);
  double uz = 3. * (n+1);
  
  double tol = sub->tol;
  int failed = 0;

  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, PARTICLE_TYPE, 0);

  for (int p = 0; p < psc->nr_patches; p++) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    unsigned int n_prts = particle_range_size(prts);

    for (int n = 0; n < n_prts; n++) {
      particle_t *p = particle_iter_at(prts.begin, n);
      if (fabs(p->pxi - ux) > tol ||
          fabs(p->pyi - uy) > tol ||
          fabs(p->pzi - uz) > tol) {
        failed++;
	mprintf("n %d: xi [%g %g %g] pxi [%g %g %g] qni_wni %g kind %d delta %g %g %g\n", n,
		p->xi, p->yi, p->zi,
		p->pxi, p->pyi, p->pzi, p->qni_wni, p->kind, 
		p->pxi - ux, p->pyi - uy, p->pzi - uz);
      }
    }
  }

  psc_mparticles_put_as(mprts, psc->particles, MP_DONT_COPY);

  assert(!failed);
}

// ======================================================================
// psc_test_em_wave_ops

struct psc_ops psc_test_em_wave_ops = {
  .name             = "test_em_wave",
  .size             = sizeof(struct psc_test_em_wave),
  .param_descr      = psc_test_em_wave_descr,
  .create           = psc_test_em_wave_create,
  .init_field       = psc_test_em_wave_init_field,
  //  .step             = psc_test_em_wave_step,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_em_wave_ops);
}
