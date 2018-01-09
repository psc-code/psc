
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_collision.h>
#include <psc_balance.h>

#include <psc_particles_as_single.h>

#include <mrc_params.h>

#include <math.h>

// ----------------------------------------------------------------------

void *
rng(int i)
{
  return NULL;
}

double
uniform(void *dummy, double lo, double hi)
{
  return lo + (hi - lo) * random() / ((float) RAND_MAX + 1);
}

double
normal(void *dummy, double mu, double sigma)
{
  float ran1, ran2;
  do {
    ran1 = random() / ((float) RAND_MAX + 1);
    ran2 = random() / ((float) RAND_MAX + 1);
  } while (ran1 >= 1.f || ran2 >= 1.f);
	      
  return mu + sigma * sqrtf(-2.f * logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
}

// ----------------------------------------------------------------------

enum {
  KIND_TEST,
  N_KINDS,
};

struct psc_test_accel {
  // params
  int n_step;
  int n_prts;
  double L;
  double tol;
  
  // state
};

#define psc_test_accel(psc) mrc_to_subobj(psc, struct psc_test_accel)

#define VAR(x) (void *)offsetof(struct psc_test_accel, x)
static struct param psc_test_accel_descr[] = {
  { "n_step"            , VAR(n_step)              , PARAM_INT(64)         },
  { "n_prts"            , VAR(n_prts)              , PARAM_INT(131)        },
  { "L"                 , VAR(L)                   , PARAM_DOUBLE(1e10)    },

  // tolerance -- 6th order method gets to < 1e-6, but not the regular 2nd order one
  { "tol"               , VAR(tol)                 , PARAM_DOUBLE(1e-2)    },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_test_accel_create

static void
psc_test_accel_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nicell = 1;
  psc->prm.cfl = 1.;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1;
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

  psc_push_particles_set_type(psc->push_particles, "1vbec_double");
}

// ----------------------------------------------------------------------
// psc_test_accel_setup

static void
psc_test_accel_setup(struct psc *psc)
{
  struct psc_test_accel *sub = psc_test_accel(psc);

  psc->prm.nmax = sub->n_step;

  psc->domain.length[0] = sub->L;
  psc->domain.length[1] = sub->L;
  psc->domain.length[2] = sub->L;

  struct psc_kind kinds[N_KINDS] = {
    { .q = 1., .m = 1., .name = "test_species", },
  };
  psc_set_kinds(psc, N_KINDS, kinds);

  // initializes fields, particles, etc.
  psc_setup_super(psc);
}

// ----------------------------------------------------------------------
// psc_test_accel_init_field

static double
psc_test_accel_init_field(struct psc *psc, double crd[3], int m)
{
  //struct psc_test_accel *sub = psc_test_accel(psc);

  switch (m) {
  case EX: return 1.;
  case EY: return 2.;
  case EZ: return 3.;
  
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_accel_setup_particles

#define inject_particle(knd, x, y, z, ux, uy, uz, weight, a, b) do {	\
    double Vi = 1./(psc->patch[0].dx[0] * psc->patch[0].dx[1] * psc->patch[0].dx[2]); \
    Vi = 1.; /* FIXME HACK */						\
    particle_t prt;							\
    prt.xi = x - xmin;							\
    prt.yi = y - ymin;							\
    prt.zi = z - zmin;							\
    prt.pxi = ux;							\
    prt.pyi = uy;							\
    prt.pzi = uz;							\
    prt.qni_wni = psc->kinds[knd].q * weight * Vi;			\
    prt.kind_ = knd;							\
    mprts[p].push_back(prt);						\
  } while (0)

static void
psc_test_accel_setup_particles(struct psc *psc, int *nr_particles_by_patch, bool count_only)
{
  struct psc_test_accel *sub = psc_test_accel(psc);

  for (int p = 0; p < psc->nr_patches; p++) {
    nr_particles_by_patch[p] = sub->n_prts;
  }

  if (count_only) {
    return;
  }

  int test_species = KIND_TEST;

  mparticles_t mprts = psc->particles->get_as<mparticles_t>(MP_DONT_COPY);

  for (int p = 0; p < psc->nr_patches; p++) {
    struct psc_patch *patch = &psc->patch[p];
    double xmin = patch->xb[0], xmax = patch->xb[0] + patch->ldims[0] * patch->dx[0];
    double ymin = patch->xb[1], ymax = patch->xb[1] + patch->ldims[1] * patch->dx[1];
    double zmin = patch->xb[2], zmax = patch->xb[2] + patch->ldims[2] * patch->dx[2];

    for (int n = 0; n < sub->n_prts; n++) {
      double x = uniform(rng(0), xmin, xmax);
      double y = uniform(rng(0), ymin, ymax);
      double z = uniform(rng(0), zmin, zmax);
      double weight = 1.;
     
      inject_particle(test_species,
		      x, y, z,
		      0., 0., 0.,
		      weight, 0, 0 );
    }
  }

  mprts.put_as(psc->particles);
}

// ----------------------------------------------------------------------
// psc_test_accel_step

static void
psc_test_accel_step(struct psc *psc)
{
  struct psc_test_accel *sub = psc_test_accel(psc);

  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);

  int n = psc->timestep;
  double ux = 1. * (n+1);
  double uy = 2. * (n+1);
  double uz = 3. * (n+1);
  
  double tol = sub->tol;
  int failed = 0;

  mparticles_t mprts = psc->particles->get_as<mparticles_t>();

  for (int p = 0; p < psc->nr_patches; p++) {
    particle_range_t prts = mprts[p].range();
    unsigned int n_prts = prts.size();

    for (int n = 0; n < n_prts; n++) {
      particle_t *p = particle_iter_at(prts.begin, n);
      if (fabs(p->pxi - ux) > tol ||
          fabs(p->pyi - uy) > tol ||
          fabs(p->pzi - uz) > tol) {
        failed++;
	mprintf("n %d: xi [%g %g %g] pxi [%g %g %g] qni_wni %g kind %d delta %g %g %g\n", n,
		p->xi, p->yi, p->zi,
		p->pxi, p->pyi, p->pzi, p->qni_wni, p->kind_, 
		p->pxi - ux, p->pyi - uy, p->pzi - uz);
      }
    }
  }

  mprts.put_as(psc->particles, MP_DONT_COPY);

  assert(!failed);
}

// ======================================================================
// psc_test_accel_ops

struct psc_ops psc_test_accel_ops = {
  .name             = "test_accel",
  .size             = sizeof(struct psc_test_accel),
  .param_descr      = psc_test_accel_descr,
  .create           = psc_test_accel_create,
  .setup            = psc_test_accel_setup,
  .init_field       = psc_test_accel_init_field,
  .setup_particles  = psc_test_accel_setup_particles,
  .step             = psc_test_accel_step,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_accel_ops);
}
