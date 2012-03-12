
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_moments.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_c.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>


struct psc_es1 {
  // parameters
  double om_pe;
  double om_ce;

  // calculated from the above
};

#define to_psc_es1(psc) mrc_to_subobj(psc, struct psc_es1)

#define VAR(x) (void *)offsetof(struct psc_es1, x)
static struct param psc_es1_descr[] = {
  { "om_pe"         , VAR(om_pe)           , PARAM_DOUBLE(1.)            },
  { "om_ce"         , VAR(om_ce)           , PARAM_DOUBLE(2.)            },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_es1_create

static void
psc_es1_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.nmax = 16000;
  psc->prm.cpum = 5*24.0*60*60;
  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc->prm.nr_kinds = 1;
  psc->prm.nicell = 4;
  psc->prm.cfl = 0.98;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 1.; // no y-dependence
  psc->domain.length[2] = 2. * M_PI;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 32;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part[2] = BND_PART_PERIODIC;

  psc_moments_set_type(psc->moments, "1st_cc");
}

// ----------------------------------------------------------------------
// psc_es1_setup
//
// the parameters are now set, calculate quantities to initialize fields,
// particles

static void
psc_es1_setup(struct psc *psc)
{
  //  struct psc_es1 *es1 = to_psc_es1(psc);

  psc_setup_default(psc);
}

// ----------------------------------------------------------------------
// psc_es1_init_field

static double
psc_es1_init_field(struct psc *psc, double x[3], int m)
{
  //  struct psc_es1 *es1 = to_psc_es1(psc);
  const double mode = 1;
  const double x1 = .001;
  double theta = 2. * M_PI * mode / psc->domain.length[2] * x[2];

  switch (m) {
  case EZ: return - x1 * cos(theta);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_es1_setup_particles

void
psc_es1_setup_particles(struct psc *psc, int *nr_particles_by_patch,
		       bool count_only)
{
  if (count_only) {
    psc_foreach_patch(psc, p) {
      int *ldims = psc->patch[p].ldims;
      int n = 0;
      for (int kind = 0; kind < psc->prm.nr_kinds; kind++) {
	n += ldims[0] * ldims[1] * ldims[2] * psc->prm.nicell;
      }
      nr_particles_by_patch[p] = n;
    }
    return;
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (psc->prm.seed_by_time) {
    srandom(10*rank + time(NULL));
  } else {
    srandom(rank);
  }

  mparticles_t *particles = psc_mparticles_get_cf(psc->particles, MP_DONT_COPY);

  const double l = psc->domain.length[2];
  const int nlg = 1; // number of loading groups
  const double v0 = 0;
  const double q = -1.;
  const double m = 1.;
  const double x1 = 0.001;
  const double v1 = 0.;
  const double thetax = 0.;
  const double thetav = 0.;
  const int mode = 1;

  psc_foreach_patch(psc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    int il1 = 0;

    int *ldims = psc->patch[p].ldims;
    int n = ldims[0] * ldims[1] * ldims[2] * psc->prm.nicell;
    int ngr = n / nlg;
    //double lg = l / nlg;
    double ddx = l / n;
    for (int i = 0; i < ngr; i++) {
      particle_t *p = particles_get_one(pp, il1 + i);
      double x0 = (i + .5) * ddx;

      p->zi = x0;
      p->pzi = v0;
      p->qni = q;
      p->mni = m;
      p->wni = 1.;
    }
    for (int i = 0; i < n; i++) {
      particle_t *p = particles_get_one(pp, il1 + i);
      double theta = 2. * M_PI * mode * p->zi / l;
      p->zi  += x1 * cos(theta + thetax);
      p->pzi += v1 * sin(theta + thetav);
    }
    pp->n_part = ngr;
    assert(pp->n_part == nr_particles_by_patch[p]);
  }
  psc_mparticles_put_cf(particles, psc->particles);
}

// ======================================================================
// psc_es1_ops

struct psc_ops psc_es1_ops = {
  .name             = "es1",
  .size             = sizeof(struct psc_es1),
  .param_descr      = psc_es1_descr,
  .create           = psc_es1_create,
  .setup            = psc_es1_setup,
  .init_field       = psc_es1_init_field,
  .setup_particles  = psc_es1_setup_particles,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_es1_ops);
}
