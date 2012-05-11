
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

// plasma oscillation
// src/psc_es1 --psc_output_particles_type ascii --psc_output_particles_every_step 1 --mrc_io_type ascii --pfield_step 1 --psc_diag_every_step 1

// two-stream
// src/psc_es1 --psc_output_particles_type ascii --psc_output_particles_every_step 10 --mrc_io_type ascii --write_tfield no --pfield_step 10 --psc_diag_every_step 10 --nr_kinds 2 --v0_1 1. --v0_2 -1. --cc 5. --nicell 32

struct psc_es1_species {
  int nlg; // number of loading groups
  double q;
  double m;
  double mode;
  double v0;
  double x1;
  double v1;
  double thetax;
  double thetav;
};

#define MAX_KINDS (2)

struct psc_es1 {
  // parameters
  double om_pe;
  double om_ce;
  struct psc_es1_species species[MAX_KINDS];
};

#define to_psc_es1(psc) mrc_to_subobj(psc, struct psc_es1)

#define VAR(x) (void *)offsetof(struct psc_es1, x)
static struct param psc_es1_descr[] = {
  { "om_pe"         , VAR(om_pe)            , PARAM_DOUBLE(1.)            },
  { "om_ce"         , VAR(om_ce)            , PARAM_DOUBLE(2.)            },

  { "nlg_1"         , VAR(species[0].nlg)   , PARAM_INT(1)                },
  { "q_1"           , VAR(species[0].q)     , PARAM_DOUBLE(-1.)           },
  { "m_1"           , VAR(species[0].m)     , PARAM_DOUBLE(1.)            },
  { "mode_1"        , VAR(species[0].mode)  , PARAM_DOUBLE(1.)            },
  { "v0_1"          , VAR(species[0].v0)    , PARAM_DOUBLE(0.)            },
  { "x1_1"          , VAR(species[0].x1)    , PARAM_DOUBLE(.001)          },
  { "v1_1"          , VAR(species[0].v1)    , PARAM_DOUBLE(0.)            },
  { "thetax_1"      , VAR(species[0].thetax), PARAM_DOUBLE(0.)            },
  { "thetav_1"      , VAR(species[0].thetav), PARAM_DOUBLE(0.)            },

  { "nlg_2"         , VAR(species[1].nlg)   , PARAM_INT(1)                },
  { "q_2"           , VAR(species[1].q)     , PARAM_DOUBLE(-1.)           },
  { "m_2"           , VAR(species[1].m)     , PARAM_DOUBLE(1.)            },
  { "mode_2"        , VAR(species[1].mode)  , PARAM_DOUBLE(1.)            },
  { "v0_2"          , VAR(species[1].v0)    , PARAM_DOUBLE(0.)            },
  { "x1_2"          , VAR(species[1].x1)    , PARAM_DOUBLE(.001)          },
  { "v1_2"          , VAR(species[1].v1)    , PARAM_DOUBLE(0.)            },
  { "thetax_2"      , VAR(species[1].thetax), PARAM_DOUBLE(0.)            },
  { "thetav_2"      , VAR(species[1].thetav), PARAM_DOUBLE(0.)            },

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
  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc_set_kinds(psc, 1, NULL);
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

  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;

  psc_moments_set_type(psc->moments, "c_1st_cc");
}

// ----------------------------------------------------------------------
// psc_es1_init_field

static double
psc_es1_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_es1 *es1 = to_psc_es1(psc);


  switch (m) {
  case EZ: ;
    double ez = 0;
    for (int kind = 0; kind < psc->nr_kinds; kind++) {
      struct psc_es1_species *s = &es1->species[kind];
      double theta = 2. * M_PI * s->mode / psc->domain.length[2] * x[2];
      ez += s->q * s->x1 * cos(theta + s->thetax);
    }
    return ez;

  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_es1_setup_particles

void
psc_es1_setup_particles(struct psc *psc, int *nr_particles_by_patch,
		       bool count_only)
{
  struct psc_es1 *es1 = to_psc_es1(psc);

  if (count_only) {
    psc_foreach_patch(psc, p) {
      int *ldims = psc->patch[p].ldims;
      int n = 0;
      for (int kind = 0; kind < psc->nr_kinds; kind++) {
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

  double l = psc->domain.length[2];

  psc_foreach_patch(psc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);

    int il1 = 0;
    for (int kind = 0; kind < psc->nr_kinds; kind++) {
      struct psc_es1_species *s = &es1->species[kind];
      int *ldims = psc->patch[p].ldims;
      int n = ldims[0] * ldims[1] * ldims[2] * psc->prm.nicell;
      int ngr = n / s->nlg;
      //double lg = l / nlg;
      double ddx = l / n;
      for (int i = 0; i < ngr; i++) {
	particle_t *p = particles_get_one(pp, il1 + i);
	double x0 = (i + .5) * ddx;
	
	p->zi = x0;
	p->pzi = s->v0 / psc->prm.cc;
	p->qni = s->q;
	p->mni = s->m;
	p->wni = 1.;
      }
      for (int i = 0; i < n; i++) {
	particle_t *p = particles_get_one(pp, il1++);
	double theta = 2. * M_PI * s->mode * p->zi / l;
	p->zi  += s->x1 * cos(theta + s->thetax);
	p->pzi += s->v1 * sin(theta + s->thetav);
      }
    }
    pp->n_part = il1;
    assert(pp->n_part == nr_particles_by_patch[p]);
  }
  psc_mparticles_put_cf(particles, psc->particles, 0);
}

// ======================================================================
// psc_es1_ops

struct psc_ops psc_es1_ops = {
  .name             = "es1",
  .size             = sizeof(struct psc_es1),
  .param_descr      = psc_es1_descr,
  .create           = psc_es1_create,
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
