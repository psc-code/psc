
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_double.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

// ======================================================================
// psc_es1
//
// This code replicates, to some extent, the ES1 code from Birdsall/Langdon
//
// However, PSC is an electromagnetic code, so changes in charge density don't
// cause instantaneous changes in electric potential/field everywhere, information
// propagates at the speed of light (which is a parameter)


/* plasma oscillation
   
   src/psc_es1					\
   --psc_output_particles_type ascii		\
   --psc_output_particles_every_step 1		\
   --mrc_io_type ascii				\
   --pfield_step 1				\
   --psc_diag_every_step 1
*/

/* two-stream instability
   
   src/psc_es1								\
   --psc_output_particles_type ascii					\
   --psc_output_particles_every_step 10					\
   --mrc_io_type ascii							\
   --write_tfield no --pfield_step 10					\
   --psc_diag_every_step 10						\
   --particle_kinds e1,e2						\
   --v0_1 1. --v0_2 -1. --cc 6. --nicell 32
*/

struct psc_es1_species {
  double vt1; // thermal velocity for random velocities
  double vt2; // thermal velocity for order velocities
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
  struct psc_es1_species species[MAX_KINDS];
};

#define to_psc_es1(psc) mrc_to_subobj(psc, struct psc_es1)

#define VAR(x) (void *)offsetof(struct psc_es1, x)
static struct param psc_es1_descr[] = {
  { "vt1_1"         , VAR(species[0].vt1)   , PARAM_DOUBLE(0.)            },
  { "vt2_1"         , VAR(species[0].vt2)   , PARAM_DOUBLE(0.)            },
  { "nlg_1"         , VAR(species[0].nlg)   , PARAM_INT(1)                },
  { "q_1"           , VAR(species[0].q)     , PARAM_DOUBLE(-1.)           },
  { "m_1"           , VAR(species[0].m)     , PARAM_DOUBLE(1.)            },
  { "mode_1"        , VAR(species[0].mode)  , PARAM_DOUBLE(1.)            },
  { "v0_1"          , VAR(species[0].v0)    , PARAM_DOUBLE(0.)            },
  { "x1_1"          , VAR(species[0].x1)    , PARAM_DOUBLE(.001)          },
  { "v1_1"          , VAR(species[0].v1)    , PARAM_DOUBLE(0.)            },
  { "thetax_1"      , VAR(species[0].thetax), PARAM_DOUBLE(0.)            },
  { "thetav_1"      , VAR(species[0].thetav), PARAM_DOUBLE(0.)            },

  { "vt1_2"         , VAR(species[1].vt1)   , PARAM_DOUBLE(0.)            },
  { "vt2_2"         , VAR(species[1].vt2)   , PARAM_DOUBLE(0.)            },
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
  psc_default_dimensionless(psc);

  psc_set_kinds(psc, 1, NULL);
  psc->prm.nicell = 4;
  psc->prm.cfl = 0.98;
  psc->prm.nmax = 16000;

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
}

// ----------------------------------------------------------------------
// psc_es1_init_field

static double
psc_es1_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_es1 *es1 = to_psc_es1(psc);

  switch (m) {
  case EZ: {
    double ez = 0;
    for (int kind = 0; kind < psc->nr_kinds; kind++) {
      struct psc_es1_species *s = &es1->species[kind];
      double theta = 2. * M_PI * s->mode / psc->domain.length[2] * x[2];
      ez -= s->q * s->x1 * cos(theta + s->thetax);
    }
    return ez;
  }
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ranf
//
// return random number between 0 and 1

static double
ranf()
{
  return (double) random() / RAND_MAX;

}

// ----------------------------------------------------------------------
// psc_es1_init_species

static void
psc_es1_init_species(struct psc *psc, int kind, struct psc_es1_species *s,
		     struct psc_mparticles *mprts, int p, int *p_il1)
{
  int *ldims = psc->patch[p].ldims;
  int n = ldims[0] * ldims[1] * ldims[2] * psc->prm.nicell; // FIXME
  double l = psc->domain.length[2];

  // load first of nlg groups of particles.
  int ngr = n / s->nlg;
  double lg = l / s->nlg;
  double ddx = l / n;

  double *x = (double*) calloc(n, sizeof(*x));
  double *vx = (double*) calloc(n, sizeof(*vx));

  // load evenly spaced, with drift.
  // also does cold case.
  for (int i = 0; i < ngr; i++) {
    double x0 = (i + .5) * ddx;
    int i1 = i; // + il1;
    x [i1] = x0;
    vx[i1] = s->v0;
  }

  // load order velocities in vx ("quiet start", or at least subdued).
  // is set up for maxwellian*v*nv2, but can do any smooth distribution.
  // hereafter, ngr is prferably a power of 2
  if (s->vt2 != 0.) {
    // first store indefinite integral of distribution function in x array.
    // use midpoint rule -simple and quite accurate.
    double vmax = 5. * s->vt2;
    double dv = 2.*vmax / (n-1);
    double vvnv2 = 1.;
    x[0] = 0.;
    for (int i = 2; i <= n; i++) {
      double vv = ((i-1.5) * dv - vmax) / s->vt2;
#if 0
      if (s->nv2 != 0.) {
	vvnv2 = pow(vv, s->nv2);
      }
#endif
      double fv = vvnv2*exp(-.5 * sqr(vv));
      int i1 = i-1;
      x[i1] = x[i1-1] + fmax(fv, 0.);
      printf("i1 %d x %g\n", i1, x[i1]);
    }
    
    // for evenly spaced (half-integer multiples) valus of the integral,
    // find corresponding velocities by inverse linear interpolation
    double df = x[n-1] / ngr;
    int i1 = 0;
    int j = 0;
    for (int i = 1;i <= ngr; i++) {
      double fv = (i - .5) * df;
      while (fv >= x[j+1]) {
	j++;
	assert(j < n - 1);
      }
      double vv = dv * (j + (fv - x[j])/(x[j+1] - x[j])) - vmax;
      vx[i1] += vv;
      i1++;
    }

    // for ordered velocities, scramble positions to reduce correlations
    float xs = 0.0;
    for (int i=1; i <= ngr; i++)
    {
      i1 = i - 1;
      x[i1] = xs * lg + 0.5*ddx;
      float xsi = 1.0;
      do {
	xsi *= 0.5;
	xs -= xsi;
      } while (xs >= 0.0);
      xs += 2.0*xsi;
    }
  }

  // if magnetized, rotate (vx, 0) into (vx, vy)
#if 0
  if (s->wc != 0.) {
    assert(0);
  }
#endif

  // copy first group into rest of groups.
  if (s->nlg > 1) {
    int xs = 0.;
    for (int i = ngr + 1; i <= n; i += ngr) {
      xs += lg;
      for (int j = 1; j <= ngr; j++) {
	int i1 = j - 1 - 1;
	int i2 = i1 + i - 1 - 1;
	x [i2] = x [i1] + xs;
	vx[i2] = vx[i1];
#if 0
	if (s->wc != 0.) {
	  vy[i2] = vy[i1];
	}
#endif
      }
    }
  }

  // add random maxwellian.
  if (s->vt1 != 0.) {
    for (int i = 1; i < n; i++) {
      int i1 = i - 1;
      for (int j = 1; j < 12.; j++) {
#if 0
	if (s->wc != 0.) {
	} 
#endif
	vx[i1] += s->vt1 * (ranf() - .5);
      }
    }
  }

  // add perturbation.p5
  for (int i = 0; i < n; i++) {
    double theta = 2. * M_PI * s->mode * x[i] / l;
    x [i] += s->x1 * cos(theta + s->thetax);
    vx[i] += s->v1 * sin(theta + s->thetav);
  }

  assert(0);
#if 0
  // copy to PSC data structure
  for (int i = 0; i < n; i++) {
    if (x[i] < 0.) x[i] += l;
    if (x[i] >= l) x[i] -= l;

    particle_t prt = {};
    prt.zi  = x[i];
    prt.pzi = vx[i] / psc->prm.cc;
    prt.qni = s->q;
    prt.mni = s->m;
    prt.wni = 1.;
    prt.kind = kind;
    mprts[p].push_back(prt);
  }
#endif
  
  free(x);
  free(vx);
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

  mparticles_t mprts = psc->particles->get_as<mparticles_t>();
  assert(mprts.n_patches() == 1);
  psc_foreach_patch(psc, p) {
    int il1 = 0;
    for (int kind = 0; kind < psc->nr_kinds; kind++) {
      psc_es1_init_species(psc, kind, &es1->species[kind], mprts.mprts(), p, &il1);
    }
  }
  mprts.put_as(psc->particles);
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
