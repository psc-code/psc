
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_push_particles.h>
#include <psc_bnd.h>
#include <psc_bnd_particles.h>
#include <psc_bnd_photons.h>
#include <psc_particles_as_c.h>

#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>

struct psc_photon_test {
  // parameters
  int n_in_cell;
  double rho_0;
  double tau_0;
  double k_L;
};

#define to_psc_photon_test(psc) mrc_to_subobj(psc, struct psc_photon_test)

#define VAR(x) (void *)offsetof(struct psc_photon_test, x)
static struct param psc_photon_test_descr[] = {
  { "n_in_cell"       , VAR(n_in_cell)         , PARAM_INT(1000)      },
  { "rho_0"           , VAR(rho_0)             , PARAM_DOUBLE(1e-6)   },
  { "tau_0"           , VAR(tau_0)             , PARAM_DOUBLE(1.5e-6) },
  { "k_L"             , VAR(k_L)               , PARAM_DOUBLE(1.)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_photon_test_create

static void
psc_photon_test_create(struct psc *psc)
{
  psc->prm.nmax = 201;
  psc->prm.nicell = 1;

  psc->domain.length[0] = 10 * 1e-6;
  psc->domain.length[1] = 10 * 1e-6;
  psc->domain.length[2] = 10 * 1e-6;

  psc->domain.gdims[0] = 64;
  psc->domain.gdims[1] = 64;
  psc->domain.gdims[2] = 64;
}

// ----------------------------------------------------------------------
// psc_photon_test_init_npt

static void
psc_photon_test_init_npt(struct psc *psc, int kind, double x[3],
			 struct psc_particle_npt *npt)
{
  double *dx = psc->patch[0].dx;
  double xc[3];
  for (int d = 0; d < 3; d++) {
    xc[d] = psc->domain.length[d] / 2. / psc->coeff.ld;
  }

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->p[0] = 10.;
    if (x[0] >= xc[0] && x[0] < xc[0] + dx[0] &&
	x[1] >= xc[1] && x[1] < xc[1] + dx[1] &&
	x[2] >= xc[2] && x[2] < xc[2] + dx[2]) {
      npt->n = 1.;
    }
    break;
  }
}

static void
psc_event_generator_emission(struct psc *psc)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, "c", 0);
  
  psc_foreach_patch(psc, p) {
    // get array of particles on this patch
    particle_range_t prts = particle_range_mprts(mprts, p);
  
    photons_t *photons = &psc->mphotons->p[p];
    // and iterative over all particles in the array
    PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
      particle_t *part = particle_iter_deref(prt_iter);

      if (part->qni >= 0. || part->pxi < 1.5) {
	continue;
      }
      // only electrons with pxi > 1.5
      long r = random();
      if (r < .2 * RAND_MAX) {
	photons_realloc(photons, photons->nr + 1);
	int new_nr = photons->nr++;
	photon_t *ph = &photons->photons[new_nr];
	ph->x[0] = part->xi;
	ph->x[1] = part->yi;
	ph->x[2] = part->zi;
	ph->p[0] = part->pxi / 2.;
	ph->wni = 1.;
	part->pxi /= 2.;
      }
    }
  }

  psc_mparticles_put_as(mprts, psc->particles, 0);
}

static void
psc_photon_test_step(struct psc *psc)
{
  psc_output(psc);
 
  psc_event_generator_emission(psc);
  
  // field propagation n*dt -> (n+0.5)*dt
  psc_push_fields_step_a(psc->push_fields, psc->flds);
  
  // particle propagation n*dt -> (n+1.0)*dt
  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);
  
  psc_push_photons_run(psc->mphotons);
  psc_bnd_photons_exchange(psc->bnd_photons, psc->mphotons);
  
  psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
  psc_bnd_fill_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
  
  // field propagation (n+0.5)*dt -> (n+1.0)*dt
  psc_push_fields_step_b1(psc->push_fields, psc->flds);
  psc_push_fields_step_b2(psc->push_fields, psc->flds);
}

// ======================================================================
// psc_photon_test_ops

struct psc_ops psc_photon_test_ops = {
  .name             = "photon_test",
  .size             = sizeof(struct psc_photon_test),
  .param_descr      = psc_photon_test_descr,
  .create           = psc_photon_test_create,
  .init_npt         = psc_photon_test_init_npt,
  .step             = psc_photon_test_step,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_photon_test_ops);
}
