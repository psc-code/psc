
#include "psc_diag_item_private.h"
#include "psc_particles_double.h"

#include <math.h>

static void
psc_diag_item_particle_energy_run(struct psc_diag_item *item,
				  struct psc *psc, PscMparticlesBase mprts_base,
				  double *result)
{
  auto& mprts = mprts_base->get_as<MparticlesDouble>();

  const Grid_t& grid = psc->grid();
  double fnqs = grid.fnqs;
  double fac = grid.domain.dx[0] * grid.domain.dx[1] * grid.domain.dx[2];

  for (int p = 0; p < mprts.n_patches(); p++) {
    for (auto& prt : mprts[p]) {
      double gamma = sqrt(1.f + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
      double Ekin = (gamma - 1.) * mprts.prt_mni(prt) * mprts.prt_wni(prt) * fnqs;
      double qni = mprts.prt_qni(prt);
      if (qni < 0.) {
	result[0] += Ekin * fac;
      } else if (qni > 0.) {
	result[1] += Ekin * fac;
      } else {
	assert(0);
      }
    }
  }

  mprts_base->put_as(mprts, MP_DONT_COPY);
}

// ======================================================================
// psc_diag_item_particle_energy

struct psc_diag_item_ops_p : psc_diag_item_ops {
  psc_diag_item_ops_p() {
    name      = "particle_energy";
    run       = psc_diag_item_particle_energy_run;
    nr_values = 2;
    title[0]  = "E_electron";
    title[1]  = "E_ion";
  }
} psc_diag_item_particle_energy_ops;

