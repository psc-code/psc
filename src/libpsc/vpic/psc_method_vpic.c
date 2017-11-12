
#include "psc_method_private.h"

#include <psc_fields_vpic.h>
#include <psc_particles_vpic.h>
#include <psc_push_particles_vpic.h>

#include <vpic_iface.h>

// ======================================================================
// psc_method "vpic"

// ----------------------------------------------------------------------
// psc_method_vpic_do_setup

static void
psc_method_vpic_do_setup(struct psc_method *method, struct psc *psc)
{
}

// ----------------------------------------------------------------------
// psc_method_vpic_initialize

static void
psc_method_vpic_initialize(struct psc_method *method, struct psc *psc)
{
  struct psc_mparticles *mprts = psc->particles;
  psc_mfields_view(psc->flds);
  struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, "vpic", 0, 0);
  
  // Do some consistency checks on user initialized fields

  mpi_printf(psc_comm(psc), "Checking interdomain synchronization\n");
  double err = psc_mfields_synchronize_tang_e_norm_b(mflds);
  mpi_printf(psc_comm(psc), "Error = %g (arb units)\n", err);
  
  mpi_printf(psc_comm(psc), "Checking magnetic field divergence\n");
  psc_mfields_compute_div_b_err(mflds);
  err = psc_mfields_compute_rms_div_b_err(mflds);
  mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
  psc_mfields_clean_div_b(mflds);
  
  // Load fields not initialized by the user

  mpi_printf(psc_comm(psc), "Initializing radiation damping fields\n");
  psc_mfields_compute_curl_b(mflds);

  mpi_printf(psc_comm(psc), "Initializing bound charge density\n");
  psc_mfields_clear_rhof(mflds);
  psc_mfields_accumulate_rho_p(mflds, mprts);
  psc_mfields_synchronize_rho(mflds);
  psc_mfields_compute_rhob(mflds);

  // Internal sanity checks

  mpi_printf(psc_comm(psc), "Checking electric field divergence\n");
  psc_mfields_compute_div_e_err(mflds);
  err = psc_mfields_compute_rms_div_e_err(mflds);
  mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
  psc_mfields_clean_div_e(mflds);

  mpi_printf(psc_comm(psc), "Rechecking interdomain synchronization\n");
  err = psc_mfields_synchronize_tang_e_norm_b(mflds);
  mpi_printf(psc_comm(psc), "Error = %e (arb units)\n", err);

  mpi_printf(psc_comm(psc), "Uncentering particles\n");
  psc_push_particles_stagger(psc->push_particles, mprts, mflds);

  psc_mfields_put_as(mflds, psc->flds, 0, 9);

  // First output / stats
  
  mpi_printf(psc_comm(psc), "Performing initial diagnostics.\n");
  vpic_diagnostics();
  psc_method_default_output(method, psc);

  vpic_print_status();
  psc_stats_log(psc);
  psc_print_profiling(psc);
}

// ----------------------------------------------------------------------
// psc_method_vpic_output

static void
psc_method_vpic_output(struct psc_method *method, struct psc *psc)
{
  // FIXME, a hacky place to do this
  vpic_inc_step(psc->timestep);

  vpic_diagnostics();
  
  if (psc->prm.stats_every > 0 && psc->timestep % psc->prm.stats_every == 0) {
    vpic_print_status();
  }
  
  psc_method_default_output(NULL, psc);
}

// ----------------------------------------------------------------------
// psc_method "vpic"

struct psc_method_ops psc_method_ops_vpic = {
  .name                = "vpic",
  .do_setup            = psc_method_vpic_do_setup,
  .initialize          = psc_method_vpic_initialize,
  .output              = psc_method_vpic_output,
};
