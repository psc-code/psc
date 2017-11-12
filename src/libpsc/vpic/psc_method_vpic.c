
#include "psc_method_private.h"

#include <psc_fields_vpic.h>
#include <psc_particles_vpic.h>
#include <psc_push_particles_vpic.h>


// ======================================================================
// psc_method "vpic"

// ----------------------------------------------------------------------
// psc_method_vpic_initialize

static void
psc_method_vpic_initialize(struct psc_method *method, struct psc *psc)
{
  // Do some consistency checks on user initialized fields

  mpi_printf(psc_comm(psc), "Checking interdomain synchronization\n");
  double err = psc_mfields_synchronize_tang_e_norm_b(psc->flds);
  mpi_printf(psc_comm(psc), "Error = %g (arb units)\n", err);
  
  mpi_printf(psc_comm(psc), "Checking magnetic field divergence\n");
  psc_mfields_compute_div_b_err(psc->flds);
  err = psc_mfields_compute_rms_div_b_err(psc->flds);
  mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
  psc_mfields_clean_div_b(psc->flds);
  
  // Load fields not initialized by the user

  mpi_printf(psc_comm(psc), "Initializing radiation damping fields\n");
  psc_mfields_compute_curl_b(psc->flds);

  mpi_printf(psc_comm(psc), "Initializing bound charge density\n");
  psc_mfields_clear_rhof(psc->flds);
  psc_mfields_accumulate_rho_p(psc->flds, psc->particles);
  psc_mfields_synchronize_rho(psc->flds);
  psc_mfields_compute_rhob(psc->flds);

  // Internal sanity checks

  mpi_printf(psc_comm(psc), "Checking electric field divergence\n");
  psc_mfields_compute_div_e_err(psc->flds);
  err = psc_mfields_compute_rms_div_e_err(psc->flds);
  mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
  psc_mfields_clean_div_e(psc->flds);

  mpi_printf(psc_comm(psc), "Rechecking interdomain synchronization\n");
  err = psc_mfields_synchronize_tang_e_norm_b(psc->flds);
  mpi_printf(psc_comm(psc), "Error = %e (arb units)\n", err);

  mpi_printf(psc_comm(psc), "Uncentering particles\n");
  psc_push_particles_stagger(psc->push_particles, psc->particles, psc->flds);

  // First output / stats
  
  mpi_printf(psc_comm(psc), "Performing initial diagnostics.\n");
  vpic_diagnostics();
  psc_output_default(psc);

  vpic_print_status();
  psc_stats_log(psc);
  psc_print_profiling(psc);
}

// ----------------------------------------------------------------------
// psc_method "vpic"

struct psc_method_ops psc_method_ops_vpic = {
  .name                = "vpic",
  .initialize          = psc_method_vpic_initialize,
};
