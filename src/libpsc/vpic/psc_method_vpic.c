
#include "psc_method_private.h"

#include <psc_fields_vpic.h>
#include <psc_particles_vpic.h>
#include <psc_push_particles_vpic.h>

#include <vpic_iface.h>


// ======================================================================
// psc_method "vpic"

static void
psc_method_vpic_initialize(struct psc_method *method, struct psc *psc)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(psc->flds)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(psc->particles)->vmprts;
  struct vpic_push_particles *vpushp = psc_push_particles_vpic(psc->push_particles)->vpushp;

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
  vpic_simulation_init2(vpushp, vmflds, vmprts);

  mpi_printf(psc_comm(psc), "Initialization complete.\n");
}

// ----------------------------------------------------------------------
// psc_method "vpic"

struct psc_method_ops psc_method_ops_vpic = {
  .name                = "vpic",
  .initialize          = psc_method_vpic_initialize,
};
