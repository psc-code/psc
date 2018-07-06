
#include "psc_method_private.h"

#include <psc_fields_vpic.h>
#include <psc_particles_vpic.h>

#include <psc_balance.h>
#include <psc_marder.h>
#include <setup_particles.hxx>
#include <psc_particles_double.h>

#include <vpic_iface.h>

// ======================================================================
// psc_method "vpic"

struct psc_method_vpic {
  // state
  Simulation *sim;
};

#define psc_method_vpic(method) mrc_to_subobj(method, struct psc_method_vpic)

#define VAR(x) (void *)offsetof(struct psc_method_vpic, x)
static struct param psc_method_vpic_descr[] = {
  { "sim"                   , VAR(sim)                            , PARAM_PTR(NULL)   },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_method_vpic_initialize

static void
psc_method_vpic_initialize(struct psc_method *method, struct psc *psc,
			   PscMparticlesBase mprts_base)
{
  struct psc_method_vpic *sub = psc_method_vpic(method);

  auto mflds_base = PscMfieldsBase{psc->flds};
  auto& mflds = mflds_base->get_as<MfieldsVpic>(0, VPIC_MFIELDS_N_COMP);
  auto& mprts = mprts_base->get_as<MparticlesVpic>();
  
  // Do some consistency checks on user initialized fields

  mpi_printf(psc_comm(psc), "Checking interdomain synchronization\n");
  double err = mflds.synchronize_tang_e_norm_b();
  mpi_printf(psc_comm(psc), "Error = %g (arb units)\n", err);
  
  mpi_printf(psc_comm(psc), "Checking magnetic field divergence\n");
  mflds.compute_div_b_err();
  err = mflds.compute_rms_div_b_err();
  mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
  mflds.clean_div_b();
  
  // Load fields not initialized by the user

  mpi_printf(psc_comm(psc), "Initializing radiation damping fields\n");
  mflds.compute_curl_b();

  mpi_printf(psc_comm(psc), "Initializing bound charge density\n");
  mflds.clear_rhof();
  mflds.accumulate_rho_p(mprts.vmprts);
  mflds.synchronize_rho();
  mflds.compute_rhob();

  // Internal sanity checks

  mpi_printf(psc_comm(psc), "Checking electric field divergence\n");
  mflds.compute_div_e_err();
  err = mflds.compute_rms_div_e_err();
  mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
  mflds.clean_div_e();

  mpi_printf(psc_comm(psc), "Rechecking interdomain synchronization\n");
  err = mflds.synchronize_tang_e_norm_b();
  mpi_printf(psc_comm(psc), "Error = %e (arb units)\n", err);

  FieldArray *vmflds = mflds.vmflds_fields;
  Simulation_initialize(sub->sim, mprts.vmprts, vmflds);

  mprts_base->put_as(mprts);
  mflds_base->put_as(mflds, 0, VPIC_MFIELDS_N_COMP);

  // First output / stats
  
  mpi_printf(psc_comm(psc), "Performing initial diagnostics.\n");
  Simulation_diagnostics_run(sub->sim);
  psc_method_default_output(method, psc, mprts_base);

  Simulation_print_status(sub->sim);
  psc_stats_log(psc);
  psc_print_profiling(psc);
}

// ----------------------------------------------------------------------
// psc_method_vpic_output

static void
psc_method_vpic_output(struct psc_method *method, struct psc *psc,
		       PscMparticlesBase mprts)
{
  struct psc_method_vpic *sub = psc_method_vpic(method);

  // FIXME, a hacky place to do this
  Simulation_inc_step(sub->sim, psc->timestep);

  Simulation_diagnostics_run(sub->sim);
  
  if (psc->prm.stats_every > 0 && psc->timestep % psc->prm.stats_every == 0) {
    Simulation_print_status(sub->sim);
  }
  
  psc_method_default_output(NULL, psc, mprts);
}

// ----------------------------------------------------------------------
// psc_method "vpic"

struct psc_method_ops_vpic : psc_method_ops {
  psc_method_ops_vpic() {
    name                          = "vpic";
    size                          = sizeof(struct psc_method_vpic);
    param_descr                   = psc_method_vpic_descr;
    initialize                    = psc_method_vpic_initialize;
    output                        = psc_method_vpic_output;
  }
} psc_method_ops_vpic;
