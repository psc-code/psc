
#include "psc_method_private.h"

#include <psc_fields_vpic.h>
#include <psc_particles_vpic.h>

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

void
psc_method_vpic_initialize(struct psc_method *method, struct psc *psc,
			   MfieldsBase& mflds_base, MparticlesBase& mprts_base)
{
  struct psc_method_vpic *sub = psc_method_vpic(method);

  auto& mflds = mflds_base.get_as<MfieldsVpic>(0, VPIC_MFIELDS_N_COMP);
  auto& mprts = mprts_base.get_as<MparticlesVpic>();
  
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
  mpi_printf(psc_comm(psc), "Uncentering particles\n");
  sub->sim->uncenter_p(mprts.vmprts, vmflds);

  mprts_base.put_as(mprts);
  mflds_base.put_as(mflds, 0, VPIC_MFIELDS_N_COMP);
}

// ----------------------------------------------------------------------
// psc_method_vpic_print_status

void
psc_method_vpic_print_status(struct psc_method *method)
{
  struct psc_method_vpic *sub = psc_method_vpic(method);

#ifdef HAVE_VPIC
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  update_profile(rank == 0);
#endif
}

// ----------------------------------------------------------------------
// psc_method_vpic_inc_step

void
psc_method_vpic_inc_step(struct psc_method *method, int timestep)
{
  struct psc_method_vpic *sub = psc_method_vpic(method);
  Simulation* sim = sub->sim;

  sim->grid_->step++;
  assert(sim->grid_->step == timestep);
}

// ----------------------------------------------------------------------
// psc_method_vpic_diagnostics_run

void
psc_method_vpic_diagnostics_run(struct psc_method *method, struct psc *psc,
				MfieldsBase& mflds, MparticlesBase& mprts)
{
  struct psc_method_vpic *sub = psc_method_vpic(method);

  auto& particles = sub->sim->particles_;
  sub->sim->runDiag(particles);
}

// ----------------------------------------------------------------------
// psc_method "vpic"

struct psc_method_ops_vpic : psc_method_ops {
  psc_method_ops_vpic() {
    name                          = "vpic";
    size                          = sizeof(struct psc_method_vpic);
    param_descr                   = psc_method_vpic_descr;
  }
} psc_method_ops_vpic;
