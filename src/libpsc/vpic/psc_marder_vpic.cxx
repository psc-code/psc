
#include "psc_marder_vpic.h"

#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_marder_vpic_clean_div_e

static void
psc_marder_vpic_clean_div_e(struct psc_marder *marder, struct psc_mfields *mflds,
			    struct psc_mparticles *mprts)
{
  mpi_printf(psc_marder_comm(marder), "Divergence cleaning electric field\n");
  
  psc_mfields_clear_rhof(mflds);
  psc_mfields_accumulate_rho_p(mflds, mprts);
  psc_mfields_synchronize_rho(mflds);
  
  for (int round = 0; round < marder->num_div_e_round; round++ ) {
    psc_mfields_compute_div_e_err(mflds);
    if (round == 0 || round == marder->num_div_e_round - 1) {
      double err = psc_mfields_compute_rms_div_e_err(mflds);
      mpi_printf(psc_marder_comm(marder), "%s rms error = %e (charge/volume)\n",
		 round == 0 ? "Initial" : "Cleaned", err);
    }
    psc_mfields_clean_div_e(mflds);
  }
}

// ----------------------------------------------------------------------
// psc_marder_vpic_clean_div_b

static void
psc_marder_vpic_clean_div_b(struct psc_marder *marder, struct psc_mfields *mflds)
{
  mpi_printf(psc_marder_comm(marder), "Divergence cleaning magnetic field\n");
  
  for (int round = 0; round < marder->num_div_b_round; round++) {
    psc_mfields_compute_div_b_err(mflds);
    if (round == 0 || round == marder->num_div_b_round - 1) {
      double err = psc_mfields_compute_rms_div_b_err(mflds);
      mpi_printf(psc_marder_comm(marder), "%s rms error = %e (charge/volume)\n",
		 round == 0 ? "Initial" : "Cleaned", err);
    }
    psc_mfields_clean_div_b(mflds);
  }
}

// ----------------------------------------------------------------------
// psc_marder_vpic_run

static void
psc_marder_vpic_run(struct psc_marder *marder,
		    struct psc_mfields *mflds_base,
		    struct psc_mparticles *mprts_base)
{
  struct psc *psc = ppsc; // FIXME

  int clean_div_e_interval = marder->clean_div_e_interval;
  int clean_div_b_interval = marder->clean_div_b_interval;
  int sync_shared_interval = marder->sync_shared_interval;
  bool clean_div_e = (clean_div_e_interval > 0 && psc->timestep % clean_div_e_interval == 0);
  bool clean_div_b = (clean_div_b_interval > 0 && psc->timestep % clean_div_b_interval == 0);
  bool sync_shared = (sync_shared_interval > 0 && psc->timestep % sync_shared_interval == 0);

  if (!(clean_div_e || clean_div_b || sync_shared)) {
    return;
  }
  
  mfields_vpic_t mf = mflds_base->get_as<mfields_vpic_t>(EX, VPIC_MFIELDS_N_COMP);

  mparticles_vpic_t mprts = mprts_base->get_as<mparticles_vpic_t>();

  // Divergence clean e
  if (clean_div_e) {
    // needs E, rhof, rhob, material
    psc_marder_vpic_clean_div_e(marder, mf.mflds(), mprts.mprts());
    // upates E, rhof, div_e_err
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);

  // Divergence clean b
  if (clean_div_b) {
    // needs B
    psc_marder_vpic_clean_div_b(marder, mf.mflds());
    // updates B, div_b_err
  }
  
  // Synchronize the shared faces
  if (sync_shared) {
    // needs E, B, TCA
    mpi_printf(psc_marder_comm(marder), "Synchronizing shared tang e, norm b\n");
    double err = psc_mfields_synchronize_tang_e_norm_b(mf.mflds());
    mpi_printf(psc_marder_comm(marder), "Domain desynchronization error = %e (arb units)\n", err);
    // updates E, B, TCA
  }

  mf.put_as(mflds_base, EX, 16);
}

// ----------------------------------------------------------------------
// psc_marder: subclass "vpic"

struct psc_marder_ops psc_marder_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_marder_vpic),
  .run                   = psc_marder_vpic_run,
};

