
#pragma once

#include "marder.hxx"

// ======================================================================
// psc_marder "vpic"

struct MarderVpic : MarderBase
{
  using real_t = MfieldsState::real_t;
  
  MarderVpic(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : comm_{grid.comm()}
  {}
  
  // ----------------------------------------------------------------------
  // psc_marder_vpic_clean_div_e

  void psc_marder_vpic_clean_div_e(MfieldsState& mflds, MparticlesVpic& mprts)
  {
    mpi_printf(comm_, "Divergence cleaning electric field\n");
  
    TIC CleanDivOps::clear_rhof(mflds); TOC(clear_rhof, 1);
    ParticlesOps::accumulate_rho_p(mprts, mflds);
    CleanDivOps::synchronize_rho(mflds);

    for (int round = 0; round < num_div_e_round_; round++ ) {
      TIC CleanDivOps::compute_div_e_err(mflds); TOC(compute_div_e_err, 1);
      if (round == 0 || round == num_div_e_round_ - 1) {
	double err;
	TIC err = CleanDivOps::compute_rms_div_e_err(mflds); TOC(compute_rms_div_e_err, 1);
	mpi_printf(comm_, "%s rms error = %e (charge/volume)\n",
		   round == 0 ? "Initial" : "Cleaned", err);
      }
      TIC CleanDivOps::clean_div_e(mflds); TOC(clean_div_e, 1);
    }
  }

  // ----------------------------------------------------------------------
  // psc_marder_vpic_clean_div_b

  void psc_marder_vpic_clean_div_b(MfieldsState& mflds)
  {
    mpi_printf(comm_, "Divergence cleaning magnetic field\n");
  
    for (int round = 0; round < num_div_b_round_; round++) {
      // FIXME, having the TIC/TOC stuff with every use is not pretty
      TIC CleanDivOps::compute_div_b_err(mflds); TOC(compute_div_b_err, 1);
      if (round == 0 || round == num_div_b_round_ - 1) {
	double err;
	TIC err = CleanDivOps::compute_rms_div_b_err(mflds); TOC(compute_rms_div_b_err, 1);
	mpi_printf(comm_, "%s rms error = %e (charge/volume)\n",
		   round == 0 ? "Initial" : "Cleaned", err);
      }
      TIC CleanDivOps::clean_div_b(mflds); TOC(clean_div_b, 1);
    }
  }

  // ----------------------------------------------------------------------
  // run

  void run(MfieldsStateBase& mflds_base, MparticlesBase& mprts_base) override
  {
    assert(0);
  }
  
  void operator()(MfieldsState& mflds, MparticlesVpic& mprts)
  {
    const auto& grid = mflds.grid();
    bool clean_div_e = (clean_div_e_interval_ > 0 && grid.timestep() % clean_div_e_interval_ == 0);
    bool clean_div_b = (clean_div_b_interval_ > 0 && grid.timestep() % clean_div_b_interval_ == 0);
    bool sync_shared = (sync_shared_interval_ > 0 && grid.timestep() % sync_shared_interval_ == 0);

    if (!(clean_div_e || clean_div_b || sync_shared)) {
      return;
    }
  
    // Divergence clean e
    if (clean_div_e) {
      // needs E, rhof, rhob, material
      psc_marder_vpic_clean_div_e(mflds, mprts);
      // upates E, rhof, div_e_err
    }

    // Divergence clean b
    if (clean_div_b) {
      // needs B
      psc_marder_vpic_clean_div_b(mflds);
      // updates B, div_b_err
    }
  
    // Synchronize the shared faces
    if (sync_shared) {
      // needs E, B, TCA
      mpi_printf(comm_, "Synchronizing shared tang e, norm b\n");
      double err;
      TIC err = CleanDivOps::synchronize_tang_e_norm_b(mflds); TOC(synchronize_tang_e_norm_b, 1);
      mpi_printf(comm_, "Domain desynchronization error = %e (arb units)\n", err);
      // updates E, B, TCA
    }
  }

private:
  MPI_Comm comm_;
  int clean_div_e_interval_ = 10; // FIXME, hardcoded...
  int clean_div_b_interval_ = 10;
  int sync_shared_interval_ = 10;
  int num_div_e_round_ = 3;
  int num_div_b_round_ = 3;
};

