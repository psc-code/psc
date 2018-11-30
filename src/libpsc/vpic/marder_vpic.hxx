
#pragma once

#include "marder.hxx"

// ======================================================================
// MarderVpicOps

template<typename Mparticles, typename MfieldsState, typename CleanDivOps>
struct MarderVpicOps
{
  void clear_rhof(MfieldsState& mflds)
  {
    TIC CleanDivOps::clear_rhof(mflds); TOC(clear_rhof, 1);
  }

  void accumulate_rho_p(Mparticles& mprts, MfieldsState& mflds)
  {
    CleanDivOps::accumulate_rho_p(mprts, mflds);
  }

  void synchronize_rho(MfieldsState& mflds)
  {
    CleanDivOps::synchronize_rho(mflds);
  }

  void compute_div_e_err(MfieldsState& mflds)
  {
    TIC CleanDivOps::compute_div_e_err(mflds); TOC(compute_div_e_err, 1);
  }

  double compute_rms_div_e_err(MfieldsState& mflds)
  {
    double err;
    TIC err = CleanDivOps::compute_rms_div_e_err(mflds); TOC(compute_rms_div_e_err, 1);
    return err;
  }

  void clean_div_e(MfieldsState& mflds)
  {
    TIC CleanDivOps::clean_div_e(mflds); TOC(clean_div_e, 1);
  }
  
  void compute_div_b_err(MfieldsState& mflds)
  {
    TIC CleanDivOps::compute_div_b_err(mflds); TOC(compute_div_e_err, 1);
  }

  double compute_rms_div_b_err(MfieldsState& mflds)
  {
    double err;
    TIC err = CleanDivOps::compute_rms_div_b_err(mflds); TOC(compute_rms_div_e_err, 1);
    return err;
  }

  void clean_div_b(MfieldsState& mflds)
  {
    TIC CleanDivOps::clean_div_b(mflds); TOC(clean_div_e, 1);
  }
};

// ======================================================================
// MarderVpic

template<typename Mparticles, typename _CleanDivOps>
struct MarderVpic : MarderBase, MarderVpicOps<Mparticles, typename _CleanDivOps::MfieldsState, _CleanDivOps>
{
  using CleanDivOps = _CleanDivOps;
  using MfieldsState = typename CleanDivOps::MfieldsState;
  using real_t = typename MfieldsState::real_t;
  
  MarderVpic(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : comm_{grid.comm()}
  {}
  
  // ----------------------------------------------------------------------
  // psc_marder_vpic_clean_div_e

  void psc_marder_vpic_clean_div_e(MfieldsState& mflds, Mparticles& mprts)
  {
    mpi_printf(comm_, "Divergence cleaning electric field\n");

    this->clear_rhof(mflds);
    this->accumulate_rho_p(mprts, mflds);
    this->synchronize_rho(mflds);

    for (int round = 0; round < num_div_e_round_; round++ ) {
      this->compute_div_e_err(mflds);
      if (round == 0 || round == num_div_e_round_ - 1) {
	double err = this->compute_rms_div_e_err(mflds);
	mpi_printf(comm_, "%s rms error = %e (charge/volume)\n",
		   round == 0 ? "Initial" : "Cleaned", err);
      }
      this->clean_div_e(mflds);
    }
  }

  // ----------------------------------------------------------------------
  // psc_marder_vpic_clean_div_b

  void psc_marder_vpic_clean_div_b(MfieldsState& mflds)
  {
    mpi_printf(comm_, "Divergence cleaning magnetic field\n");
  
    for (int round = 0; round < num_div_b_round_; round++) {
      // FIXME, having the TIC/TOC stuff with every use is not pretty
      this->compute_div_b_err(mflds);
      if (round == 0 || round == num_div_b_round_ - 1) {
	double err = this->compute_rms_div_b_err(mflds);
	mpi_printf(comm_, "%s rms error = %e (charge/volume)\n",
		   round == 0 ? "Initial" : "Cleaned", err);
      }
      this->clean_div_b(mflds);
    }
  }

  // ----------------------------------------------------------------------
  // run

  void run(MfieldsStateBase& mflds_base, MparticlesBase& mprts_base) override
  {
    assert(0);
  }
  
  void operator()(MfieldsState& mflds, Mparticles& mprts)
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

