
#pragma once

#include "marder.hxx"

#ifdef USE_VPIC

// ======================================================================
// MarderVpicOpsWrap

template<typename _Mparticles, typename _MfieldsState>
struct MarderVpicOpsWrap
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  
  void clear_rhof(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    TIC fa->kernel->clear_rhof(fa); TOC(clear_rhof, 1);
  }

  void synchronize_rho(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    TIC fa->kernel->synchronize_rho(fa); TOC(synchronize_rho, 1);
  }

  void compute_div_e_err(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    TIC fa->kernel->compute_div_e_err(fa); TOC(compute_div_e_err, 1);
  }

  double compute_rms_div_e_err(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    double err;
    TIC err = fa->kernel->compute_rms_div_e_err(fa); TOC(compute_rms_div_e_err, 1);
    return err;
  }

  void clean_div_e(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    TIC fa->kernel->clean_div_e(fa); TOC(clean_div_e, 1);
  }
  
  void compute_div_b_err(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    TIC fa->kernel->compute_div_b_err(fa); TOC(compute_div_e_err, 1);
  }

  double compute_rms_div_b_err(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    double err;
    TIC err = fa->kernel->compute_rms_div_b_err(mflds); TOC(compute_rms_div_e_err, 1);
    return err;
  }

  void clean_div_b(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    TIC fa->kernel->clean_div_b(fa); TOC(clean_div_e, 1);
  }

  double synchronize_tang_e_norm_b(MfieldsState& mflds)
  {
    field_array_t* fa = mflds;
    double err;
    TIC err = fa->kernel->synchronize_tang_e_norm_b(fa); TOC(synchronize_tang_e_norm_b, 1);
    return err;
  }

  void accumulate_rho_p(Mparticles& mprts, MfieldsState& mflds)
  {
    auto& vmprts = mprts.vmprts();
    field_array_t* fa = mflds;
    for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
      TIC ::accumulate_rho_p(fa, &*sp); TOC(accumulate_rho_p, 1);
    }
  }
};

#endif

// ======================================================================
// MarderVpicOps

template<typename _Mparticles, typename _MfieldsState>
struct MarderVpicOps
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using LocalOps = PscFieldArrayLocalOps<MfieldsState>;
  using RemoteOps = PscFieldArrayRemoteOps<MfieldsState>;
  using CleanDivOps = PscCleanDivOps<Mparticles, MfieldsState, LocalOps, RemoteOps>;
  
  void clear_rhof(MfieldsState& mflds)
  {
    TIC CleanDivOps::clear_rhof(mflds); TOC(clear_rhof, 1);
  }

  void synchronize_rho(MfieldsState& mflds)
  {
    TIC CleanDivOps::synchronize_rho(mflds); TOC(synchronize_rho, 1);
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

  double synchronize_tang_e_norm_b(MfieldsState& mflds)
  {
    double err;
    TIC err = CleanDivOps::synchronize_tang_e_norm_b(mflds); TOC(synchronize_tang_e_norm_b, 1);
    return err;
  }

  void accumulate_rho_p(Mparticles& mprts, MfieldsState& mflds)
  {
    CleanDivOps::accumulate_rho_p(mprts, mflds);
  }
};

// ======================================================================
// MarderVpic_

template<typename Ops>
struct MarderVpic_ : MarderBase, Ops
{
  using Mparticles = typename Ops::Mparticles;
  using MfieldsState = typename Ops::MfieldsState;
  using real_t = typename MfieldsState::real_t;
  
  MarderVpic_(const Grid_t& grid, real_t diffusion, int loop, bool dump)
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
      double err = synchronize_tang_e_norm_b(mflds);
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

template<typename Mparticles, typename MfieldsState>
using MarderVpicWrap = MarderVpic_<MarderVpicOpsWrap<Mparticles, MfieldsState>>;

template<typename Mparticles, typename MfieldsState>
using MarderVpic = MarderVpic_<MarderVpicOps<Mparticles, MfieldsState>>;

