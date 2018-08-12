
#pragma once

#include <psc_method.h>
#include <mrc_profile.h>
#include <psc_method.h>
#include <psc_diag.h>

#include <particles.hxx>

#include <push_particles.hxx>
#include <checks.hxx>
#include <output_particles.hxx>
#include <output_fields_c.hxx>

#ifdef VPIC
#include <psc_fields_vpic.h>
#include <psc_particles_vpic.h>

void psc_method_vpic_inc_step(struct psc_method *method, int timestep); // FIXME
void psc_method_vpic_print_status(struct psc_method *method); // FIXME
#endif

// ======================================================================
// PscParams

struct PscParams
{
  double cfl = { .75 };            // CFL number used to determine time step
  int nmax;                        // Number of timesteps to run
  double wallclock_limit = { 0. }; // Maximum wallclock time to run
  bool write_checkpoint = { false };
  int write_checkpoint_every_step = { 0 };

  bool detailed_profiling; // output profiling info for each process separately
  int stats_every;         // output timing and other info every so many steps
};
  
// ======================================================================
// Psc

template<typename PscConfig>
struct Psc
{
  using Mparticles_t = typename PscConfig::Mparticles_t;
  using MfieldsState = typename PscConfig::MfieldsState;
  using Balance_t = typename PscConfig::Balance_t;
  using Sort_t = typename PscConfig::Sort_t;
  using Collision_t = typename PscConfig::Collision_t;
  using PushParticles_t = typename PscConfig::PushParticles_t;
  using PushFields_t = typename PscConfig::PushFields_t;
  using Bnd_t = typename PscConfig::Bnd_t;
  using BndFields_t = typename PscConfig::BndFields_t;
  using BndParticles_t = typename PscConfig::BndParticles_t;
  using Checks_t = typename PscConfig::Checks_t;
  using Marder_t = typename PscConfig::Marder_t;
  using Simulation = typename PscConfig::Simulation;

  // ----------------------------------------------------------------------
  // ctor

  Psc()
  {
    time_start_ = MPI_Wtime();
    
    psc_ = psc_create(MPI_COMM_WORLD);
    psc_set_from_options(psc_);

#ifdef VPIC
    vgrid_ = Grid::create();
#endif
  }
  

  // ----------------------------------------------------------------------
  // init

  void init()
  {
#ifdef VPIC
    mflds_.reset(new MfieldsState{grid(), *material_list_});
#else
    mflds_.reset(new MfieldsState{grid()});
#endif
      
#ifdef VPIC
    hydro_.reset(new MfieldsHydroVpic{grid(), 16, psc_->ibn});
    interpolator_.reset(new Interpolator{vgrid_});
    accumulator_.reset(new Accumulator{vgrid_});
#endif
    mprts_.reset(new Mparticles_t{grid()});
    sort_.reset(new Sort_t{});
    pushp_.reset(new PushParticles_t{});
    pushf_.reset(new PushFields_t{});
    bnd_.reset(new Bnd_t{psc_->grid(), psc_->mrc_domain_, psc_->ibn});
    bndf_.reset(new BndFields_t{});
    bndp_.reset(new BndParticles_t{psc_->mrc_domain_, psc_->grid()});
  }

  // ----------------------------------------------------------------------
  // dtor

  ~Psc()
  {
    psc_destroy(psc_);
    delete sim_;
  }
  
  // ----------------------------------------------------------------------
  // initialize_stats
  
  void initialize_stats()
  {
    st_nr_particles = psc_stats_register("nr particles");
    st_time_step = psc_stats_register("time entire step");

    // generic stats categories
    st_time_particle = psc_stats_register("time particle update");
    st_time_field = psc_stats_register("time field update");
    st_time_comm = psc_stats_register("time communication");
    st_time_output = psc_stats_register("time output");
  }
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    psc_view(psc_);
    mprts_->view();

    if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
#ifdef VPIC
      initialize_vpic();
#endif
    } else {
      initialize_default(psc_->method, psc_, *mflds_, *mprts_);
    }

    // initial output / stats
    mpi_printf(psc_comm(psc_), "Performing initial diagnostics.\n");
    diagnostics();
    print_status();

    mpi_printf(psc_comm(psc_), "Initialization complete.\n");
  }

  // ----------------------------------------------------------------------
  // integrate

  void integrate()
  {
    static int pr;
    if (!pr) {
      pr = prof_register("psc_step", 1., 0, 0);
    }

    mpi_printf(psc_comm(psc_), "*** Advancing\n");
    double elapsed = MPI_Wtime();

    bool first_iteration = true;
    while (psc_->timestep < p_.nmax) {
      prof_start(pr);
      psc_stats_start(st_time_step);

      if (!first_iteration &&
	  p_.write_checkpoint_every_step > 0 &&
	  psc_->timestep % p_.write_checkpoint_every_step == 0) {
	psc_write_checkpoint(psc_);
      }
      first_iteration = false;

      mpi_printf(psc_comm(psc_), "**** Step %d / %d, Code Time %g, Wall Time %g\n", psc_->timestep + 1,
		 p_.nmax, psc_->timestep * dt(), MPI_Wtime() - time_start_);

      prof_start(pr_time_step_no_comm);
      prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
    
      psc_->timestep++; // FIXME, too hacky
#ifdef VPIC
      if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
	psc_method_vpic_inc_step(psc_->method, psc_->timestep);
      }
#endif
      
      diagnostics();

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts_->get_n_prts();

      if (psc_->timestep % p_.stats_every == 0) {
	print_status();
      }

      if (p_.wallclock_limit > 0.) {
	double wallclock_elapsed = MPI_Wtime() - time_start_;
	double wallclock_elapsed_max;
	MPI_Allreduce(&wallclock_elapsed, &wallclock_elapsed_max, 1, MPI_DOUBLE, MPI_MAX,
		      MPI_COMM_WORLD);
      
	if (wallclock_elapsed_max > p_.wallclock_limit) {
	  mpi_printf(MPI_COMM_WORLD, "WARNING: Max wallclock time elapsed!\n");
	  break;
	}
      }
    }

    if (p_.write_checkpoint) {
      psc_write_checkpoint(psc_);
    }

    // FIXME, merge with existing handling of wallclock time
    elapsed = MPI_Wtime() - elapsed;

    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    mpi_printf(psc_comm(psc_), "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
	       elapsed, w, d, h, m, s );
  }

  virtual void step() = 0;

  // ----------------------------------------------------------------------
  // set_dt
  
  static double set_dt(const PscParams& p, const Grid_t::Domain& domain)
  {
    double inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!domain.isInvar(d)) {
	inv_sum += 1. / sqr(domain.dx[d]);
      }
    }
    if (!inv_sum) { // simulation has 0 dimensions
      inv_sum = 1.;
    }
    return p.cfl * sqrt(1./inv_sum);
  }
  
  // ----------------------------------------------------------------------
  // setup_diagnostics

  void setup_diagnostics()
  {
#ifdef VPIC
    sim_->setupDiag();
#endif
  }

protected:
  double dt() const { return psc_->grid().dt; }

private:

  // ----------------------------------------------------------------------
  // print_profiling

  void print_profiling()
  {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!p_.detailed_profiling) {
      prof_print_mpi(MPI_COMM_WORLD);
    } else {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      for (int i = 0; i < size; i++) {
	if (i == rank) {
	  mprintf("profile\n");
	  prof_print();
	}
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
  

  // ----------------------------------------------------------------------
  // initialize_default
  
  static void initialize_default(struct psc_method *method, struct psc *psc,
				 MfieldsState& mflds, Mparticles_t& mprts)
  {
    //pushp_.stagger(mprts, mflds); FIXME, vpic does it
  }

#ifdef VPIC

  // ----------------------------------------------------------------------
  // initialize_vpic
  
  void initialize_vpic()
  {
    // FIXME, just change the uses
    auto sim = sim_;
    auto psc = psc_;
    auto& mprts = *mprts_;
    auto& mflds = mflds_;
    // Do some consistency checks on user initialized fields
    
    mpi_printf(psc_comm(psc), "Checking interdomain synchronization\n");
    double err;
    TIC err = CleanDivOps::synchronize_tang_e_norm_b(mflds->vmflds()); TOC(synchronize_tang_e_norm_b, 1);
    mpi_printf(psc_comm(psc), "Error = %g (arb units)\n", err);
    
    mpi_printf(psc_comm(psc), "Checking magnetic field divergence\n");
    TIC CleanDivOps::compute_div_b_err(mflds->vmflds()); TOC(compute_div_b_err, 1);
    TIC err = CleanDivOps::compute_rms_div_b_err(mflds->vmflds()); TOC(compute_rms_div_b_err, 1);
    mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
    TIC CleanDivOps::clean_div_b(mflds->vmflds()); TOC(clean_div_b, 1);
    
    // Load fields not initialized by the user
    
    mpi_printf(psc_comm(psc), "Initializing radiation damping fields\n");
    TIC AccumulateOps::compute_curl_b(mflds->vmflds()); TOC(compute_curl_b, 1);
    
    mpi_printf(psc_comm(psc), "Initializing bound charge density\n");
    TIC CleanDivOps::clear_rhof(mflds->vmflds()); TOC(clear_rhof, 1);
    ParticlesOps::accumulate_rho_p(mprts.vmprts_, mflds->vmflds());
    CleanDivOps::synchronize_rho(mflds->vmflds());
    TIC AccumulateOps::compute_rhob(mflds->vmflds()); TOC(compute_rhob, 1);
    
    // Internal sanity checks
    
    mpi_printf(psc_comm(psc), "Checking electric field divergence\n");
    TIC CleanDivOps::compute_div_e_err(mflds->vmflds()); TOC(compute_div_e_err, 1);
    TIC err = CleanDivOps::compute_rms_div_e_err(mflds->vmflds()); TOC(compute_rms_div_e_err, 1);
    mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
    TIC CleanDivOps::clean_div_e(mflds->vmflds()); TOC(clean_div_e, 1);
    
    mpi_printf(psc_comm(psc), "Rechecking interdomain synchronization\n");
    TIC err = CleanDivOps::synchronize_tang_e_norm_b(mflds->vmflds()); TOC(synchronize_tang_e_norm_b, 1);
    mpi_printf(psc_comm(psc), "Error = %e (arb units)\n", err);
    
    mpi_printf(psc_comm(psc), "Uncentering particles\n");
    auto& vmprts = mprts.vmprts_;
    if (!vmprts.empty()) {
      TIC InterpolatorOps::load(*interpolator_, mflds->vmflds()); TOC(load_interpolator, 1);
      
      for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
	TIC ParticlesOps::uncenter_p(&*sp, *interpolator_); TOC(uncenter_p, 1);
      }
    }
  }

#endif
  
  // ----------------------------------------------------------------------
  // diagnostics

  virtual void diagnostics()
  {
#ifdef VPIC
    if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
      sim_->runDiag(mprts_->vmprts_, mflds_->vmflds(), *interpolator_, *hydro_->vmflds_hydro, grid().domain.np);
    }
#else
    // FIXME
    psc_diag_run(psc_->diag, psc_, *mprts_, *mflds_);
    // FIXME
    (*outf_)(*mflds_, *mprts_);
#endif
    PscOutputParticlesBase{psc_->output_particles}.run(*mprts_);
  }

  // ----------------------------------------------------------------------
  // print_status

  void print_status()
  {
#ifdef VPIC
    if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
      psc_method_vpic_print_status(psc_->method);
    }
#endif
    psc_stats_log(psc_->timestep);
    print_profiling();
  }

public:
  const Grid_t& grid() { return *grid_; }

protected:
  double time_start_;

  PscParams p_;
  const Grid_t* grid_;
  Simulation* sim_;
#ifdef VPIC
  Grid* vgrid_;
#endif
  psc* psc_;

#ifdef VPIC
  // FIXME? MaterialList would be perfectly fine as actual member
  std::unique_ptr<MaterialList> material_list_;
#endif
  std::unique_ptr<MfieldsState> mflds_;
#ifdef VPIC
  std::unique_ptr<MfieldsHydroVpic> hydro_;
  std::unique_ptr<Interpolator> interpolator_;
  std::unique_ptr<Accumulator> accumulator_;
  ParticleBcList particle_bc_list_;
#endif
  std::unique_ptr<Mparticles_t> mprts_;

  std::unique_ptr<Balance_t> balance_;
  std::unique_ptr<Sort_t> sort_;
  std::unique_ptr<Collision_t> collision_;
  std::unique_ptr<PushParticles_t> pushp_;
  std::unique_ptr<PushFields_t> pushf_;
  std::unique_ptr<Bnd_t> bnd_;
  std::unique_ptr<BndFields_t> bndf_;
  std::unique_ptr<BndParticles_t> bndp_;
  std::unique_ptr<Checks_t> checks_;
  std::unique_ptr<Marder_t> marder_;
  std::unique_ptr<OutputFieldsC> outf_;

  // FIXME, maybe should be private
  // need to make sure derived class sets these (? -- or just leave them off by default)
  int balance_interval;
  int sort_interval;
  int marder_interval;

  int num_comm_round = {3};
  
  int st_nr_particles;
  int st_time_step;
};
