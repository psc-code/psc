
#include <psc_config.h>

#include <psc.h>

#include <mrc_profile.h>
#include <psc_diag.h>

#include <particles.hxx>

#include <push_particles.hxx>
#include <checks.hxx>
#include <output_particles.hxx>
#include <output_fields_c.hxx>

// ======================================================================
// PscParams

struct PscParams
{
  double cfl = .75;            // CFL number used to determine time step
  int nmax;                    // Number of timesteps to run
  double wallclock_limit = 0.; // Maximum wallclock time to run
  bool write_checkpoint = false;
  int write_checkpoint_every_step = 0;

  bool detailed_profiling = false; // output profiling info for each process separately
  int stats_every = 10;    // output timing and other info every so many steps
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

  // ----------------------------------------------------------------------
  // ctor

  Psc()
    : grid_{ggrid}
  {
    time_start_ = MPI_Wtime();

    // FIXME, we should use RngPool consistently throughout
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srandom(rank);

    diag_ = psc_diag_create(MPI_COMM_WORLD);
    psc_diag_set_from_options(diag_);

    outp_ = psc_output_particles_create(MPI_COMM_WORLD);
    psc_output_particles_set_from_options(outp_);
  }

  // ----------------------------------------------------------------------
  // define_grid

  void define_grid(Grid_t::Domain& domain, GridBc& bc, Grid_t::Kinds& kinds,
		   double dt, Grid_t::NormalizationParams& norm_params)
  {
    auto coeff = Grid_t::Normalization{norm_params};
    grid_ = Grid_t::psc_make_grid(domain, bc, kinds, coeff, dt, ibn);
  }
  
  // ----------------------------------------------------------------------
  // define_field_array

  void define_field_array(double damp = 0.)
  {
    mflds_.reset(new MfieldsState{grid()});
  }
  
  // ----------------------------------------------------------------------
  // init

  void init()
  {
    sort_.reset(new Sort_t{});
    pushp_.reset(new PushParticles_t{});
    pushf_.reset(new PushFields_t{});
    bnd_.reset(new Bnd_t{grid(), ibn});
    bndf_.reset(new BndFields_t{});
    bndp_.reset(new BndParticles_t{grid()});

    psc_diag_setup(diag_);
    psc_output_particles_setup(outp_);

    initialize_stats();
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

    // FIXME not quite the right place
    pr_time_step_no_comm = prof_register("time step w/o comm", 1., 0, 0);
  }
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    initialize_default();

    // initial output / stats
    mpi_printf(grid().comm(), "Performing initial diagnostics.\n");
    diagnostics();
    print_status();

    mpi_printf(grid().comm(), "Initialization complete.\n");
  }

  // ----------------------------------------------------------------------
  // integrate

  void integrate()
  {
    static int pr;
    if (!pr) {
      pr = prof_register("psc_step", 1., 0, 0);
    }

    mpi_printf(grid().comm(), "*** Advancing\n");
    double elapsed = MPI_Wtime();

    bool first_iteration = true;
    while (grid().timestep() < p_.nmax) {
      prof_start(pr);
      psc_stats_start(st_time_step);

      if (!first_iteration &&
	  p_.write_checkpoint_every_step > 0 &&
	  grid().timestep() % p_.write_checkpoint_every_step == 0) {
	//psc_write_checkpoint(psc_);
      }
      first_iteration = false;

      mpi_printf(grid().comm(), "**** Step %d / %d, Code Time %g, Wall Time %g\n", grid().timestep() + 1,
		 p_.nmax, grid().timestep() * grid().dt, MPI_Wtime() - time_start_);

      prof_start(pr_time_step_no_comm);
      prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
      grid_->timestep_++; // FIXME, too hacky
      
      diagnostics();

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts_->get_n_prts();

      if (grid().timestep() % p_.stats_every == 0) {
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
      //psc_write_checkpoint(psc_);
    }

    // FIXME, merge with existing handling of wallclock time
    elapsed = MPI_Wtime() - elapsed;

    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    mpi_printf(grid().comm(), "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
	       elapsed, w, d, h, m, s );
  }

  // ----------------------------------------------------------------------
  // step_psc

  void step_psc()
  {
    using DIM = typename PscConfig::dim_t;

    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject_prts;
    if (!pr_sort) {
      pr_sort = prof_register("step_sort", 1., 0, 0);
      pr_collision = prof_register("step_collision", 1., 0, 0);
      pr_push_prts = prof_register("step_push_prts", 1., 0, 0);
      pr_inject_prts = prof_register("step_inject_prts", 1., 0, 0);
      pr_push_flds = prof_register("step_push_flds", 1., 0, 0);
      pr_bndp = prof_register("step_bnd_prts", 1., 0, 0);
      pr_bndf = prof_register("step_bnd_flds", 1., 0, 0);
      pr_checks = prof_register("step_checks", 1., 0, 0);
      pr_marder = prof_register("step_marder", 1., 0, 0);
    }

    // state is at: x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}
    MPI_Comm comm = grid().comm();
    int timestep = grid().timestep();

    auto& mprts = *mprts_;
    auto& mflds = *mflds_;

    if (balance_interval > 0 && timestep % balance_interval == 0) {
      (*balance_)(mprts);
    }

    if (sort_interval > 0 && timestep % sort_interval == 0) {
      mpi_printf(comm, "***** Sorting...\n");
      prof_start(pr_sort);
      (*sort_)(mprts);
      prof_stop(pr_sort);
    }
    
    if (collision_->interval() > 0 && timestep % collision_->interval() == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      prof_start(pr_collision);
      (*collision_)(mprts);
      prof_stop(pr_collision);
    }
    
    // === particle injection
    prof_start(pr_inject_prts);
    inject_particles();
    prof_stop(pr_inject_prts);

    if (checks_->continuity_every_step > 0 && timestep % checks_->continuity_every_step == 0) {
      mpi_printf(comm, "***** Checking continuity...\n");
      prof_start(pr_checks);
      checks_->continuity_before_particle_push(mprts);
      prof_stop(pr_checks);
    }

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    prof_start(pr_push_prts);
    pushp_->push_mprts(mprts, mflds);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

    // === field propagation B^{n+1/2} -> B^{n+1}
    prof_start(pr_push_flds);
    pushf_->push_H(mflds, .5, DIM{});
    prof_stop(pr_push_flds);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    prof_start(pr_bndp);
    (*bndp_)(mprts);
    prof_stop(pr_bndp);

    // === field propagation E^{n+1/2} -> E^{n+3/2}
    prof_start(pr_bndf);
#if 1
    bndf_->fill_ghosts_H(mflds);
    bnd_->fill_ghosts(mflds, HX, HX + 3);
#endif
    
    bndf_->add_ghosts_J(mflds);
    bnd_->add_ghosts(mflds, JXI, JXI + 3);
    bnd_->fill_ghosts(mflds, JXI, JXI + 3);
    prof_stop(pr_bndf);
    
    prof_restart(pr_push_flds);
    pushf_->push_E(mflds, 1., DIM{});
    prof_stop(pr_push_flds);
    
#if 1
    prof_restart(pr_bndf);
    bndf_->fill_ghosts_E(mflds);
    bnd_->fill_ghosts(mflds, EX, EX + 3);
    prof_stop(pr_bndf);
#endif
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}
      
    // === field propagation B^{n+1} -> B^{n+3/2}
    prof_restart(pr_push_flds);
    pushf_->push_H(mflds, .5, DIM{});
    prof_stop(pr_push_flds);

#if 1
    prof_start(pr_bndf);
    bndf_->fill_ghosts_H(mflds);
    bnd_->fill_ghosts(mflds, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}
#endif

    if (checks_->continuity_every_step > 0 && timestep % checks_->continuity_every_step == 0) {
      prof_restart(pr_checks);
      checks_->continuity_after_particle_push(mprts, mflds);
      prof_stop(pr_checks);
    }
    
    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not its natural time,
    // but div B should be == 0 at any time...)
    if (marder_interval > 0 && timestep % marder_interval == 0) {
      mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      (*marder_)(mflds, mprts);
      prof_stop(pr_marder);
    }
    
    if (checks_->gauss_every_step > 0 && timestep % checks_->gauss_every_step == 0) {
      prof_restart(pr_checks);
      checks_->gauss(mprts, mflds);
      prof_stop(pr_checks);
    }
    
    //psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
  }

  virtual void step()
  {
    step_psc();
  }

  // ----------------------------------------------------------------------
  // inject_particles

  virtual void inject_particles()
  {}
  
  // ----------------------------------------------------------------------
  // define_periodic_grid
  
  void define_periodic_grid(double xl[3], double xh[3], const int gdims[3], const int np[3])
  {
  }
  
  // ----------------------------------------------------------------------
  // set_domain_field_bc
  
  void set_domain_field_bc(Int3 bnd, int bc)
  {
  }

  // ----------------------------------------------------------------------
  // set_domain_particle_bc
  
  void set_domain_particle_bc(Int3 bnd, int bc)
  {
  }

  void grid_setup_communication()
  {
  }

  // ----------------------------------------------------------------------
  // courant length
  
  double courant_length(const Grid_t::Domain& domain)
  {
    double inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!domain.isInvar(d)) {
	inv_sum += 1. / sqr(domain.dx[d]);
      }
    }
    if (!inv_sum) { // simulation has 0 dimensions (happens in some test?)
      inv_sum = 1.;
    }
    return sqrt(1. / inv_sum);
  }

  // ----------------------------------------------------------------------
  // create_diagnotics

  void create_diagnostics(int interval)
  {
  }
  
  // ----------------------------------------------------------------------
  // setup_diagnostics

  void setup_diagnostics()
  {
  }

  // ----------------------------------------------------------------------
  // run_diagnostics
  
  void run_diagnostics()
  {
  }

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
  
  void initialize_default()
  {
    //pushp_.stagger(mprts, mflds); FIXME, vpic does it

    checks_->gauss(*mprts_, *mflds_);
  }

  // ----------------------------------------------------------------------
  // diagnostics

  virtual void diagnostics()
  {
    // FIXME
    psc_diag_run(diag_, *mprts_, *mflds_);
    // FIXME
    (*outf_)(*mflds_, *mprts_);
    PscOutputParticlesBase{outp_}.run(*mprts_);
  }

  // ----------------------------------------------------------------------
  // print_status

  void print_status()
  {
    psc_stats_log(grid().timestep());
    print_profiling();
  }

public:
  const Grid_t& grid() { return *grid_; }

protected:
  double time_start_;

  PscParams p_;
  Grid_t*& grid_;

  std::unique_ptr<MfieldsState> mflds_;
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

  psc_diag* diag_;             ///< timeseries diagnostics
  psc_output_particles* outp_; ///< particle output

  // FIXME, maybe should be private
  // need to make sure derived class sets these (? -- or just leave them off by default)
  int balance_interval;
  int sort_interval;
  int marder_interval;

  int num_comm_round = {3};
  Int3 ibn = {2, 2, 2}; // FIXME!!! need to factor in invar dims (but not in vpic...)
  
  int st_nr_particles;
  int st_time_step;
};

#include <balance.hxx>
#include <particles.hxx>
#include <fields3d.hxx>
#include <push_particles.hxx>
#include <push_fields.hxx>
#include <sort.hxx>
#include <collision.hxx>
#include <bnd_particles.hxx>
#include <bnd.hxx>
#include <bnd_fields.hxx>
#include <setup_particles.hxx>
#include <setup_fields.hxx>

#include "../libpsc/psc_checks/checks_impl.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#include "../libpsc/cuda/setup_particles_cuda.hxx"
#endif

#include "psc_config.hxx"

enum {
  MY_ION,
  MY_ELECTRON,
  N_MY_KINDS,
};

// EDIT to change order / floating point type / cuda / 2d/3d
using dim_t = dim_xyz;
using PscConfig = PscConfig1vbecSingle<dim_t>;
//using PscConfig = PscConfig1vbecCuda<dim_t>;

// ======================================================================
// PscTestIo

struct PscTestIo : Psc<PscConfig>
{
  using DIM = PscConfig::dim_t;

  // ----------------------------------------------------------------------
  // ctor
  
  PscTestIo()
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    p_.nmax = 5001;
    p_.cfl = 0.75;

    BB_ = 0.;
    Zi_ = 1.;
    
    // --- for background plasma
    background_n_  = .002;
    background_Te_ = .001;
    background_Ti_ = .001;
    
    // -- setup particle kinds
    // last population ("e") is neutralizing
    // FIXME, hardcoded mass ratio 100
    Grid_t::Kinds kinds = {{Zi_, 100.*Zi_, "i"}, { -1., 1., "e"}};
    
    double d_i = sqrt(kinds[MY_ION].m / kinds[MY_ION].q);
    
    mpi_printf(comm, "d_e = %g, d_i = %g\n", 1., d_i);
    mpi_printf(comm, "lambda_De (background) = %g\n", sqrt(background_Te_));
    
    // --- setup domain
    Grid_t::Real3 LL = { 400., 800., 400.*6 }; // domain size (in d_e)
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 40, 80, 4 }; // division into patches

    if (dim::InvarX::value) { ibn[0] = 0; } // FIXME, wrong place, not for VPIC...
    
    auto grid_domain = Grid_t::Domain{gdims, LL, -.5 * LL, np};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 5;

    double dt = p_.cfl * courant_length(grid_domain);
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles_t{grid()});

    // -- Balance
    balance_interval = 50;
    balance_.reset(new Balance_t{balance_interval, .1, true});

    // -- Checks
    ChecksParams checks_params{};
    checks_params.continuity_every_step = 20;
    checks_params.continuity_threshold = 1e-5;
    checks_.reset(new Checks_t{grid(), comm, checks_params});

    // -- output fields
    OutputFieldsCParams outf_params;
    outf_params.output_fields = "e,h,j,n_1st_single,v_1st_single,T_1st_single";
    outf_params.pfield_step = 100;
    outf_.reset(new OutputFieldsC{grid(), outf_params});

    mpi_printf(comm, "**** Setting up fields...\n");
    setup_initial_fields(*mflds_);

    init();
  }

  // ----------------------------------------------------------------------
  // setup_initial_fields
  
  void setup_initial_fields(MfieldsState& mflds)
  {
    SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
	switch (m) {
	case HY: return BB_;
	default: return 0.;
	}
      });
  }

private:
  double BB_;
  double Zi_;

  double background_n_;
  double background_Te_;
  double background_Ti_;
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);
  
  auto psc = new PscTestIo;

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
