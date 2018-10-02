
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
    // initial output / stats
    mpi_printf(grid().comm(), "Performing initial diagnostics.\n");
    diagnostics();

    mpi_printf(grid().comm(), "Initialization complete.\n");
  }

  // ----------------------------------------------------------------------
  // inject_particles

  virtual void inject_particles()
  {}
  
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

private:

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

public:
  const Grid_t& grid() { return *grid_; }

protected:
  double time_start_;

  PscParams p_;
  Grid_t*& grid_;

  std::unique_ptr<MfieldsState> mflds_;
  std::unique_ptr<Mparticles_t> mprts_;

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

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
