
#include <psc_config.h>

#include <psc.h>
#include <psc.hxx>

// small 3d box (heating)
#define TEST_1_HEATING_3D 1
#define TEST_2_FLATFOIL_3D 2
#define TEST_3_NILSON_3D 3
#define TEST_4_SHOCK_3D 4

// EDIT to change test we're running (if TEST is not defined, default is regular 2d flatfoil)
//#define TEST TEST_1_HEATING_3D
//#define TEST TEST_2_FLATFOIL_3D
//#define TEST TEST_3_NILSON_3D
#define TEST TEST_4_SHOCK_3D

#ifdef USE_VPIC
#include "../libpsc/vpic/vpic_iface.h"
#endif

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
#include <marder.hxx>
#include <inject.hxx>
#include <heating.hxx>
#include <setup_particles.hxx>
#include <setup_fields.hxx>

#include "../libpsc/psc_inject/psc_inject_impl.hxx"
#include "../libpsc/psc_heating/psc_heating_impl.hxx"
#include "../libpsc/psc_checks/checks_impl.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#include "../libpsc/cuda/setup_particles_cuda.hxx"
#endif

#include "psc_config.hxx"

#include "heating_spot_foil.hxx"

enum {
  MY_ION,
  MY_ELECTRON,
  N_MY_KINDS,
};

// ======================================================================
// InjectFoil

struct InjectFoilParams
{
  double yl, yh;
  double zl, zh;
  double n;
  double Te, Ti;
};

struct InjectFoil : InjectFoilParams
{
  InjectFoil() = default;
  
  InjectFoil(const InjectFoilParams& params)
    : InjectFoilParams(params)
  {}

  bool is_inside(double crd[3])
  {
    return (crd[1] >= yl && crd[1] <= yh &&
	    crd[2] >= zl && crd[2] <= zh);
  }

  void init_npt(int pop, double crd[3], struct psc_particle_npt *npt)
  {
    if (!is_inside(crd)) {
      npt->n = 0;
      return;
    }
    
    switch (pop) {
    case MY_ION:
      npt->n    = n;
      npt->T[0] = Ti;
      npt->T[1] = Ti;
      npt->T[2] = Ti;
      break;
    case MY_ELECTRON:
      npt->n    = n;
      npt->T[0] = Te;
      npt->T[1] = Te;
      npt->T[2] = Te;
      break;
    default:
      assert(0);
    }
  }
};

// ======================================================================
// PscFlatfoilParams

struct PscFlatfoilParams
{
};

// EDIT to change order / floating point type / cuda / 2d/3d
#if TEST == TEST_1_HEATING_3D || TEST == TEST_2_FLATFOIL_3D || TEST == TEST_3_NILSON_3D || TEST == TEST_4_SHOCK_3D
using dim_t = dim_xyz;
#else
using dim_t = dim_yz;
#endif
using PscConfig = PscConfig1vbecSingle<dim_t>;

// ======================================================================
// PscFlatfoil
//
// eventually, a Psc replacement / derived class, but for now just
// pretending to be something like that

struct PscFlatfoil : Psc<PscConfig>, PscFlatfoilParams
{
  using DIM = PscConfig::dim_t;
  using Heating_t = typename HeatingSelector<Mparticles_t>::Heating;
  using Inject_t = typename InjectSelector<Mparticles_t, MfieldsState, InjectFoil, DIM>::Inject;

  // ----------------------------------------------------------------------
  // ctor
  
  PscFlatfoil()
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    p_.nmax = 5001;
    p_.cfl = 0.75;

    PscFlatfoilParams params;

    BB = 0.;
    Zi = 1.;
    
    // --- for background plasma
    background_n  = .002;
    background_Te = .001;
    background_Ti = .001;
    
    // -- setup particle kinds
    // last population ("e") is neutralizing
    // FIXME, hardcoded mass ratio 100
#if TEST == TEST_4_SHOCK_3D
    Grid_t::Kinds kinds = {{Zi, 100.*Zi, "i"}, { -1., 1., "e"}};
#else
    Grid_t::Kinds kinds = {{Zi, 100.*Zi, "i"}, { -1., 1., "e"}};
#endif
    
    double d_i = sqrt(kinds[MY_ION].m / kinds[MY_ION].q);
    
    mpi_printf(comm, "d_e = %g, d_i = %g\n", 1., d_i);
    mpi_printf(comm, "lambda_De (background) = %g\n", sqrt(background_Te));
    
#if TEST == TEST_4_SHOCK_3D
    BB = 0.02;
    background_n = .01;
    background_Te = .002;
    background_Ti = .002;
    p_.nmax = 100002;
#endif

#if TEST == TEST_3_NILSON_3D
    background_n = .02;
    p_.nmax = 101;
#endif
    
#if TEST == TEST_1_HEATING_3D
    background_n  = 1.0;
    
    params.checks_params.continuity_every_step = 1;
    params.checks_params.continuity_threshold = 1e-12;
    params.checks_params.continuity_verbose = true;
    
    params.checks_params.gauss_every_step = 1;
    // eventually, errors accumulate above 1e-10, but it should take a long time
    params.checks_params.gauss_threshold = 1e-10;
    params.checks_params.gauss_verbose = true;

    p_.collision_interval = 0;
#endif

    static_cast<PscFlatfoilParams&>(*this) = params;
  
    // --- setup domain
    Grid_t::Real3 LL = { 1., 400.*4, 400. }; // domain size (in d_e)
    Int3 gdims = { 1, 4096, 1024 }; // global number of grid points
    Int3 np = { 1, 64, 16 }; // division into patches
    
#if TEST == TEST_4_SHOCK_3D
    LL = { 256., 256., 512. }; // domain size (in d_e)
    gdims = { 256, 256, 512 }; // global number of grid points
    np = { 16, 16, 32 }; // division into patches
#endif
    
#if TEST == TEST_3_NILSON_3D
    LL = { 400., 200., 800. }; // domain size (in d_e)
    gdims = { 64, 32, 128 }; // global number of grid points
    np = { 2, 1, 4 }; // division into patches
#endif
  
#if TEST == TEST_2_FLATFOIL_3D
    LL = { 400., 400.*4, 400. }; // domain size (in d_e)
    gdims = { 16, 64, 16 }; // global number of grid points
    np = { 1, 4, 1 }; // division into patches
#endif
  
#if TEST == TEST_1_HEATING_3D
    LL = { 2., 2., 2. }; // domain size (in d_e)
    gdims = { 8, 8, 8 }; // global number of grid points
    np = { 1, 1, 1 }; // division into patches
#endif
  
    auto grid_domain = Grid_t::Domain{gdims, LL, -.5 * LL, np};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 100;
#if TEST == TEST_4_SHOCK_3D
    norm_params.nicell = 100;
#endif
#if TEST == TEST_3_NILSON_3D
    norm_params.nicell = 50;
#endif

    double dt = p_.cfl * courant_length(grid_domain);
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles_t{grid()});

    // -- Balance
    balance_interval = 0;
    double balance_factor_fields = .1;
    bool balance_print_loads = true;
    bool balance_write_loads = false;
    balance_.reset(new Balance_t{balance_interval, balance_factor_fields,
	  balance_print_loads, balance_write_loads});

    // -- Sort
    // FIXME, needs a way to make it gets set?
    sort_interval = 10;

    // -- Collision
    int collision_interval = 10;
    double collision_nu = .1;
#if TEST == TEST_4_SHOCK_3D
    collision_interval = 0;
#endif
#if TEST == TEST_3_NILSON_3D
    collision_interval = 0;
#endif
#if TEST == TEST_2_FLATFOIL_3D
    collision_interval = 0;
#endif
#if TEST == TEST_1_HEATING_3D
    collision_interval = 0;
#endif
    collision_.reset(new Collision_t{grid(), collision_interval, collision_nu});

    // -- Checks
    ChecksParams checks_params;

    checks_params.continuity_every_step = -1;
    checks_params.continuity_threshold = 1e-6;
    checks_params.continuity_verbose = true;
    checks_params.continuity_dump_always = false;
    
    checks_params.gauss_every_step = -1;
    checks_params.gauss_threshold = 1e-6;
    checks_params.gauss_verbose = true;
    checks_params.gauss_dump_always = false;

    checks_.reset(new Checks_t{grid(), comm, checks_params});

    // -- Marder correction
    double marder_diffusion = 0.9;
    int marder_loop = 3;
    bool marder_dump = false;
    marder_interval = 0*5;
    marder_.reset(new Marder_t(grid(), marder_diffusion, marder_loop, marder_dump));

    // -- Heating
    auto heating_foil_params = HeatingSpotFoilParams{};
    heating_foil_params.zl = -1. * d_i;
    heating_foil_params.zh =  1. * d_i;
    heating_foil_params.xc =  0. * d_i;
    heating_foil_params.yc =  0. * d_i;
#if TEST == TEST_4_SHOCK_3D
    heating_foil_params.rH = 100000 * d_i;
    heating_foil_params.T  = .06;
#else
    heating_foil_params.rH = 3. * d_i;
    heating_foil_params.T  = .04;
#endif
    heating_foil_params.Mi = kinds[MY_ION].m;
    auto heating_spot = HeatingSpotFoil{heating_foil_params};

    heating_interval = 20;
#if TEST == TEST_1_HEATING_3D || TEST == TEST_2_FLATFOIL_3D
    heating_interval = 0;
#endif
    heating_begin = 0;
    heating_end = 10000000;
    heating_.reset(new Heating_t{grid(), heating_interval, MY_ELECTRON, heating_spot});

    // -- Particle injection
    // -- setup injection
    auto inject_foil_params = InjectFoilParams{};
    inject_foil_params.yl = -100000. * d_i;
    inject_foil_params.yh =  100000. * d_i;
    double target_zwidth =  1.;
    inject_foil_params.zl = - target_zwidth * d_i;
    inject_foil_params.zh =   target_zwidth * d_i;
#if TEST == TEST_4_SHOCK_3D
    inject_foil_params.n  = 2.5;
    inject_foil_params.Te = .002;
    inject_foil_params.Ti = .002;
#else
    inject_foil_params.n  = 1.;
    inject_foil_params.Te = .001;
    inject_foil_params.Ti = .001;
#endif
    inject_target = InjectFoil{inject_foil_params};
    inject_interval = 20;
#if TEST == TEST_4_SHOCK_3D || TEST == TEST_3_NILSON_3D || TEST == TEST_1_HEATING_3D
    inject_interval = 0;
#endif
#if TEST == TEST_2_FLATFOIL_3D
    inject_interval = 5;
#endif
    
    
    int inject_tau = 40;
    inject_.reset(new Inject_t{grid(), inject_interval, inject_tau, MY_ELECTRON, inject_target});
    
    // -- output fields
    OutputFieldsCParams outf_params;
    outf_.reset(new OutputFieldsC{grid(), outf_params});

    // --- partition particles and initial balancing
    mpi_printf(comm, "**** Partitioning...\n");
    auto n_prts_by_patch_old = setup_initial_partition();
    auto n_prts_by_patch_new = balance_->initial(n_prts_by_patch_old);
    // balance::initial does not rebalance particles, because the old way of doing this
    // does't even have the particle data structure created yet -- FIXME?
    mprts_->reset(grid());
    
    mpi_printf(comm, "**** Setting up particles...\n");
    setup_initial_particles(*mprts_, n_prts_by_patch_new);
    
    mpi_printf(comm, "**** Setting up fields...\n");
    setup_initial_fields(*mflds_);

    init();
  }

  void init_npt(int kind, double crd[3], psc_particle_npt& npt)
  {
    switch (kind) {
    case MY_ION:
      npt.n    = background_n;
      npt.T[0] = background_Ti;
      npt.T[1] = background_Ti;
      npt.T[2] = background_Ti;
      break;
    case MY_ELECTRON:
      npt.n    = background_n;
      npt.T[0] = background_Te;
      npt.T[1] = background_Te;
      npt.T[2] = background_Te;
      break;
    default:
      assert(0);
    }
      
    if (inject_target.is_inside(crd)) {
      // replace values above by target values
      inject_target.init_npt(kind, crd, &npt);
    }
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_partition
  
  std::vector<uint> setup_initial_partition()
  {
    SetupParticles<Mparticles_t> setup_particles;
    setup_particles.fractional_n_particles_per_cell = true; // FIXME, should use same setup_particles for partition/setup
    setup_particles.neutralizing_population = MY_ELECTRON;
    return setup_particles.setup_partition(grid(), [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_particles
  
  void setup_initial_particles(Mparticles_t& mprts, std::vector<uint>& n_prts_by_patch)
  {
    SetupParticles<Mparticles_t> setup_particles; // FIXME, injection uses another setup_particles, which won't have those settings
    setup_particles.fractional_n_particles_per_cell = true;
    setup_particles.neutralizing_population = MY_ELECTRON;
    setup_particles.setup_particles(mprts, n_prts_by_patch, [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }

  // ----------------------------------------------------------------------
  // setup_initial_fields
  
  void setup_initial_fields(MfieldsState& mflds)
  {
    SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
	switch (m) {
	case HY: return BB;
	default: return 0.;
	}
      });
  }

  void inject_particles() override
  {
    static int pr_inject, pr_heating;
    if (!pr_inject) {
      pr_inject = prof_register("inject", 1., 0, 0);
      pr_heating = prof_register("heating", 1., 0, 0);
    }

    auto comm = grid().comm();
    auto timestep = grid().timestep();
    
    if (inject_interval > 0 && timestep % inject_interval == 0) {
      mpi_printf(comm, "***** Performing injection...\n");
      prof_start(pr_inject);
      (*inject_)(*mprts_);
      prof_stop(pr_inject);
    }
      
    // only heating between heating_tb and heating_te
    if (timestep >= heating_begin && timestep < heating_end &&
	heating_interval > 0 && timestep % heating_interval == 0) {
      mpi_printf(comm, "***** Performing heating...\n");
      prof_start(pr_heating);
      (*heating_)(*mprts_);
      prof_stop(pr_heating);
    }
  }

protected:
  std::unique_ptr<Heating_t> heating_;
  std::unique_ptr<Inject_t> inject_;

public: // these don't need to be public, but oh well...
  double BB;
  double Zi;

  double background_n;
  double background_Te;
  double background_Ti;

  int inject_interval;
  InjectFoil inject_target;

  int heating_begin;
  int heating_end;
  int heating_interval;
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
#ifdef USE_VPIC
  vpic_base_init(&argc, &argv);
#else
  MPI_Init(&argc, &argv);
#endif
  libmrc_params_init(argc, argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);

  auto psc = new PscFlatfoil;

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
