
#include <psc_config.h>

#include <psc.h>
#include <psc.hxx>

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

#include "inject_impl.hxx"
#include "../libpsc/psc_heating/psc_heating_impl.hxx"
#include "../libpsc/psc_checks/checks_impl.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
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

// EDIT to change order / floating point type / cuda / 2d/3d
using dim_t = dim_xyz;
using PscConfig = PscConfig1vbecSingle<dim_t>;
//using PscConfig = PscConfig1vbecCuda<dim_t>;

// ======================================================================
// PscFlatfoil

struct PscFlatfoil : Psc<PscConfig>
{
  using DIM = PscConfig::dim_t;
  using Heating_t = typename HeatingSelector<Mparticles_t>::Heating;
  using Inject_t = typename InjectSelector<Mparticles_t, InjectFoil, DIM>::Inject;

  // ----------------------------------------------------------------------
  // ctor
  
  PscFlatfoil()
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
#if 1
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 40, 80, 4 }; // division into patches
#else
    Int3 gdims = { 40, 10, 20}; // global number of grid points
    Int3 np = { 4, 1, 2 }; // division into patches
#endif

    if (dim::InvarX::value) { ibn[0] = 0; } // FIXME, wrong place, not for VPIC...
    
    auto grid_domain = Grid_t::Domain{gdims, LL, -.5 * LL, np};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 50;

    double dt = p_.cfl * courant_length(grid_domain);
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles_t{grid()});

    // -- Balance
    balance_interval = 50;
    balance_.reset(new Balance_t{balance_interval, .1, false});

    // -- Sort
    // FIXME, needs a way to make sure it gets set?
    sort_interval = 10;

    // -- Collision
    int collision_interval = 10;
    double collision_nu = .1;
    collision_.reset(new Collision_t{grid(), collision_interval, collision_nu});

    // -- Checks
    ChecksParams checks_params{};
    checks_params.continuity_every_step = 50;
    checks_params.continuity_threshold = 1e-5;
    checks_params.continuity_verbose = false;
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
    heating_foil_params.rH = 3. * d_i;
    heating_foil_params.T  = .04;
    heating_foil_params.Mi = kinds[MY_ION].m;
    auto heating_spot = HeatingSpotFoil{heating_foil_params};

    heating_interval_ = 20;
    heating_begin_ = 0;
    heating_end_ = 10000000;
    heating_.reset(new Heating_t{grid(), heating_interval_, MY_ELECTRON, heating_spot});

    // -- Particle injection
    auto inject_foil_params = InjectFoilParams{};
    inject_foil_params.yl = -100000. * d_i;
    inject_foil_params.yh =  100000. * d_i;
    double target_zwidth =  1.;
    inject_foil_params.zl = - target_zwidth * d_i;
    inject_foil_params.zh =   target_zwidth * d_i;
    inject_foil_params.n  = 1.;
    inject_foil_params.Te = .001;
    inject_foil_params.Ti = .001;
    inject_target_ = InjectFoil{inject_foil_params};
    inject_interval_ = 20;
    
    int inject_tau = 40;
    inject_.reset(new Inject_t{grid(), inject_interval_, inject_tau, MY_ELECTRON, inject_target_});
    
    // -- output fields
    OutputFieldsCParams outf_params;
    outf_params.output_fields = "e,h,j,n_1st_single,v_1st_single,T_1st_single";
    outf_params.pfield_step = 200;
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
      npt.n    = background_n_;
      npt.T[0] = background_Ti_;
      npt.T[1] = background_Ti_;
      npt.T[2] = background_Ti_;
      break;
    case MY_ELECTRON:
      npt.n    = background_n_;
      npt.T[0] = background_Te_;
      npt.T[1] = background_Te_;
      npt.T[2] = background_Te_;
      break;
    default:
      assert(0);
    }
      
    if (inject_target_.is_inside(crd)) {
      // replace values above by target values
      inject_target_.init_npt(kind, crd, &npt);
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
	case HY: return BB_;
	default: return 0.;
	}
      });
  }

  // ----------------------------------------------------------------------
  // inject_particles

  void inject_particles() override
  {
    static int pr_inject, pr_heating;
    if (!pr_inject) {
      pr_inject = prof_register("inject", 1., 0, 0);
      pr_heating = prof_register("heating", 1., 0, 0);
    }

    auto comm = grid().comm();
    auto timestep = grid().timestep();
    
    if (inject_interval_ > 0 && timestep % inject_interval_ == 0) {
      mpi_printf(comm, "***** Performing injection...\n");
      prof_start(pr_inject);
      (*inject_)(*mprts_);
      prof_stop(pr_inject);
    }
      
    // only heating between heating_tb and heating_te
    if (timestep >= heating_begin_ && timestep < heating_end_ &&
	heating_interval_ > 0 && timestep % heating_interval_ == 0) {
      mpi_printf(comm, "***** Performing heating...\n");
      prof_start(pr_heating);
      (*heating_)(*mprts_);
      prof_stop(pr_heating);
    }
  }

protected:
  std::unique_ptr<Heating_t> heating_;
  std::unique_ptr<Inject_t> inject_;

private:
  double BB_;
  double Zi_;

  double background_n_;
  double background_Te_;
  double background_Ti_;

  int inject_interval_;
  InjectFoil inject_target_;

  int heating_begin_;
  int heating_end_;
  int heating_interval_;
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);
  
  auto psc = new PscFlatfoil;

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
