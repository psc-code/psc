
#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "psc_config.hxx"

#include "../libpsc/psc_heating/psc_heating_impl.hxx"
#include "heating_spot_foil.hxx"
#include "inject_impl.hxx"

// ======================================================================
// Particle kinds
//
// Particle kinds can be used to define different species, or different
// populations of the same species
//
// Here, we only enumerate the types, the actual information gets set up later.
// The last kind (MY_ELECTRON) will be used as "neutralizing kind", ie, in the
// initial setup, the code will add as many electrons as there are ions in a
// cell, at the same position, to ensure the initial plasma is neutral
// (so that Gauss's Law is satisfied).

enum
{
  MY_ION,
  MY_ELECTRON,
  N_MY_KINDS,
};

// ======================================================================
// InjectFoil
//
// This class describes the actual foil, specifically, where the foil
// is located and what its plasma parameters are

struct InjectFoilParams
{
  double yl, yh;
  double zl, zh;
  double n;
  double Te, Ti;
};

class InjectFoil : InjectFoilParams
{
public:
  InjectFoil() = default;

  InjectFoil(const InjectFoilParams& params) : InjectFoilParams(params) {}

  bool is_inside(double crd[3])
  {
    return (crd[1] >= yl && crd[1] <= yh && crd[2] >= zl && crd[2] <= zh);
  }

  void init_npt(int pop, double crd[3], psc_particle_npt& npt)
  {
    if (!is_inside(crd)) {
      npt.n = 0;
      return;
    }

    switch (pop) {
      case MY_ION:
        npt.n = n;
        npt.T[0] = Ti;
        npt.T[1] = Ti;
        npt.T[2] = Ti;
        break;
      case MY_ELECTRON:
        npt.n = n;
        npt.T[0] = Te;
        npt.T[1] = Te;
        npt.T[2] = Te;
        break;
      default: assert(0);
    }
  }
};

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_yz;

#if 1
#ifdef USE_CUDA
using PscConfig = PscConfig1vbecCuda<Dim>;
#else
using PscConfig = PscConfig1vbecSingle<Dim>;
#endif

#else

#include "particle_with_id.h"

using PscConfig =
  PscConfig_<Dim, MparticlesSimple<ParticleWithId<float>>, MfieldsStateSingle,
             MfieldsSingle, PscConfigPushParticles1vbec>;

#endif

// ----------------------------------------------------------------------

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;
using Inject = typename InjectSelector<Mparticles, InjectFoil, Dim>::Inject;
using Heating = typename HeatingSelector<Mparticles>::Heating;

// ======================================================================
// PscFlatfoilParams

struct PscFlatfoilParams
{
  double BB;
  double Zi;
  double mass_ratio;

  double background_n;
  double background_Te;
  double background_Ti;

  int inject_interval;

  int heating_begin;
  int heating_end;
  int heating_interval;

  // The following parameters are calculated from the above / and other
  // information

  double d_i;
};

// ======================================================================
// Global parameters
//
// I'm not a big fan of global parameters, but they're only for
// this particular case and they help make things simpler.

// An "anonymous namespace" makes these variables visible in this source file
// only
namespace
{

// Parameters specific to this case. They don't really need to be collected in a
// struct, but maybe it's nice that they are

PscFlatfoilParams g;

std::string read_checkpoint_filename;

// This is a set of generic PSC params (see include/psc.hxx),
// like number of steps to run, etc, which also should be set by the case
PscParams psc_params;

} // namespace

// ======================================================================
// setupParameters

void setupParameters()
{
  // -- set some generic PSC parameters
  psc_params.nmax = 2001; // 5001;
  psc_params.cfl = 0.75;
  psc_params.write_checkpoint_every_step = 500;

  // -- start from checkpoint:
  //
  // Uncomment when wanting to start from a checkpoint, ie.,
  // instead of setting up grid, particles and state fields here,
  // they'll be read from a file
  // FIXME: This parameter would be a good candidate to be provided
  // on the command line, rather than requiring recompilation when change.

  // read_checkpoint_filename = "checkpoint_500.bp";

  // -- Set some parameters specific to this case
  g.BB = 0.;
  g.Zi = 1.;
  g.mass_ratio = 100.; // 25.

  g.background_n = .002;
  g.background_Te = .001;
  g.background_Ti = .001;
}

// ======================================================================
// setupGrid
//
// This helper function is responsible for setting up the "Grid",
// which is really more than just the domain and its decomposition, it
// also encompasses PC normalization parameters, information about the
// particle kinds, etc.

Grid_t* setupGrid()
{
  // --- setup domain
#if 0
  Grid_t::Real3 LL = {384., 384.*2., 384.*6}; // domain size (in d_e)
  Int3 gdims = {384, 384*2, 384*6}; // global number of grid points
  Int3 np = {12, 24, 72}; // division into patches
#endif
#if 0
  Grid_t::Real3 LL = {192., 192.*2, 192.*6}; // domain size (in d_e)
  Int3 gdims = {192, 192*2, 192*6}; // global number of grid points
  Int3 np = {6, 12, 36}; // division into patches
  // -> patch size 32 x 32 x 32
#endif
#if 0
  Grid_t::Real3 LL = {32., 32.*2., 32.*6 }; // domain size (in d_e)
  Int3 gdims = {32, 32*2, 32*6}; // global number of grid points
  Int3 np = { 1, 2, 6 }; // division into patches
#endif
#if 0
  Grid_t::Real3 LL = {1., 1600., 400.}; // domain size (in d_e)
  // Int3 gdims = {40, 10, 20}; // global number of grid points
  // Int3 np = {4, 1, 2; // division into patches
  Int3 gdims = {1, 2048, 512}; // global number of grid points
  Int3 np = {1, 64, 16}; // division into patches
#endif
#if 1
  Grid_t::Real3 LL = {1., 800., 200.}; // domain size (in d_e)
  Int3 gdims = {1, 1024, 256};         // global number of grid points
  Int3 np = {1, 8, 2};                 // division into patches
#endif

  Grid_t::Domain domain{gdims, LL, -.5 * LL, np};

  psc::grid::BC bc{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                   {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                   {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                   {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  // -- setup particle kinds
  // last population ("e") is neutralizing
  Grid_t::Kinds kinds(N_MY_KINDS);
  kinds[MY_ION] = {g.Zi, g.mass_ratio * g.Zi, "i"};
  kinds[MY_ELECTRON] = {-1., 1., "e"};

  g.d_i = sqrt(kinds[MY_ION].m / kinds[MY_ION].q);

  mpi_printf(MPI_COMM_WORLD, "d_e = %g, d_i = %g\n", 1., g.d_i);
  mpi_printf(MPI_COMM_WORLD, "lambda_De (background) = %g\n",
             sqrt(g.background_Te));

  // -- setup normalization
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 50;

  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  Int3 ibn = {2, 2, 2};
  if (Dim::InvarX::value) {
    ibn[0] = 0;
  }
  if (Dim::InvarY::value) {
    ibn[1] = 0;
  }
  if (Dim::InvarZ::value) {
    ibn[2] = 0;
  }

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// ======================================================================
// initializeParticles

void initializeParticles(SetupParticles<Mparticles>& setup_particles,
                         Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts,
                         InjectFoil& inject_target)
{
  // -- set particle initial condition
  partitionAndSetupParticles(
    setup_particles, balance, grid_ptr, mprts,
    [&](int kind, double crd[3], psc_particle_npt& npt) {
      switch (kind) {
        case MY_ION:
          npt.n = g.background_n;
          npt.T[0] = g.background_Ti;
          npt.T[1] = g.background_Ti;
          npt.T[2] = g.background_Ti;
          break;
        case MY_ELECTRON:
          npt.n = g.background_n;
          npt.T[0] = g.background_Te;
          npt.T[1] = g.background_Te;
          npt.T[2] = g.background_Te;
          break;
        default: assert(0);
      }

      if (inject_target.is_inside(crd)) {
        // replace values above by target values
        inject_target.init_npt(kind, crd, npt);
      }
    });
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case HY: return g.BB;
      default: return 0.;
    }
  });
}

// ======================================================================
// run
//
// This is basically the main function of this run,
// which sets up everything and then uses PscIntegrator to run the
// simulation

void run()
{
  mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

  // ----------------------------------------------------------------------
  // setup various parameters first

  setupParameters();

  // ----------------------------------------------------------------------
  // Set up grid, state fields, particles

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  // -- Balance
  psc_params.balance_interval = 300;
  Balance balance{psc_params.balance_interval, 3, true};

  // -- Sort
  psc_params.sort_interval = 10;

  // -- Collision
  int collision_interval = 10;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.continuity_every_step = 50;
  checks_params.continuity_threshold = 1e-5;
  checks_params.continuity_verbose = false;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = 0 * 5;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output
  //
  // FIXME, this really is too complicated and not very flexible

  // -- output fields
  OutputFieldsCParams outf_params{};
  outf_params.pfield_step = 200;
  auto outf = defaultOutputFieldsC<Mparticles>(grid, outf_params);

  // -- output particles
  OutputParticlesParams outp_params{};
  outp_params.every_step = 0;
  outp_params.data_dir = ".";
  outp_params.basename = "prt";
  OutputParticles outp{grid, outp_params};

  int oute_interval = 100;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // Set up objects specific to the flatfoil case

  // -- Heating
  HeatingSpotFoilParams heating_foil_params{};
  heating_foil_params.zl = -1. * g.d_i;
  heating_foil_params.zh = 1. * g.d_i;
  heating_foil_params.xc = 0. * g.d_i;
  heating_foil_params.yc = 0. * g.d_i;
  heating_foil_params.rH = 3. * g.d_i;
  heating_foil_params.T = .04;
  heating_foil_params.Mi = grid.kinds[MY_ION].m;
  HeatingSpotFoil heating_spot{heating_foil_params};

  g.heating_interval = 20;
  g.heating_begin = 0;
  g.heating_end = 10000000;
  auto& heating =
    *new Heating{grid, g.heating_interval, MY_ELECTRON, heating_spot};

  // -- Particle injection
  InjectFoilParams inject_foil_params;
  inject_foil_params.yl = -100000. * g.d_i;
  inject_foil_params.yh = 100000. * g.d_i;
  double target_zwidth = 1.;
  inject_foil_params.zl = -target_zwidth * g.d_i;
  inject_foil_params.zh = target_zwidth * g.d_i;
  inject_foil_params.n = 1.;
  inject_foil_params.Te = .001;
  inject_foil_params.Ti = .001;
  InjectFoil inject_target{inject_foil_params};

  g.inject_interval = 20;
  int inject_tau = 40;

  SetupParticles<Mparticles> setup_particles(grid);
  setup_particles.fractional_n_particles_per_cell = true;
  setup_particles.neutralizing_population = MY_ELECTRON;

  Inject inject{grid,        g.inject_interval, inject_tau,
                MY_ELECTRON, inject_target,     setup_particles};

  auto lf_inject = [&](const Grid_t& grid, Mparticles& mprts) {
    static int pr_inject, pr_heating;
    if (!pr_inject) {
      pr_inject = prof_register("inject", 1., 0, 0);
      pr_heating = prof_register("heating", 1., 0, 0);
    }

    auto comm = grid.comm();
    auto timestep = grid.timestep();

    if (g.inject_interval > 0 && timestep % g.inject_interval == 0) {
      mpi_printf(comm, "***** Performing injection...\n");
      prof_start(pr_inject);
      inject(mprts);
      prof_stop(pr_inject);
    }

    // only heating between heating_tb and heating_te
    if (timestep >= g.heating_begin && timestep < g.heating_end &&
        g.heating_interval > 0 && timestep % g.heating_interval == 0) {
      mpi_printf(comm, "***** Performing heating...\n");
      prof_start(pr_heating);
      heating(mprts);
      prof_stop(pr_heating);
    }
  };

  // ----------------------------------------------------------------------
  // setup initial conditions

  if (read_checkpoint_filename.empty()) {
    initializeParticles(setup_particles, balance, grid_ptr, mprts,
                        inject_target);
    initializeFields(mflds);
  } else {
    read_checkpoint(read_checkpoint_filename, *grid_ptr, mprts, mflds);
  }
  
  // ----------------------------------------------------------------------
  // hand off to PscIntegrator to run the simulation

  auto psc = makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts,
                                          balance, collision, checks, marder,
                                          diagnostics, lf_inject);

  psc.integrate();
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  psc_init(argc, argv);

  run();

  psc_finalize();
  return 0;
}
