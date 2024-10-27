
#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"

// ======================================================================
// PscFlatfoilParams

struct PscFlatfoilParams
{
  double theta = 0;
  double k = 2 * 2. * M_PI / 10.;
  double amplitude_s = 1.;
  double amplitude_p = 0.;
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
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_xyz;

using PscConfig = PscConfig1vbecSingle<Dim>;

using Writer = WriterDefault; // can choose WriterMrc, WriterAdios2

// ======================================================================
// Moment_n_Selector
//
// FIXME, should go away eventually

template <typename Mparticles, typename Dim, typename Enable = void>
struct Moment_n_Selector
{
  using type = Moment_n_1st<MfieldsSingle::Storage, Dim>;
};

// ----------------------------------------------------------------------

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;
using Moment_n = typename Moment_n_Selector<Mparticles, Dim>::type;

// ======================================================================
// setupParameters

void setupParameters()
{
  // -- set some generic PSC parameters
  psc_params.nmax = 10000001; // 5001;
  psc_params.cfl = 0.75;
  psc_params.write_checkpoint_every_step = 1000;
  psc_params.stats_every = 1;

  // -- start from checkpoint:
  //
  // Uncomment when wanting to start from a checkpoint, ie.,
  // instead of setting up grid, particles and state fields here,
  // they'll be read from a file
  // FIXME: This parameter would be a good candidate to be provided
  // on the command line, rather than requiring recompilation when change.

  // read_checkpoint_filename = "checkpoint_500.bp";

  // -- Set some parameters specific to this case
  g.theta = 0;
  g.k = 2 * 2. * M_PI / 10.;
  g.amplitude_s = 1.;
  g.amplitude_p = 0.;
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
  Grid_t::Real3 LL = {10, 10., 10.}; // domain size (in d_e)
  Int3 gdims = {50, 50, 50};         // global number of grid points
  Int3 np = {5, 5, 5};               // division into patches

  Grid_t::Domain domain{gdims, LL, -.5 * LL, np};

  psc::grid::BC bc{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                   {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                   {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                   {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  // -- setup particle kinds
  Grid_t::Kinds kinds;

  // -- setup normalization
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 100;

  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  mpi_printf(MPI_COMM_WORLD, "dt = %g\n", dt);

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
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  double theta_rad = g.theta * 2. * M_PI / 360.;
  double ky = g.k * sin(theta_rad), kz = g.k * cos(theta_rad);
  double hat_y = cos(theta_rad);
  double hat_z = -sin(theta_rad);

  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case EX: return g.amplitude_s * sin(ky * crd[1] + kz * crd[2]);
      case HY: return g.amplitude_s * hat_y * sin(ky * crd[1] + kz * crd[2]);
      case HZ: return g.amplitude_s * hat_z * sin(ky * crd[1] + kz * crd[2]);

      case HX: return g.amplitude_p * sin(ky * crd[1] + kz * crd[2]);
      case EY: return g.amplitude_p * -hat_y * sin(ky * crd[1] + kz * crd[2]);
      case EZ: return g.amplitude_p * -hat_z * sin(ky * crd[1] + kz * crd[2]);
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

  Mparticles mprts(grid);
  MfieldsState mflds(grid);
  if (!read_checkpoint_filename.empty()) {
    read_checkpoint(read_checkpoint_filename, grid, mprts, mflds);
  }

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  // -- Balance
  psc_params.balance_interval = 500;
  Balance balance{3};

  // -- Sort
  psc_params.sort_interval = 10;

  // -- Collision
  int collision_interval = -1;
  double collision_nu = 1.;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = -1;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output
  //
  // FIXME, this really is too complicated and not very flexible

  // -- output fields
  OutputFieldsItemParams outf_item_params{};
  OutputFieldsParams outf_params{};
  outf_item_params.pfield.out_interval = 1;

  outf_params.fields = outf_item_params;
  OutputFields<MfieldsState, Mparticles, Dim, Writer> outf{grid, outf_params};

  // -- output particles
  OutputParticlesParams outp_params{};
  OutputParticles outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // setup initial conditions

  if (read_checkpoint_filename.empty()) {
    initializeFields(mflds);
  }

  // ----------------------------------------------------------------------
  // hand off to PscIntegrator to run the simulation

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

  MEM_STATS();
  psc.integrate();
  MEM_STATS();
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  psc_init(argc, argv);

  run();

  MEM_STATS();

  psc_finalize();
  return 0;
}
