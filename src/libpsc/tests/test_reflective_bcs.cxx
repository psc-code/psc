#include "gtest/gtest.h"

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "../psc_config.hxx"

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

// ----------------------------------------------------------------------

using BgkMfields = PscConfig::Mfields;
using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

// ======================================================================
// Global parameters

namespace
{
// General PSC parameters
PscParams psc_params;

int nx = 1;
int ny = 4;
int nz = 2;
int nt = 10;

int n_patches_x = 1;
int n_patches_y = 1;
int n_patches_z = 1;

double dx = 10;
double dy = 10;
double dz = 10;

double len_x = nx * dx;
double len_y = ny * dy;
double len_z = nz * dz;

int n_writes = nt;
int out_interval = nt / n_writes;

} // namespace

// ======================================================================
// setupParameters

void setupParameters(int argc, char** argv)
{
  psc_params.nmax = nt;
  psc_params.stats_every = out_interval;
  psc_params.cfl = .75;

  psc_params.write_checkpoint_every_step = 0;
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
  // FIXME add a check to catch mismatch between Dim and n grid points early
  auto domain =
    Grid_t::Domain{{nx, ny, nz},          // n grid points
                   {len_x, len_y, len_z}, // physical lengths
                   {0, 0, 0},             // location of lower corner
                   {n_patches_x, n_patches_y, n_patches_z}}; // n patches

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  // --- generic setup
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 100;

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

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts)
{
  auto injector = mprts.injector();
  injector[0]({{0.0, 36.0, 5.0}, {0.0, 3.0, 0.0}, 1.0, 0});
}

// ======================================================================
// run

static void run(int argc, char** argv)
{
  mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

  // ----------------------------------------------------------------------
  // setup various parameters first

  setupParameters(argc, argv);

  // ----------------------------------------------------------------------
  // Set up grid, state fields, particles

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};
  BgkMfields phi{grid, 1, mflds.ibn()};
  BgkMfields gradPhi{grid, 3, mflds.ibn()};
  BgkMfields divGradPhi{grid, 1, mflds.ibn()};

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  // -- Balance
  psc_params.balance_interval = 0;
  Balance balance{.1};

  // -- Sort
  psc_params.sort_interval = 100;

  // -- Collision
  int collision_interval = 0;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.gauss.check_interval = 1;
  checks_params.continuity.check_interval = 1;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = -1;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output

  // -- output fields
  OutputFieldsParams outf_params{};
  outf_params.fields.pfield.out_interval = out_interval;
  outf_params.moments.pfield.out_interval = out_interval;
  OutputFields<MfieldsState, Mparticles, Dim> outf{grid, outf_params};

  // -- output particles
  OutputParticlesParams outp_params{};
  outp_params.every_step = out_interval;
  outp_params.data_dir = ".";
  outp_params.basename = "prt";
  OutputParticles outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // set up initial conditions

  initializeParticles(balance, grid_ptr, mprts);

  // ----------------------------------------------------------------------
  // run the simulation

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

  psc.integrate();
}

// ======================================================================
// main

// ======================================================================
// Test1

TEST(ReflectiveBcsTest, Test1)
{
  // psc_init(argc, argv);
  // FIXME restore whatever previous functionality there was with options
  int argc = 1;
  char** argv = nullptr;
  psc_init(argc, argv);

  run(argc, argv);

  psc_finalize();
}

int main(int argc, char** argv)
{
  //   MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  //   MPI_Finalize();
  return rc;
}
