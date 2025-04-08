#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"

#include "psc_bgk_util/params_parser.hxx"

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_yz;
#ifdef USE_CUDA
using PscConfig = PscConfig1vbecCuda<Dim>;
#else
using PscConfig = PscConfig1vbecDouble<Dim>;
#endif

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

double electron_temperature;
double ion_temperature;
double electron_mass;
double ion_mass;

double v_upstream_x;
double v_upstream_y;
double v_upstream_z;

double b_x;
double b_y;
double b_z;

double e_x;
double e_y;
double e_z;

int nx;
int ny;
int nz;

int n_patches_x;
int n_patches_y;
int n_patches_z;

double len_x;
double len_y;
double len_z;

int out_interval;

} // namespace

// ======================================================================
// setupParameters

void setupParameters(int argc, char** argv)
{
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " path/to/params\nExiting."
              << std::endl;
    exit(1);
  }
  std::string path_to_params(argv[1]);
  ParsedParams parsedParams(path_to_params);

  psc_params.stats_every = 1000;
  psc_params.cfl = parsedParams.getOrDefault<double>("cfl", .75);
  psc_params.write_checkpoint_every_step = 0;

  electron_temperature = parsedParams.get<double>("electron_temperature");
  ion_temperature = parsedParams.get<double>("ion_temperature");
  electron_mass = parsedParams.get<double>("electron_mass");
  ion_mass = parsedParams.get<double>("ion_mass");

  v_upstream_x = parsedParams.get<double>("v_upstream_x");
  v_upstream_y = parsedParams.get<double>("v_upstream_y");
  v_upstream_z = parsedParams.get<double>("v_upstream_z");

  double b_angle_y_to_x_rad = parsedParams.get<double>("b_angle_y_to_x_rad");
  double b_mag = parsedParams.get<double>("b_mag");
  b_x = b_mag * sin(b_angle_y_to_x_rad);
  b_y = b_mag * cos(b_angle_y_to_x_rad);
  b_z = 0.0;

  e_x = -(v_upstream_y * b_z - v_upstream_z * b_y);
  e_y = -(v_upstream_z * b_x - v_upstream_x * b_z);
  e_z = -(v_upstream_x * b_y - v_upstream_y * b_x);

  nx = parsedParams.get<int>("nx");
  ny = parsedParams.get<int>("ny");
  nz = parsedParams.get<int>("nz");
  psc_params.nmax = parsedParams.get<int>("nt");

  n_patches_x = parsedParams.get<int>("n_patches_x");
  n_patches_y = parsedParams.get<int>("n_patches_y");
  n_patches_z = parsedParams.get<int>("n_patches_z");

  double dx = parsedParams.get<double>("dx");
  double dy = parsedParams.get<double>("dy");
  double dz = parsedParams.get<double>("dz");

  len_x = nx * dx;
  len_y = ny * dy;
  len_z = nz * dz;

  int n_writes = parsedParams.getOrDefault<int>("n_writes", 100);
  out_interval = psc_params.nmax / n_writes;

  std::ifstream src(path_to_params, std::ios::binary);
  std::ofstream dst("params_record.txt", std::ios::binary);
  dst << src.rdbuf();
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
  kinds[KIND_ELECTRON] = {-1.0, electron_mass, "e"};
  kinds[KIND_ION] = {1.0, ion_mass, "i"};

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
  SetupParticles<Mparticles> setup_particles(*grid_ptr);
  setup_particles.centerer = Centering::Centerer(Centering::CC);
  setup_particles.initial_momentum_gamma_correction = true;

  auto init_np = [&](int kind, Double3 crd, int p, Int3 idx,
                     psc_particle_np& np) {
    double temperature =
      np.kind == KIND_ION ? ion_temperature : electron_temperature;
    np.n = 1.0;
    np.p = setup_particles.createMaxwellian(
      {np.kind,
       np.n,
       {v_upstream_x, v_upstream_y, v_upstream_z},
       {temperature, temperature, temperature},
       np.tag});
  };

  partitionAndSetupParticles(setup_particles, balance, grid_ptr, mprts,
                             init_np);
}

// ======================================================================
// fillGhosts

template <typename MF>
void fillGhosts(MF& mfld, int compBegin, int compEnd)
{
  Bnd_ bnd{};
  bnd.fill_ghosts(mfld, compBegin, compEnd);
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  setupFields(mflds, [&](int component, double coords[3]) {
    switch (component) {
      case EX: return e_x;
      case EY: return e_y;
      case EZ: return e_z;
      case HX: return b_x;
      case HY: return b_y;
      case HZ: return b_z;
      default: return 0.0;
    }
  });
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
  checks_params.gauss_every_step = out_interval;
  // checks_params.gauss_dump_always = true;
  checks_params.gauss_threshold = 1e-5;

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
  initializeFields(mflds);

  // ----------------------------------------------------------------------
  // run the simulation

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

  psc.integrate();
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  // psc_init(argc, argv);
  // FIXME restore whatever previous functionality there was with options
  int temp = 1;
  psc_init(temp, argv);

  run(argc, argv);

  psc_finalize();
  return 0;
}
