#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "OutputFieldsDefault.h"
#include "psc_config.hxx"
#include "input_params.hxx"

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

double kmin;
double kmax;
int nk;

double turb_gamma;
double turb_energy_density;
double turb_correlation_length;

int out_interval;

bool mirror_domain;

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
  InputParams parsedParams(path_to_params);

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

  kmin = parsedParams.get<double>("kmin");
  kmax = parsedParams.get<double>("kmax");
  nk = parsedParams.get<int>("nk");

  turb_gamma = parsedParams.get<double>("turb_gamma");
  turb_energy_density = parsedParams.get<double>("turb_energy_density");
  turb_correlation_length = parsedParams.get<double>("turb_correlation_length");

  int n_writes = parsedParams.getOrDefault<int>("n_writes", 100);
  out_interval = psc_params.nmax / n_writes;

  mirror_domain = parsedParams.getOrDefault<bool>("mirror_domain", false);

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
  int y_mult = mirror_domain ? 2 : 1;
  auto domain = Grid_t::Domain{
    {nx, y_mult * ny, nz},          // n grid points
    {len_x, y_mult * len_y, len_z}, // physical lengths
    {0, 0, 0},                      // location of lower corner
    {n_patches_x, y_mult * n_patches_y, n_patches_z}}; // n patches

  auto BND_FLD_Y = mirror_domain ? BND_FLD_PERIODIC : BND_FLD_CONDUCTING_WALL;
  auto BND_PRT_Y = mirror_domain ? BND_PRT_PERIODIC : BND_PRT_REFLECTING;
  auto bc = psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_Y, BND_FLD_PERIODIC},
                          {BND_FLD_PERIODIC, BND_FLD_Y, BND_FLD_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_Y, BND_PRT_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_Y, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, electron_mass, "e"};
  kinds[KIND_ION] = {1.0, ion_mass, "i"};

  // --- generic setup
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 100;

  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// ======================================================================
// initializeParticles

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts)
{
  SetupParticles<Mparticles> setup_particles(*grid_ptr);
  setup_particles.centerer = centering::Centerer(centering::CC);
  setup_particles.initial_momentum_gamma_correction = true;

  auto init_np = [&](int kind, Double3 crd, int p, Int3 idx,
                     psc_particle_np& np) {
    double temperature =
      np.kind == KIND_ION ? ion_temperature : electron_temperature;
    np.n = 1.0;
    double vy_coef = mirror_domain && crd[1] > len_y ? -1 : 1;
    np.p = setup_particles.createMaxwellian(
      {np.kind,
       np.n,
       {v_upstream_x, vy_coef * v_upstream_y, v_upstream_z},
       {temperature, temperature, temperature},
       np.tag});
  };

  partitionAndSetupParticles(setup_particles, balance, grid_ptr, mprts,
                             init_np);
}

// ======================================================================
// initializeFields

double round_to_periodic_k(double len, double k)
{
  // want: k = 2*pi*n/len for an integer n
  double n = round(k * len / (2 * M_PI));
  return 2 * M_PI * n / len;
}

void initializeFields(MfieldsState& mflds)
{
  // turbulence: superposition of Alfven waves, each with a random direction,
  // phase, and polarity. Only one wave of each k-shell is chosen, and that wave
  // is given the entire shell's power.
  auto ks = std::vector<double>(nk);
  auto kxs = std::vector<double>(nk);
  auto kys = std::vector<double>(nk);
  auto kzs = std::vector<double>(nk);

  double dk_over_k = pow(kmax / kmin, 1.0 / double(nk - 1)) - 1.0;

  auto powers = std::vector<double>(nk);

  auto cos_thetas = std::vector<double>(nk);
  auto sin_thetas = std::vector<double>(nk);

  auto phis = std::vector<double>(nk);
  auto phases = std::vector<double>(nk);
  auto alfven_polarities = std::vector<double>(nk);

  auto rng = rng::Uniform<double>();

  double total_power = 0.0;

  for (int i = 0; i < nk; i++) {
    cos_thetas[i] = rng.get(-1.0, 1.0);
    sin_thetas[i] = sin(acos(cos_thetas[i]));

    phis[i] = rng.get(0.0, 2.0 * M_PI);
    phases[i] = rng.get(0.0, 2.0 * M_PI);
    alfven_polarities[i] = rng.get(0.0, 2.0 * M_PI);

    ks[i] = kmin * pow(kmax / kmin, i / double(nk - 1));
    kxs[i] = ks[i] * cos(phis[i]) * sin_thetas[i];
    kys[i] = ks[i] * sin(phis[i]) * sin_thetas[i];
    kzs[i] = ks[i] * cos_thetas[i];

    // force periodicity by "snapping" to it
    kxs[i] = round_to_periodic_k(len_x, kxs[i]);
    kys[i] = round_to_periodic_k(len_y, kys[i]);
    kzs[i] = round_to_periodic_k(len_z, kzs[i]);
    ks[i] = sqrt(sqr(kxs[i]) + sqr(kys[i]) + sqr(kzs[i]));
    LOG_INFO("k[%d] = %f\n", i, ks[i]);

    powers[i] = 1.0 / (1.0 + pow(ks[i] * turb_correlation_length, turb_gamma));
    total_power += powers[i] * pow(ks[i], 3.0) * dk_over_k;
  }

  double db1 = sqrt(turb_energy_density / total_power);
  auto dbs = std::vector<double>(nk);

  auto aaxs = std::vector<double>(nk);
  auto aays = std::vector<double>(nk);
  auto aazs = std::vector<double>(nk);

  auto bbxs = std::vector<double>(nk);
  auto bbys = std::vector<double>(nk);

  for (int i = 0; i < nk; i++) {
    dbs[i] = sqrt(powers[i] * pow(ks[i], 3.0) * dk_over_k);

    aaxs[i] = dbs[i] * cos(alfven_polarities[i]) * cos_thetas[i] * cos(phis[i]);
    aays[i] = dbs[i] * cos(alfven_polarities[i]) * cos_thetas[i] * sin(phis[i]);
    aazs[i] = -dbs[i] * cos(alfven_polarities[i]) * sin_thetas[i];

    bbxs[i] = dbs[i] * sin(alfven_polarities[i]) * sin(phis[i]);
    bbys[i] = -dbs[i] * sin(alfven_polarities[i]) * cos(phis[i]);
  }

  setupFields(mflds, [&](int component, double coords[3]) {
    double e_coef = mirror_domain && coords[1] > len_y ? -1 : 1;
    switch (component) {
      case EX: return e_coef * e_x;
      case EY: return e_coef * e_y;
      case EZ: return e_coef * e_z;

      case HX: {
        double dbx = 0.0;
        for (int i = 0; i < nk; i++) {
          double arg = kxs[i] * coords[0] + kys[i] * coords[1] +
                       kzs[i] * coords[2] + phases[i];

          dbx += aaxs[i] * cos(arg) + bbxs[i] * sin(arg);
        }

        return b_x + dbx;
      }
      case HY: {
        double dby = 0.0;
        for (int i = 0; i < nk; i++) {
          double arg = kxs[i] * coords[0] + kys[i] * coords[1] +
                       kzs[i] * coords[2] + phases[i];

          dby += aays[i] * cos(arg) + bbys[i] * sin(arg);
        }

        return b_y + dby;
      }
      case HZ: {
        double dbz = 0.0;
        for (int i = 0; i < nk; i++) {
          double arg = kxs[i] * coords[0] + kys[i] * coords[1] +
                       kzs[i] * coords[2] + phases[i];

          dbz += aazs[i] * cos(arg);
        }

        return b_z + dbz;
      }
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
  checks_params.gauss.check_interval = out_interval;
  checks_params.continuity.check_interval = out_interval;

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
  DiagEnergies<Mparticles, MfieldsState> oute{grid.comm(), oute_interval};

  // ----------------------------------------------------------------------
  // set up initial conditions

  initializeParticles(balance, grid_ptr, mprts);
  initializeFields(mflds);

  // ----------------------------------------------------------------------
  // run the simulation

  auto psc = makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts,
                                          balance, collision, checks, marder);

  psc.add_diagnostic(&outf);
  psc.add_diagnostic(&outp);
  psc.add_diagnostic(&oute);

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
