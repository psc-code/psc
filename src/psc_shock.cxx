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

double turb_db2;
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

  turb_db2 = parsedParams.get<double>("turb_dB^2");
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
  // literally the power of each shell, such that
  // 1/2 <dB^2> = int_{k=0}^\infty shellpower(k) dk.
  // Note that regardless of dimensionality, shellpower(k) ~ k^{-5/3},
  // i.e., shellpower(k) includes the Jacobian.
  auto shell_power = [&](double k) {
    return 1.0 / (1.0 + pow(k * turb_correlation_length, 5.0 / 3.0));
  };

  // 0. constant (k=0) part

  setupFields(mflds, [&](int component, double coords[3]) {
    switch (component) {
      case HX: return b_x;
      case HY: return b_y;
      case HZ: return b_z;
      default: return 0.0;
    }
  });

  // 1. compute values of |k|

  double dkx = 2.0 * M_PI / len_x;
  double dky = 2.0 * M_PI / len_y;
  double dkz = 2.0 * M_PI / len_z;

  double dk = pow(dkx * dky * dkz, 1.0 / 3.0);
  double kmax = sqrt(sqr(dkx * nx) + sqr(dky * ny) + sqr(dkz * nz)) / 2.0;
  int nk = kmax / dk + 1;

  auto k_mins = std::vector<double>(nk);
  for (int i = 0; i < nk; i++) {
    k_mins[i] = dk * i;
  }

  // 2. compute number of degeneracies of each k-shell

  int ix_min = -nx / 2;
  int iy_min = -ny / 2;
  int iz_min = -nz / 2;

  // inject in only half of k-space, since +k and -k modes are indistinguishable
  if (nx != 1) {
    ix_min = 0;
  } else if (ny != 1) {
    iy_min = 0;
  } else {
    iz_min = 0;
  }

  int ix_max = nx / 2 + 1;
  int iy_max = ny / 2 + 1;
  int iz_max = nz / 2 + 1;

  auto n_cells_per_shell = std::vector<int>(nk);

  for (int ix = ix_min; ix < ix_max; ix++) {
    for (int iy = iy_min; iy < iy_max; iy++) {
      for (int iz = iz_min; iz < iz_max; iz++) {
        double kx = ix * dkx;
        double ky = iy * dky;
        double kz = iz * dkz;
        double k = sqrt(sqr(kx) + sqr(ky) + sqr(kz));

        if (k == 0) {
          continue;
        }

        int idx = k / dk;
        n_cells_per_shell[idx] += 1;
      }
    }
  }

  // 3. compute powers for each shell

  auto shell_energy_densities = std::vector<double>(nk);
  double initial_total_energy_density = 0.0;

  for (int i = 0; i < nk; i++) {
    double k_shell = k_mins[i] + 0.5 * dky;
    shell_energy_densities[i] = shell_power(k_shell) * dky;
    initial_total_energy_density += shell_energy_densities[i];
  }

  // scale shell powers to match desired energy density

  double target_total_energy_density = 0.5 * turb_db2;
  for (int i = 0; i < nk; i++) {
    shell_energy_densities[i] *=
      target_total_energy_density / initial_total_energy_density;
  }

  // 4. inject each mode at a random phase and polarization

  // TODO randomize the seed based on e.g. time
  int seed = 5; // all processes must use same seed to ensure B is continuous
  auto rng = rng::Uniform<double>(0.0, 1.0, seed);
  const auto& grid = mflds.grid();

  for (int ix = ix_min; ix < ix_max; ix++) {
    for (int iy = iy_min; iy < iy_max; iy++) {
      for (int iz = iz_min; iz < iz_max; iz++) {
        double kx = ix * dkx;
        double ky = iy * dky;
        double kz = iz * dkz;

        double kxy = sqrt(sqr(kx) + sqr(ky));
        double k = sqrt(sqr(kx) + sqr(ky) + sqr(kz));

        if (k == 0.0) {
          // no power for the constant mode, and avoid division-by-0
          continue;
        }

        int idx = k / dk;

        double proportion_of_shell_energy_density_in_cell =
          1.0 / n_cells_per_shell[idx];
        double shell_energy_density = shell_energy_densities[idx];
        double lambda = 4; // experimentally obtained
        double cell_energy_density =
          lambda * proportion_of_shell_energy_density_in_cell *
          shell_energy_density;
        double db = sqrt(cell_energy_density);

        double cos_theta = kz / k;
        double sin_theta = kxy / k;

        double cos_phi = kx / kxy;
        double sin_phi = ky / kxy;

        if (kxy == 0.0) {
          // shouldn't matter given random phase and polarization, just want to
          // avoid nans
          cos_phi = 1.0;
          sin_phi = 0.0;
        }

        double phase = rng.get(0.0, 2.0 * M_PI);
        double polarization = rng.get(0.0, 2.0 * M_PI);

        double aax = db * cos(polarization) * cos_theta * cos_phi;
        double aay = db * cos(polarization) * cos_theta * sin_phi;
        double aaz = -db * cos(polarization) * sin_theta;

        double bbx = db * sin(polarization) * sin_phi;
        double bby = -db * sin(polarization) * cos_phi;
        double bbz = 0.0;

        Double3 k_vec = {kx, ky, kz};

        for (int p = 0; p < mflds.n_patches(); ++p) {
          auto& patch = grid.patches[p];
          auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

          int n_ghosts = mflds.ibn().max();
          grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
            Int3 index{jx, jy, jz};

            Double3 pos_fc_x =
              centering::get_pos(patch, index, centering::FC, 0);
            Double3 pos_fc_y =
              centering::get_pos(patch, index, centering::FC, 1);
            Double3 pos_fc_z =
              centering::get_pos(patch, index, centering::FC, 2);

            double arg_x = k_vec.dot(pos_fc_x) + phase;
            double arg_y = k_vec.dot(pos_fc_y) + phase;
            double arg_z = k_vec.dot(pos_fc_z) + phase;

            field_patch(HX, jx, jy, jz) += aax * cos(arg_x) + bbx * sin(arg_x);
            field_patch(HY, jx, jy, jz) += aay * cos(arg_y) + bby * sin(arg_y);
            field_patch(HZ, jx, jy, jz) += aaz * cos(arg_z) + bbz * sin(arg_z);
          });
        }
      }
    }
  }
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
