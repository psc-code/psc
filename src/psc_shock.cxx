#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "OutputFieldsDefault.h"
#include "psc_config.hxx"
#include "input_params.hxx"
#include "kg/include/kg/VecRange.hxx"

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
int marder_interval;

bool mirror_domain;
std::string turb_method;

int nicell;
int seed;

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
  InputParams inputParams(path_to_params);

  psc_params.stats_every = 1000;
  psc_params.cfl = inputParams.getOrDefault<double>("cfl", .75);
  psc_params.write_checkpoint_every_step = 0;

  electron_temperature = inputParams.get<double>("electron_temperature");
  ion_temperature = inputParams.get<double>("ion_temperature");
  electron_mass = inputParams.get<double>("electron_mass");
  ion_mass = inputParams.get<double>("ion_mass");

  v_upstream_x = inputParams.get<double>("v_upstream_x");
  v_upstream_y = inputParams.get<double>("v_upstream_y");
  v_upstream_z = inputParams.get<double>("v_upstream_z");

  double b_angle_y_to_x_rad = inputParams.get<double>("b_angle_y_to_x_rad");
  double b_mag = inputParams.get<double>("b_mag");
  b_x = b_mag * sin(b_angle_y_to_x_rad);
  b_y = b_mag * cos(b_angle_y_to_x_rad);
  b_z = 0.0;

  nx = inputParams.get<int>("nx");
  ny = inputParams.get<int>("ny");
  nz = inputParams.get<int>("nz");
  psc_params.nmax = inputParams.get<int>("nt");

  n_patches_x = inputParams.get<int>("n_patches_x");
  n_patches_y = inputParams.get<int>("n_patches_y");
  n_patches_z = inputParams.get<int>("n_patches_z");

  double dx = inputParams.get<double>("dx");
  double dy = inputParams.get<double>("dy");
  double dz = inputParams.get<double>("dz");

  len_x = nx * dx;
  len_y = ny * dy;
  len_z = nz * dz;

  if (inputParams.warnIfPresent("turb_dB^2", "set turb_dB instead")) {
    turb_db2 = inputParams.get<double>("turb_dB^2");
  } else {
    turb_db2 = sqr(inputParams.get<double>("turb_dB"));
  }
  turb_correlation_length = inputParams.get<double>("turb_correlation_length");

  int n_writes = inputParams.getOrDefault<int>("n_writes", 100);
  out_interval = psc_params.nmax / n_writes;
  marder_interval = inputParams.getOrDefault<int>("marder_interval", -1);

  mirror_domain = inputParams.getOrDefault<bool>("mirror_domain", false);
  turb_method =
    inputParams.getOrDefault<std::string>("turb_method", "alfven_dense");

  nicell = inputParams.getOrDefault<int>("nicell", 100);
  seed = inputParams.getOrDefault<int>("seed", 5);

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
  norm_params.nicell = nicell;

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

void add_background_fields(MfieldsState& mflds)
{
  const auto& grid = mflds.grid();

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto& patch = grid.patches[p];
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

    int n_ghosts = mflds.ibn().max();
    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      field_patch(HX, jx, jy, jz) += b_x;
      field_patch(HY, jx, jy, jz) += b_y;
      field_patch(HZ, jx, jy, jz) += b_z;
    });
  }
}

void boost_fields(MfieldsState& mflds, double vy)
{
  // general:

  // E_y' = E_y
  // B_y' = B_y
  // E_perp' = gamma * (E_perp + v x B)
  // B_perp' = gamma * (B_perp - v x E)

  // assume E = 0 initially:
  // E_perp' = gamma * v x B
  // B_perp' = gamma * B_perp

  const auto& grid = mflds.grid();
  double gamma = 1.0 / std::sqrt(1.0 - vy * vy);
  Double3 gamma_v{0.0, gamma * vy, 0.0};

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto& patch = grid.patches[p];
    auto mf_patch = make_Fields3d<dim_xyz>(mflds[p]);

    int n_ghosts = mflds.ibn().max();

    grid.Foreach_3d(n_ghosts - 1, n_ghosts, [&](int jx, int jy, int jz) {
      mf_patch(EX, jx, jy, jz) =
        gamma * vy * 0.5 *
        (mf_patch(HZ, jx, jy - 1, jz) + mf_patch(HZ, jx, jy, jz));
      mf_patch(EZ, jx, jy, jz) =
        -gamma * vy * 0.5 *
        (mf_patch(HX, jx, jy - 1, jz) + mf_patch(HX, jx, jy, jz));
    });

    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      mf_patch(HX, jx, jy, jz) *= gamma;
      mf_patch(HZ, jx, jy, jz) *= gamma;
    });
  }
}

void set_mean_b2(MfieldsState& mflds, double mean_b2)
{
  const auto& grid = mflds.grid();

  double sum_b2_local = 0.0;

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

    grid.Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
      // magnetic energy density is cell-centered (H^2)
      // note: this is not the most efficient calculation method
      sum_b2_local += 0.5 * (sqr(field_patch(HX, jx, jy, jz)) +
                             sqr(field_patch(HX, jx + 1, jy, jz)) +
                             sqr(field_patch(HY, jx, jy, jz)) +
                             sqr(field_patch(HY, jx, jy + 1, jz)) +
                             sqr(field_patch(HZ, jx, jy, jz)) +
                             sqr(field_patch(HZ, jx, jy, jz + 1)));
      ;
    });
  }

  double sum_b2 = 0.0;
  MPI_Allreduce(&sum_b2_local, &sum_b2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double target_sum_b2 = nx * ny * nz * mean_b2;
  double scale_factor = sqrt(target_sum_b2 / sum_b2);

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

    int n_ghosts = mflds.ibn().max();
    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      field_patch(HX, jx, jy, jz) *= scale_factor;
      field_patch(HY, jx, jy, jz) *= scale_factor;
      field_patch(HZ, jx, jy, jz) *= scale_factor;
    });
  }
}

int AX = 0, AY = 1, AZ = 2;

void inject_b_from_potential(MfieldsState& mflds,
                             PscConfig::Mfields& vector_potential)
{
  const auto& grid = mflds.grid();

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);
    auto vector_potential_patch = make_Fields3d<dim_xyz>(vector_potential[p]);

    grid.Foreach_3d(2, 1, [&](int jx, int jy, int jz) {
      field_patch(HX, jx, jy, jz) = vector_potential_patch(AZ, jx, jy + 1, jz) -
                                    vector_potential_patch(AZ, jx, jy, jz) -
                                    vector_potential_patch(AY, jx, jy, jz + 1) +
                                    vector_potential_patch(AY, jx, jy, jz);
      field_patch(HY, jx, jy, jz) = vector_potential_patch(AX, jx, jy, jz + 1) -
                                    vector_potential_patch(AX, jx, jy, jz) -
                                    vector_potential_patch(AZ, jx + 1, jy, jz) +
                                    vector_potential_patch(AZ, jx, jy, jz);
      field_patch(HZ, jx, jy, jz) = vector_potential_patch(AY, jx + 1, jy, jz) -
                                    vector_potential_patch(AY, jx, jy, jz) -
                                    vector_potential_patch(AX, jx, jy + 1, jz) +
                                    vector_potential_patch(AX, jx, jy, jz);
    });
  }
}

void inject_plane_alfven_wave(PscConfig::Mfields& vector_potential, double db,
                              Double3& k_vec, double polarization, double phase)
{
  const auto& grid = vector_potential.grid();

  double k2 = k_vec.mag2();
  double k = std::sqrt(k2);
  double kxy = std::sqrt(sqr(k_vec[0]) + sqr(k_vec[1]));

  double cos_theta = k_vec[2] / k;
  double sin_theta = kxy / k;

  double cos_phi = k_vec[0] / kxy;
  double sin_phi = k_vec[1] / kxy;

  if (kxy == 0.0) {
    // shouldn't matter given random phase and polarization, just want to
    // avoid nans
    cos_phi = 1.0;
    sin_phi = 0.0;
  }

  Double3 xp_hat{cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta};
  Double3 yp_hat{sin_phi, -cos_phi, 0};

  Double3 a_vec = db * cos(polarization) / k2 * xp_hat.cross(k_vec);
  Double3 b_vec = db * sin(polarization) / k2 * yp_hat.cross(k_vec);

  for (int p = 0; p < vector_potential.n_patches(); ++p) {
    auto& patch = grid.patches[p];
    auto vector_potential_patch = make_Fields3d<dim_xyz>(vector_potential[p]);

    int n_ghosts = vector_potential.ibn().max();
    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      Int3 index{jx, jy, jz};

      Double3 pos_ec_x = centering::get_pos(patch, index, centering::EC, 0);
      Double3 pos_ec_y = centering::get_pos(patch, index, centering::EC, 1);
      Double3 pos_ec_z = centering::get_pos(patch, index, centering::EC, 2);

      double arg_x = k_vec.dot(pos_ec_x) + phase;
      double arg_y = k_vec.dot(pos_ec_y) + phase;
      double arg_z = k_vec.dot(pos_ec_z) + phase;

      vector_potential_patch(AX, jx, jy, jz) +=
        a_vec[0] * sin(arg_x) - b_vec[0] * cos(arg_x);
      vector_potential_patch(AY, jx, jy, jz) +=
        a_vec[1] * sin(arg_y) - b_vec[1] * cos(arg_y);
      vector_potential_patch(AZ, jx, jy, jz) +=
        a_vec[2] * sin(arg_z) - b_vec[2] * cos(arg_z);
    });
  }
}

struct GaborKernel
{
  void randomize(rng::Uniform<double>& uni)
  {
    double theta = std::acos(uni.get(-1.0, 1.0));
    double phi = uni.get(0.0, 2.0 * M_PI);
    wave_dir = Double3{std::cos(phi) * std::sin(theta),
                       std::sin(phi) * std::sin(theta), std::cos(theta)};

    double wave_phase = uni.get(0.0, 2.0 * M_PI);
  }

  double sample(Double3 x)
  {
    return std::exp(-M_PI * x.mag2()) *
           std::cos(2.0 * M_PI * x.dot(wave_dir) + wave_phase);
  }

  Double3 wave_dir;
  double wave_phase;
};

using Int2 = kg::Vec<int, 2>;
using Int5 = kg::Vec<int, 5>;

int get_cell_seed(Int5 cell_hyperidx, Int5 hyperdims)
{
  while (cell_hyperidx.min() < 0) {
    cell_hyperidx += hyperdims;
  }
  int cell_seed = flatten_index(cell_hyperidx % hyperdims, hyperdims) + seed;
  if (cell_seed == 0) {
    cell_seed = 1;
  }
  return cell_seed;
}

template <typename KERNEL>
double sample_noise_cell(Int5 cell_hyperidx, Int5 hyperdims,
                         Double3 pos_in_cell, KERNEL& kernel)
{
  rng::Uniform uniform{0.0, 1.0, get_cell_seed(cell_hyperidx, hyperdims)};

  // TODO randomize number of impulses
  int n_impulses = 10;

  double sum = 0.0;
  for (int i = 0; i < n_impulses; i++) {
    Double3 kernel_pos_in_cell{uniform.get(0.0, 1.0), uniform.get(0.0, 1.0),
                               uniform.get(0.0, 1.0)};
    kernel.randomize(uniform);
    sum += kernel.sample(pos_in_cell - kernel_pos_in_cell);
    // double weight = uniform.get(-1.0, 1.0);
    // sum += weight * kernel.sample(pos_in_cell - kernel_pos_in_cell);
  }
  return sum;
}

template <typename KERNEL>
double sample_noise(Int2 d_igrid, Double3 pos, Int5 hyperdims, KERNEL& kernel)
{
  Int3 i3 = pos.fint();
  double sum = 0.0;
  for (Int3 cell_idx : VecRange(i3 - 1, i3 + 2)) {
    Int5 hyperindex{d_igrid[0], d_igrid[1], cell_idx[0], cell_idx[1],
                    cell_idx[2]};
    sum +=
      sample_noise_cell(hyperindex, hyperdims, pos - Double3(cell_idx), kernel);
  }
  return sum;
}

template <typename KERNEL>
void inject_turbulence_noise(MfieldsState& mflds, KERNEL& kernel)
{
  const auto& grid = mflds.grid();
  Int3 ibn = mflds.ibn();
  Int3 gdims(grid.domain.gdims);

  PscConfig::Mfields vector_potential{grid, 3, ibn};

  int n_grids = gdims.max() / 2;
  int n_dims = 3;

  for (int i_grid = 0; i_grid < n_grids; i_grid++) {
    int grid_size = i_grid + 1;
    Int3 grid_dims = Int3::min({grid_size, grid_size, grid_size}, gdims);
    Double3 pos_scale = Double3(grid_dims) / grid.domain.length;

    double shell_power =
      1.0 / (1.0 + pow(turb_correlation_length * pos_scale.max() * 2.0 * M_PI,
                       5.0 / 3.0));

    for (int p = 0; p < grid.n_patches(); ++p) {
      auto vp_patch = make_Fields3d<dim_xyz>(vector_potential[p]);
      auto& patch = grid.patches[p];

      for (Int3 i3 : VecRange(-ibn, grid.ldims + ibn)) {
        Double3 pos_ec_x = centering::get_pos(patch, i3, centering::EC, 0);
        Double3 pos_ec_y = centering::get_pos(patch, i3, centering::EC, 1);
        Double3 pos_ec_z = centering::get_pos(patch, i3, centering::EC, 2);

        Int5 hyperdims{n_dims, n_grids, grid_dims[0], grid_dims[1],
                       grid_dims[2]};

        vp_patch(AX, i3[0], i3[1], i3[2]) +=
          shell_power *
          sample_noise({AX, i_grid}, pos_ec_x * pos_scale, hyperdims, kernel);
        vp_patch(AY, i3[0], i3[1], i3[2]) +=
          shell_power *
          sample_noise({AY, i_grid}, pos_ec_y * pos_scale, hyperdims, kernel);
        vp_patch(AZ, i3[0], i3[1], i3[2]) +=
          shell_power *
          sample_noise({AZ, i_grid}, pos_ec_z * pos_scale, hyperdims, kernel);
      }
    }
  }

  inject_b_from_potential(mflds, vector_potential);

  set_mean_b2(mflds, turb_db2);
}

void inject_turbulence_dense(MfieldsState& mflds)
{

  // literally the power of each shell, such that
  // 1/2 <dB^2> = int_{k=0}^\infty shellpower(k) dk.
  // Note that regardless of dimensionality, shellpower(k) ~ k^{-5/3},
  // i.e., shellpower(k) includes the Jacobian.
  auto shell_power = [&](double k) {
    return 1.0 / (1.0 + pow(k * turb_correlation_length, 5.0 / 3.0));
  };

  PscConfig::Mfields vector_potential{mflds._grid(), 3, mflds.ibn()};

  auto inject_at_n = [&](Int3 n) { return n != Int3{0, 0, 0}; };

  // 1. compute values of |k|

  Double3 len_vec{len_x, len_y, len_z};
  Int3 n_vec{nx, ny, nz};

  Double3 dk_vec = 2.0 * M_PI / len_vec;

  Int3 i3_min = (1 - n_vec) / 2;
  Int3 i3_max = n_vec / 2;

  // inject in only half of k-space, since +k and -k modes are indistinguishable
  for (int d = 0; d < 3; d++) {
    if (n_vec[d] > 2) {
      i3_min[d] = 0;
      break;
    }
  }

  double dk = dk_vec.min();
  double kmax = ((Double3)Int3::max(-i3_min, i3_max) * dk_vec).mag();
  int nk = kmax / dk + 1;

  int rank;
  MPI_Comm_rank(mflds.grid().comm(), &rank);
  if (rank == 0) {
    LOG_INFO("nk = %d, kmax = %f\n", nk, kmax);
  }

  auto k_mins = std::vector<double>(nk);
  for (int i = 0; i < nk; i++) {
    k_mins[i] = dk * i;
  }

  // 2. compute number of degeneracies of each k-shell

  auto n_cells_per_shell = std::vector<int>(nk);

  for (Int3 i3 : VecRange(i3_min, i3_max)) {
    if (!inject_at_n(i3)) {
      continue;
    }

    Double3 k_vec = Double3(i3) * dk_vec;
    double k = k_vec.mag();
    int idx = k / dk;

    if (idx >= nk) {
      LOG_ERROR("k=%f > kmax at idx=(%d, %d, %d)\n", k, i3[0], i3[1], i3[2]);
    }

    n_cells_per_shell[idx] += 1;
  }

  // 3. compute powers for each shell

  auto shell_db2s = std::vector<double>(nk);

  for (int i = 0; i < nk; i++) {
    if (n_cells_per_shell[i] > 0) {
      double k_shell = k_mins[i] + 0.5 * dk;
      shell_db2s[i] = shell_power(k_shell) * dk;
    }
  }

  // 4. inject each mode at a random phase and polarization

  // all processes must use same seed to ensure B is continuous
  auto rng = rng::Uniform<double>(0.0, 1.0, seed);
  const auto& grid = mflds.grid();

  for (Int3 i3 : VecRange(i3_min, i3_max)) {
    if (!inject_at_n(i3)) {
      continue;
    }

    Double3 k_vec = Double3(i3) * dk_vec;
    int idx = k_vec.mag() / dk;

    double proportion_of_shell_db2_in_cell = 1.0 / n_cells_per_shell[idx];
    double shell_db2 = shell_db2s[idx];
    double cell_db2 = proportion_of_shell_db2_in_cell * shell_db2;
    double db = sqrt(cell_db2);

    double phase = rng.get(0.0, 2.0 * M_PI);
    double polarization = rng.get(0.0, 2.0 * M_PI);

    inject_plane_alfven_wave(vector_potential, db, k_vec, polarization, phase);
  }

  // step 5: take curl (or something proportional to it)
  inject_b_from_potential(mflds, vector_potential);

  // step 6: normalize db2
  set_mean_b2(mflds, turb_db2);
}

void initializeFields(MfieldsState& mflds)
{
  if (turb_db2 > 0.0) {
    if (turb_method == "alfven_dense") {
      inject_turbulence_dense(mflds);
    } else if (turb_method == "noise_gabor") {
      GaborKernel kernel;
      inject_turbulence_noise(mflds, kernel);
    } else {
      LOG_ERROR("Unrecognized turbulence method: %s\n", turb_method.c_str());
    }
  }

  add_background_fields(mflds);
  boost_fields(mflds, -v_upstream_y);
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
  psc_params.marder_interval = marder_interval;
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
                                          balance, collision, checks);

  psc.add_gauss_corrector(&marder);

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
