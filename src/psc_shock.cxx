#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "output_fields.hxx"
#include "psc_config.hxx"
#include "input_params.hxx"
#include "kg/include/kg/VecRange.hxx"
#include "libpsc/psc_output_particles/output_particles_adios2_impl.hxx"
#include "libpsc/psc_bnd_fields/advecting.hxx"
#include "libpsc/psc_bnd_fields/radiating.hxx"
#include "libpsc/psc_particle_injectors/dirichlet_boundary_injector.hxx"
#include "libpsc/psc_particle_injectors/neumann_boundary_injector.hxx"
#include "libpsc/psc_turbulence/alfven_turbulence.hxx"
#include "libpsc/psc_turbulence/noise_turbulence.hxx"
#include "libpsc/psc_turbulence/turbulence_io.hxx"
#include "libpsc/axis.hxx"

#include <sys/stat.h>
#include <climits>
#include <cstdlib>

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
using real_t = PscConfig::Mfields::real_t;
using Real3 = Vec3<real_t>;

template <typename PULSE>
using RadiatingBC = psc::bnd::field::Radiating<Dim, MfieldsState, PULSE>;

// ======================================================================
// Global parameters

PscParams psc_params;

double electron_mass;
double ion_mass;

double n_upstream;
Double3 v_upstream;
double te_upstream;
double ti_upstream;
Real3 h0_upstream;
Real3 e0; // e0 is constant

// ----------------------------
// standing shock params
double n_downstream;
Double3 v_downstream;
double te_downstream;
double ti_downstream;
Real3 h0_downstream;

double transition_half_width;
double transition_steepness = 2.0; // at least sqrt(3)~1.74
// ----------------------------

Int3 gdims;
Double3 lengths;
Int3 n_patches;

double turb_db2;
double turb_correlation_length;
// noise_gabor only: 0 spaces the noise levels evenly in k, slower in big boxes
double turb_levels_per_octave;

int out_interval;
int out_moments_interval;
int out_prt_interval;
int checks_interval;
int marder_interval;

std::string turb_method;
std::string shock_method;
std::string injection_method;
std::string neumann_resample_method;
int neumann_copy_depth;

PreaccelerateMethod preaccelerate_method;
double preaccelerate_time_plasma_periods;
double preaccelerate_ion_time_plasma_periods;
bool preaccelerate_sample_boundary;

int nicell;
int seed;

std::string checkpoint_filename;

// ======================================================================
// setupParameters

void setupParameters(int argc, char** argv)
{
  if (argc != 2) {
    LOG_ERROR("Usage: %s path/to/params.txt\n", argv[0]);
  }
  std::string path_to_params(argv[1]);
  InputParams inputParams(path_to_params);

  shock_method = inputParams.getOrDefault<std::string>("shock_method", "wall");
  injection_method =
    inputParams.getOrDefault<std::string>("injection_method", "dirichlet");
  if (injection_method != "dirichlet" && injection_method != "neumann" &&
      injection_method != "periodic") {
    LOG_ERROR("unknown injection_method: %s\n", injection_method.c_str());
  }
  if (injection_method == "periodic" && shock_method == "wall") {
    LOG_ERROR("injection_method=periodic is incompatible with "
              "shock_method=wall\n");
  }
  neumann_resample_method =
    inputParams.getOrDefault<std::string>("neumann_resample_method", "none");
  if (neumann_resample_method != "none" &&
      neumann_resample_method != "maxwellian" &&
      neumann_resample_method != "rotate") {
    LOG_ERROR("unknown neumann_resample_method: %s\n",
              neumann_resample_method.c_str());
  }
  // copy particles as they cross the face this many cells from the boundary
  neumann_copy_depth = inputParams.getOrDefault<int>("neumann_copy_depth", 1);

  psc_params.stats_every =
    inputParams.getOrDefault<int>("stats_interval", 1000);
  psc_params.cfl = inputParams.getOrDefault<double>("cfl", .75);
  checkpoint_filename =
    inputParams.getOrDefault<std::string>("checkpoint_filename", "");

  electron_mass = inputParams.get<double>("m_e");
  ion_mass = inputParams.get<double>("m_i");

  n_upstream = 1.0;
  te_upstream = inputParams.get<double>("T_e");
  ti_upstream = inputParams.get<double>("T_i");

  double theta_bn_deg = inputParams.get<double>("θ_Bn_deg");
  double b0 = inputParams.get<double>("B_0");

  double theta_bn = theta_bn_deg * M_PI / 180.0;
  double theta_xz = inputParams.get<double>("θ_xz_deg") * M_PI / 180.0;
  h0_upstream = Real3{sin(theta_bn) * cos(theta_xz), cos(theta_bn),
                      sin(theta_bn) * sin(theta_xz)} *
                b0;

  v_upstream = {0.0, inputParams.get<double>("v_upstream"), 0.0};
  e0 = -v_upstream.cross(h0_upstream);

  if (shock_method == "wall" || shock_method == "none") {
    // relativistic correction
    double gamma = 1 / sqrt(1 - v_upstream.mag2());
    e0 *= gamma;
    h0_upstream *= Real3{gamma, 1.0, gamma};
  } else if (shock_method == "relaxation") {
    // for perpendicular shock (2013 Balogh eq.3.36 and normalization in
    // sec.3.3.1)
    if (theta_bn_deg != 90.0) {
      LOG_ERROR("θ_Bn must be 90° for relaxation method; got %f°\n",
                theta_bn_deg);
    }

    // no relativistic correction, since these aren't relativistic RH conditions

    double b_norm = sqrt(ion_mass * n_upstream * v_upstream.mag2());
    double t_norm = 0.5 * ion_mass * v_upstream.mag2();
    double beta1 =
      2.0 * n_upstream * (te_upstream + ti_upstream) / h0_upstream.mag2();
    Double3 B1 = h0_upstream / b_norm;
    double T1 = (te_upstream + ti_upstream) / t_norm;

    double MA_sq =
      v_upstream.mag2() * ion_mass * n_upstream / h0_upstream.mag2();

    double tmp = 1.0 + (1.0 + 2.5 * beta1) * B1.mag2();
    double r = 8.0 / (tmp + sqrt(sqr(tmp) + 2 * B1.mag2()));
    r = inputParams.getOrDefault<double>("r", r);

    double heating_factor =
      1.0 +
      4.0 / (5.0 * T1) * ((sqr(r) - 1.0) / (2.0 * sqr(r)) + (1.0 - r) / MA_sq);

    n_downstream = n_upstream * r;
    v_downstream = v_upstream / r;
    h0_downstream = h0_upstream * r;

    te_downstream =
      inputParams.getOrDefault<double>("T_e2", te_upstream * heating_factor);
    ti_downstream =
      inputParams.getOrDefault<double>("T_i2", ti_upstream * heating_factor);
  }

  if (shock_method == "none") {
    n_downstream = n_upstream;
    v_downstream = v_upstream;
    h0_downstream = h0_upstream;
    te_downstream = te_upstream;
    ti_downstream = ti_upstream;
  }

  gdims[0] = inputParams.get<int>("nx");
  gdims[1] = inputParams.get<int>("ny");
  gdims[2] = inputParams.get<int>("nz");
  psc_params.nmax = inputParams.get<int>("nt");

  n_patches[0] = inputParams.get<int>("npx");
  n_patches[1] = inputParams.get<int>("npy");
  n_patches[2] = inputParams.get<int>("npz");

  if (inputParams.has("lx")) {
    lengths[0] = inputParams.get<double>("lx");
  } else {
    lengths[0] = inputParams.get<double>("dx") * gdims[0];
  }
  if (inputParams.has("ly")) {
    lengths[1] = inputParams.get<double>("ly");
  } else {
    lengths[1] = inputParams.get<double>("dy") * gdims[1];
  }
  if (inputParams.has("lz")) {
    lengths[2] = inputParams.get<double>("lz");
  } else {
    lengths[2] = inputParams.get<double>("dz") * gdims[2];
  }

  transition_half_width = lengths[1] / 8.0;

  turb_db2 = sqr(inputParams.get<double>("dB"));
  turb_correlation_length = inputParams.get<double>("L_c");

  if (inputParams.has("checkpoint_interval")) {
    psc_params.write_checkpoint_every_step =
      inputParams.get<int>("checkpoint_interval");
    inputParams.errIfPresent(
      "n_checkpoints",
      "n_checkpoints is mutually exclusive with checkpoint_interval");
  } else if (inputParams.has("n_checkpoints")) {
    int n_checkpoints = inputParams.get<int>("n_checkpoints");
    if (n_checkpoints > 0) {
      psc_params.write_checkpoint_every_step = psc_params.nmax / n_checkpoints;
    }
  }

  if (injection_method == "dirichlet") {
    preaccelerate_time_plasma_periods = inputParams.getOrDefault<double>(
      "preaccelerate_time_plasma_periods", 0.5);
    // ions may be kicked differently (in ion plasma periods); 0 disables
    preaccelerate_ion_time_plasma_periods =
      inputParams.getOrDefault<double>("preaccelerate_ion_time_plasma_periods",
                                       preaccelerate_time_plasma_periods);
    std::string pre_method_str =
      inputParams.getOrDefault<std::string>("preaccelerate_method", "normal_e");
    preaccelerate_method =
      pre_method_str == "normal_e"  ? PreaccelerateMethod::NormalE
      : pre_method_str == "all_e"   ? PreaccelerateMethod::AllE
      : pre_method_str == "all_e_h" ? PreaccelerateMethod::AllEH
                                    : PreaccelerateMethod::None;
    // "cell": at the edge cell's center; "boundary": at the boundary itself
    preaccelerate_sample_boundary =
      inputParams.getOrDefault<std::string>("preaccelerate_sample", "cell") ==
      "boundary";
  }

  int n_writes = inputParams.getOrDefault<int>("n_writes", 100);
  // n_writes_* default to n_writes; 0 disables that output
  auto interval = [&](int n) { return n > 0 ? psc_params.nmax / n : 0; };
  out_interval = interval(n_writes);
  out_moments_interval =
    interval(inputParams.getOrDefault<int>("n_writes_moments", n_writes));
  out_prt_interval =
    interval(inputParams.getOrDefault<int>("n_writes_prt", n_writes));
  checks_interval =
    interval(inputParams.getOrDefault<int>("n_writes_checks", n_writes));
  marder_interval = inputParams.getOrDefault<int>("marder_interval", -1);

  turb_method =
    inputParams.getOrDefault<std::string>("turb_method", "alfven_dense");
  turb_levels_per_octave =
    inputParams.getOrDefault<double>("turb_levels_per_octave", 4.0);

  nicell = inputParams.getOrDefault<int>("nicell", 100);
  seed = inputParams.getOrDefault<int>("seed", 5);
  rng::set_base_seed(seed);

  std::ifstream src(path_to_params, std::ios::binary);
  std::ofstream dst("params_record.txt", std::ios::binary);
  dst << src.rdbuf();
}

template <typename T>
T interpolate_across_shock(T upstream, T downstream, double y)
{
  if (transition_half_width == 0.0) {
    return y > 0.0 ? downstream : upstream;
  } else {
    // smooth transition function
    y /= transition_half_width;
    double sigmoid_val = abs(y) >= 1.0
                           ? (y > 0.0 ? 1.0 : -1.0)
                           : tanh(transition_steepness * y / (1.0 - y * y));
    double weight_downstream = (sigmoid_val + 1.0) / 2.0;
    double weight_upstream = 1.0 - weight_downstream;
    return upstream * weight_upstream + downstream * weight_downstream;
  }
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

  Double3 corner;
  int bnd_fld_lower;
  int bnd_fld_upper;
  int bnd_prt_lower;
  int bnd_prt_upper;

  if (shock_method == "wall") {
    corner = {0.0, 0.0, 0.0};
    bnd_fld_lower = BND_FLD_OPEN;
    bnd_fld_upper = BND_FLD_CONDUCTING_WALL;
    bnd_prt_lower = BND_PRT_OPEN;
    bnd_prt_upper = BND_PRT_REFLECTING;
  } else if (shock_method == "relaxation" || shock_method == "none") {
    corner = {0.0, -lengths[1] / 2.0, 0.0};
    bnd_fld_lower = BND_FLD_OPEN;
    bnd_fld_upper = BND_FLD_OPEN;
    bnd_prt_lower = BND_PRT_OPEN;
    bnd_prt_upper = BND_PRT_OPEN;
  }

  if (injection_method == "periodic") {
    bnd_fld_lower = BND_FLD_PERIODIC;
    bnd_fld_upper = BND_FLD_PERIODIC;
    bnd_prt_lower = BND_PRT_PERIODIC;
    bnd_prt_upper = BND_PRT_PERIODIC;
  }

  auto domain = Grid_t::Domain{gdims, lengths, corner, n_patches};

  auto bc = psc::grid::BC{{BND_FLD_PERIODIC, bnd_fld_lower, BND_FLD_PERIODIC},
                          {BND_FLD_PERIODIC, bnd_fld_upper, BND_FLD_PERIODIC},
                          {BND_PRT_PERIODIC, bnd_prt_lower, BND_PRT_PERIODIC},
                          {BND_PRT_PERIODIC, bnd_prt_upper, BND_PRT_PERIODIC}};

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
  setup_particles.random_offsets = true;
  setup_particles.initial_momentum_gamma_correction = true;

  if (shock_method == "wall" || shock_method == "none") {
    auto init_np = [&](int kind, Double3 pos, int p, Int3 idx,
                       psc_particle_np& np) {
      double t = np.kind == KIND_ION ? ti_upstream : te_upstream;
      np.n = 1.0;
      np.p = setup_particles.createMaxwellian(
        {np.kind, np.n, v_upstream, {t, t, t}, np.tag});
    };

    partitionAndSetupParticles(setup_particles, balance, grid_ptr, mprts,
                               init_np);
  } else if (shock_method == "relaxation") {
    auto init_np = [&](int kind, Double3 pos, int p, Int3 idx,
                       psc_particle_np& np) {
      np.n = interpolate_across_shock(n_upstream, n_downstream, pos[1]);
      // interpolate v_thermal, not T itself
      double t = sqr(interpolate_across_shock(
        sqrt(np.kind == KIND_ION ? ti_upstream : te_upstream),
        sqrt(np.kind == KIND_ION ? ti_downstream : te_downstream), pos[1]));
      Double3 v = interpolate_across_shock(v_upstream, v_downstream, pos[1]);
      np.p =
        setup_particles.createMaxwellian({np.kind, np.n, v, {t, t, t}, np.tag});
    };

    partitionAndSetupParticles(setup_particles, balance, grid_ptr, mprts,
                               init_np);
  }
}

// ======================================================================
// initializeFields

void add_background_fields(MfieldsState& mflds)
{
  const auto& grid = mflds.grid();

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto& patch = grid.patches[p];
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

    int n_ghosts = mflds.ibn().max();
    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      Real3 h0;
      if (shock_method == "wall" || shock_method == "none") {
        h0 = h0_upstream;
      } else if (shock_method == "relaxation") {
        Double3 pos = centering::get_pos(patch, {jx, jy, jz}, centering::NC, 0);
        h0 = interpolate_across_shock(h0_upstream, h0_downstream, pos[1]);
        h0[1] = h0_upstream[1]; // parallel B isn't compressed
      }

      field_patch(HX, jx, jy, jz) += h0[0];
      field_patch(HY, jx, jy, jz) += h0[1];
      field_patch(HZ, jx, jy, jz) += h0[2];

      field_patch(EX, jx, jy, jz) += e0[0];
      field_patch(EY, jx, jy, jz) += e0[1];
      field_patch(EZ, jx, jy, jz) += e0[2];
    });
  }
}

void initialize_turbulence(MfieldsState& mflds)
{
  if (turb_db2 == 0.0) {
    return;
  }

  using namespace psc::turbulence;
  TurbulenceParams params{turb_db2, turb_correlation_length, seed};
  std::unique_ptr<TurbulenceGenerator<MfieldsState>> generator;
  if (turb_method == "alfven_dense") {
    generator.reset(new AlfvenTurbulence<MfieldsState>{params});
  } else if (turb_method == "noise_gabor") {
    auto noise = new NoiseTurbulence<MfieldsState>{params};
    noise->levels_per_octave = turb_levels_per_octave;
    generator.reset(noise);
  } else {
    LOG_ERROR("Unrecognized turbulence method: %s\n", turb_method.c_str());
  }
  generator->generate(mflds);
}

// ----------------------------------------------------------------------
// saving the advected fields for restarts
//
// The lower boundary keeps advecting in the initial turbulence, so a restart
// needs it too. With checkpoints enabled, it is written to advected_fields.bp
// in the run directory; a restart reads it from the checkpoint's directory.

const std::string advected_fields_filename = "advected_fields.bp";

std::string get_dirname(std::string path)
{
  // .bp "files" are directories, so may be given with a trailing slash
  while (path.size() > 1 && path.back() == '/') {
    path.pop_back();
  }
  auto slash = path.find_last_of('/');
  if (slash == std::string::npos) {
    return ".";
  }
  return slash == 0 ? "/" : path.substr(0, slash);
}

bool path_exists(const std::string& path)
{
  struct stat buf;
  return stat(path.c_str(), &buf) == 0;
}

bool is_same_dir(const std::string& a, const std::string& b)
{
  char real_a[PATH_MAX], real_b[PATH_MAX];
  return realpath(a.c_str(), real_a) && realpath(b.c_str(), real_b) &&
         std::string(real_a) == std::string(real_b);
}

bool checkpoints_enabled()
{
  return psc_params.write_checkpoint_every_step > 0;
}

// Restores the advected fields (the initial turbulence) of the run
// `checkpoint_filename` belongs to, or regenerates them if that run didn't
// save them.
void restore_advected_fields(MfieldsState& mflds)
{
  std::string checkpoint_dir = get_dirname(checkpoint_filename);
  std::string path = checkpoint_dir + "/" + advected_fields_filename;
  if (path_exists(path)) {
    mpi_printf(mflds.grid().comm(), "**** Reading advected fields from %s\n",
               path.c_str());
    psc::turbulence::read_turbulence(path, mflds);
  } else {
    mpi_printf(mflds.grid().comm(),
               "WARNING: %s not found; regenerating turbulence, which only "
               "matches the original run if its parameters and seed do\n",
               path.c_str());
    initialize_turbulence(mflds);
  }

  // keep it next to this run's checkpoints, too, for the next restart
  if (checkpoints_enabled() && !is_same_dir(checkpoint_dir, ".")) {
    psc::turbulence::write_turbulence(advected_fields_filename, mflds);
  }
}

// ======================================================================
// makeNeumannInjector

InjectorBase<Mparticles, MfieldsState>* makeNeumannInjector(const Grid_t& grid,
                                                            LoHi lohi,
                                                            double te,
                                                            double ti)
{
  using PushParticles = PscConfig::PushParticles;
  using psc::Axis;

  if (neumann_resample_method == "maxwellian") {
    std::vector<double> temperatures(NR_KINDS);
    temperatures[KIND_ELECTRON] = te;
    temperatures[KIND_ION] = ti;
    return new NeumannBoundaryInjector<PushParticles,
                                       NeumannResamplerMaxwellian>{
      Axis::Y, lohi, NeumannResamplerMaxwellian{grid.kinds, temperatures},
      neumann_copy_depth};
  } else if (neumann_resample_method == "rotate") {
    return new NeumannBoundaryInjector<PushParticles, NeumannResamplerRotate>{
      Axis::Y, lohi, NeumannResamplerRotate{}, neumann_copy_depth};
  }

  return new NeumannBoundaryInjector<PushParticles>{
    Axis::Y, lohi, NeumannResamplerNone{}, neumann_copy_depth};
}

// ======================================================================
// makeDirichletInjector

// injects one species, Maxwellian with temperature `t` drifting at `v`, at
// density `n`, through the given y boundary; preacceleration time is given in
// plasma periods of that species
InjectorBase<Mparticles, MfieldsState>* makeDirichletInjector(
  const Grid_t& grid, LoHi lohi, int kind, double n, Double3 v, double t,
  double mass, double preaccelerate_plasma_periods)
{
  auto injector = new DirichletBoundaryInjector<ParticleGeneratorMaxwellian,
                                                PscConfig::PushParticles>(
    ParticleGeneratorMaxwellian(kind, grid.kinds[kind], v, {t, t, t}),
    psc::Axis::Y, lohi, n);
  injector->preaccelerate_method = preaccelerate_method;
  injector->preaccelerate_sample_boundary = preaccelerate_sample_boundary;
  injector->preaccelerate_time =
    preaccelerate_plasma_periods * 2.0 * M_PI / std::sqrt(n / mass);
  return injector;
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
  psc_params.balance_interval = 1000;
  Balance balance{.1};

  // -- Sort
  psc_params.sort_interval = 100;

  // -- Collision
  int collision_interval = 0;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.gauss.check_interval = checks_interval;
  checks_params.continuity.check_interval = checks_interval;

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
  OutputFields<MfieldsState, Mparticles> out_fields;
  out_fields.pfield.out_interval = out_interval;

  OutputMoments<MfieldsState, Mparticles, Dim> out_moments;
  out_moments.pfield.out_interval = out_moments_interval;

  // -- output particles
  OutputParticlesAdios2Params outp_params{};
  outp_params.every_step = out_prt_interval;
  OutputParticlesAdios2<Mparticles, float> outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies<Mparticles, MfieldsState> oute{grid.comm(), oute_interval};

  // ----------------------------------------------------------------------
  // set up initial conditions

  // the lower boundary injects the initial turbulence, advected in
  bool advect_turbulence =
    injection_method != "periodic" && turb_db2 > 0.0 && v_upstream[1] > 0.0;

  if (checkpoint_filename.empty()) {
    initializeParticles(balance, grid_ptr, mprts);
    initialize_turbulence(mflds);
    if (advect_turbulence && checkpoints_enabled()) {
      // before the background fields are added
      psc::turbulence::write_turbulence(advected_fields_filename, mflds);
    }
  } else {
    read_checkpoint(checkpoint_filename, *grid_ptr, mprts, mflds);
  }

  // ----------------------------------------------------------------------
  // run the simulation

  auto psc = makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts,
                                          balance, collision, checks);

  psc.add_gauss_corrector(&marder);

  psc.add_diagnostic(&out_fields);
  psc.add_diagnostic(&out_moments);
  psc.add_diagnostic(&outp);
  psc.add_diagnostic(&oute);

  using psc::Axis;
  using psc::bnd::LoHi;
  using ConstantPulse = psc::bnd::field::ConstantPulse<real_t>;

  if (injection_method == "dirichlet") {
    psc.add_injector(makeDirichletInjector(
      mprts.grid(), LoHi::Lo, KIND_ION, n_upstream, v_upstream, ti_upstream,
      ion_mass, preaccelerate_ion_time_plasma_periods));
    psc.add_injector(makeDirichletInjector(
      mprts.grid(), LoHi::Lo, KIND_ELECTRON, n_upstream, v_upstream,
      te_upstream, electron_mass, preaccelerate_time_plasma_periods));
  } else if (injection_method == "neumann") {
    psc.add_injector(
      makeNeumannInjector(mprts.grid(), LoHi::Lo, te_upstream, ti_upstream));
  }

  if (injection_method != "periodic") {
    if (advect_turbulence) {
      // mflds is currently just the pure, initial turbulence
      MfieldsState* mflds_advect = &mflds;
      // the pulse copies it, so a restored one need only outlive this block
      std::unique_ptr<MfieldsState> mflds_restored;
      if (!checkpoint_filename.empty()) {
        // mflds is completely unrelated; need the initial turbulence
        mflds_restored.reset(new MfieldsState{mprts.grid()});
        restore_advected_fields(*mflds_restored);
        mflds_advect = mflds_restored.get();
      }

      using AdvectingPulse =
        psc::bnd::field::AdvectingPulse<Dim, PscConfig::Mfields>;
      auto radiating_lo = new RadiatingBC<AdvectingPulse>(
        AdvectingPulse{*mflds_advect, Axis::Y, real_t(v_upstream[1]), e0,
                       h0_upstream},
        Axis::Y, LoHi::Lo);

      psc.add_field_bc(radiating_lo);
    } else {
      psc.add_field_bc(new RadiatingBC<ConstantPulse>(
        ConstantPulse{e0, h0_upstream}, Axis::Y, LoHi::Lo));
    }
  }

  if ((shock_method == "relaxation" || shock_method == "none") &&
      injection_method != "periodic") {
    if (injection_method == "dirichlet") {
      psc.add_injector(makeDirichletInjector(
        mprts.grid(), LoHi::Hi, KIND_ION, n_downstream, v_downstream,
        ti_downstream, ion_mass, preaccelerate_ion_time_plasma_periods));
      psc.add_injector(makeDirichletInjector(
        mprts.grid(), LoHi::Hi, KIND_ELECTRON, n_downstream, v_downstream,
        te_downstream, electron_mass, preaccelerate_time_plasma_periods));
    } else if (injection_method == "neumann") {
      psc.add_injector(makeNeumannInjector(mprts.grid(), LoHi::Hi,
                                           te_downstream, ti_downstream));
    }

    auto radiating_hi = new RadiatingBC<ConstantPulse>(
      ConstantPulse{e0, h0_downstream}, Axis::Y, LoHi::Hi);

    psc.add_field_bc(radiating_hi);
  }

  if (checkpoint_filename.empty()) {
    // add background after initializing radiation inflow, which only wants
    // the perturbations to B
    add_background_fields(mflds);
  }

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
