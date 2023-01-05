
#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"

#include "psc_bgk_util/bgk_params.hxx"
#include "psc_bgk_util/input_parser.hxx"
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

// for parser
enum DATA_COL
{
  COL_RHO,
  COL_NE,
  COL_V_PHI,
  COL_TE,
  COL_E_RHO,
  COL_PHI,
  n_cols
};

// ======================================================================
// Global parameters

namespace
{
ParsedData* parsedData;

std::string read_checkpoint_filename;

// Parameters specific to BGK simulation
PscBgkParams g;

// General PSC parameters (see include/psc.hxx),
PscParams psc_params;

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
  g.loadParams(parsedParams);

  parsedData = new ParsedData(parsedParams.get<std::string>("path_to_data"),
                              n_cols, COL_RHO, 1);

  psc_params.nmax = parsedParams.get<int>("nmax");
  psc_params.stats_every = parsedParams.get<int>("stats_every");
  psc_params.cfl = parsedParams.getOrDefault<double>("cfl", .75);

  psc_params.write_checkpoint_every_step =
    parsedParams.getOrDefault<int>("checkpoint_every", 0);
  if (parsedParams.getOrDefault<bool>("read_checkpoint", false))
    read_checkpoint_filename =
      parsedParams.get<std::string>("path_to_checkpoint");

  std::ifstream src(path_to_params, std::ios::binary);
  std::ofstream dst("params_record.txt", std::ios::binary);
  dst << src.rdbuf();

  if (g.n_grid_3 > 1 && typeid(Dim) != typeid(dim_xyz)) {
    LOG_ERROR("3D runs require Dim = dim_xyz\n");
    exit(1);
  } else if (g.n_grid == 1 && typeid(Dim) != typeid(dim_yz)) {
    LOG_ERROR("2D runs require Dim = dim_yz\n");
    exit(1);
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
  auto domain = Grid_t::Domain{
    {g.n_grid_3, g.n_grid, g.n_grid},           // # grid points
    {g.box_size_3, g.box_size, g.box_size},     // physical lengths
    {0., -.5 * g.box_size, -.5 * g.box_size},   // *offset* for origin
    {g.n_patches_3, g.n_patches, g.n_patches}}; // # patches

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {g.q_e, g.m_e, "e"};
  kinds[KIND_ION] = {g.q_i, g.m_i, "i"};

  mpi_printf(MPI_COMM_WORLD, "lambda_D = %g\n",
             sqrt(parsedData->get_interpolated(COL_TE, g.box_size / sqrt(2))));

  // --- generic setup
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = g.nicell;

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
// writeGT

template <typename GT>
void writeGT(const GT& gt, const Grid_t& grid, const std::string& name,
             const std::vector<std::string>& compNames)
{
  WriterMRC writer;
  writer.open(name);
  writer.begin_step(grid.timestep(), grid.timestep() * grid.dt);
  writer.write(gt, grid, name, compNames);
  writer.end_step();
  writer.close();
}

// ----------------------------------------------------------------------
// writeMF

template <typename MF>
void writeMF(MF&& mfld, const std::string& name,
             const std::vector<std::string>& compNames)
{
  writeGT(psc::mflds::interior(mfld.grid(), mfld.gt()), mfld.grid(), name,
          compNames);
}

// ======================================================================
// helper methods

inline double getCoord(double crd)
{
  if (crd < -g.box_size / 2)
    return crd + g.box_size;
  if (crd > g.box_size / 2)
    return crd - g.box_size;
  return crd;
}

template <int LEN>
inline void setAll(double (&vals)[LEN], double newval)
{
  for (int i = 0; i < LEN; i++)
    vals[i] = newval;
}

inline double getIonDensity(double rho)
{
  if (!g.do_ion)
    return g.n_i;
  double potential = parsedData->get_interpolated(COL_PHI, rho);
  return std::exp(-potential / g.T_i);
}

// ======================================================================
// initializeParticles

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts,
                         BgkMfields& divGradPhi)
{
  SetupParticles<Mparticles> setup_particles(*grid_ptr);
  setup_particles.centerer = Centering::Centerer(Centering::NC);

  auto&& qDensity = -psc::mflds::interior(divGradPhi.grid(), divGradPhi.gt());

  auto npt_init = [&](int kind, double crd[3], int p, Int3 idx,
                      psc_particle_npt& npt) {
    double y = getCoord(crd[1]);
    double z = getCoord(crd[2]);
    double rho = sqrt(sqr(y) + sqr(z));

    switch (kind) {

      case KIND_ELECTRON:
        npt.n = (qDensity(idx[0], idx[1], idx[2], 0, p) -
                 getIonDensity(rho) * g.q_i) /
                g.q_e;
        if (rho != 0) {
          double v_phi = parsedData->get_interpolated(COL_V_PHI, rho);
          double coef = g.v_e_coef * (g.reverse_v ? -1 : 1) *
                        (g.reverse_v_half && y < 0 ? -1 : 1);
          npt.p[0] = 0;
          npt.p[1] = coef * g.m_e * v_phi * -z / rho;
          npt.p[2] = coef * g.m_e * v_phi * y / rho;
        } else {
          setAll(npt.p, 0);
        }
        setAll(npt.T, g.T_e_coef * parsedData->get_interpolated(COL_TE, rho));
        break;

      case KIND_ION:
        npt.n = getIonDensity(rho);
        setAll(npt.p, 0);
        setAll(npt.T, g.T_i);
        break;

      default: assert(false);
    }
  };

  partitionParticlesGeneralInit(setup_particles, balance, grid_ptr, mprts,
                                npt_init);
  setupParticlesGeneralInit(setup_particles, mprts, npt_init);
}

// ======================================================================
// fillGhosts

template <typename MF>
void fillGhosts(MF& mfld, int compBegin, int compEnd)
{
  auto ibn = mfld.ibn();
  int ibn_arr[3] = {ibn[0], ibn[1], ibn[2]};
  Bnd_ bnd{mfld.grid(), ibn_arr};
  bnd.fill_ghosts(mfld, compBegin, compEnd);
}

// ======================================================================
// initializePhi

void initializePhi(BgkMfields& phi)
{
  setupScalarField(
    phi, Centering::Centerer(Centering::NC), [&](int m, double crd[3]) {
      double rho = sqrt(sqr(getCoord(crd[1])) + sqr(getCoord(crd[2])));
      return parsedData->get_interpolated(COL_PHI, rho);
    });

  writeMF(phi, "phi", {"phi"});
}

// ======================================================================
// initializeGradPhi

void initializeGradPhi(BgkMfields& phi, BgkMfields& gradPhi)
{
  auto&& grad = psc::item::grad_ec(phi.gt(), phi.grid());
  psc::mflds::interior(gradPhi.grid(), gradPhi.storage()) = grad;

  fillGhosts(gradPhi, 0, 3);

  writeMF(gradPhi, "grad_phi", {"gradx", "grady", "gradz"});
}

// ======================================================================
// initializeDivGradPhi

void initializeDivGradPhi(BgkMfields& gradPhi, BgkMfields& divGradPhi)
{
  auto&& divGrad = psc::item::div_nc(gradPhi.grid(), gradPhi.gt());
  psc::mflds::interior(divGradPhi.grid(), divGradPhi.storage()) = divGrad;

  fillGhosts(divGradPhi, 0, 1);

  writeMF(divGradPhi, "div_grad_phi", {"divgrad"});
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds, BgkMfields& gradPhi)
{
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case HX: return g.Hx;
      default: return 0.;
    }
  });

  // initialize E separately
  mflds.storage().view(_all, _all, _all, _s(EX, EX + 3)) = -gradPhi.gt();
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
  psc_params.sort_interval = 10;

  // -- Collision
  int collision_interval = 0;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.gauss_every_step = g.gauss_every;
  // checks_params.gauss_dump_always = true;
  checks_params.gauss_threshold = 1e-5;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = 5;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output
  //
  // FIXME, this really is too complicated and not very flexible

  // -- output fields
  OutputFieldsParams outf_params{};
  outf_params.fields.pfield_interval = g.fields_every;
  outf_params.moments.pfield_interval = g.moments_every;
  OutputFields<MfieldsState, Mparticles, Dim> outf{grid, outf_params};

  // -- output particles
  OutputParticlesParams outp_params{};
  outp_params.every_step = g.particles_every;
  outp_params.data_dir = ".";
  outp_params.basename = "prt";
  OutputParticles outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // Set up initial conditions

  if (read_checkpoint_filename.empty()) {
    initializePhi(phi);
    initializeGradPhi(phi, gradPhi);
    initializeDivGradPhi(gradPhi, divGradPhi);
    initializeParticles(balance, grid_ptr, mprts, divGradPhi);
    initializeFields(mflds, gradPhi);
  } else {
    read_checkpoint(read_checkpoint_filename, *grid_ptr, mprts, mflds);
  }

  delete parsedData;

  // ----------------------------------------------------------------------
  // Hand off to PscIntegrator to run the simulation

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
