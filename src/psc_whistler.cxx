
#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_xyz;
using PscConfig = PscConfig1vbecSingle<Dim>;

// ----------------------------------------------------------------------

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

// ======================================================================
// PscWhistlerParams

struct PscWhistlerParams
{
  double mi_over_me;
  double vA_over_c;
  double amplitude;
  double beta_e_par;
  double beta_i_par;
  double Ti_perp_over_Ti_par;
  double Te_perp_over_Te_par;

  // calculated from the above
  double B0;
  double Te_par;
  double Te_perp;
  double Ti_par;
  double Ti_perp;
  double mi;
  double me;
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
PscWhistlerParams g;

std::string read_checkpoint_filename;

// This is a set of generic PSC params (see include/psc.hxx),
// like number of steps to run, etc, which also should be set by the case
PscParams psc_params;

} // namespace

// ======================================================================
// setupParameters

void setupParameters()
{
  auto comm = MPI_COMM_WORLD;

  // -- set some generic PSC parameters
  psc_params.nmax = 16001;
  psc_params.cfl = 0.98;

  // -- start from checkpoint:
  //
  // Uncomment when wanting to start from a checkpoint, ie.,
  // instead of setting up grid, particles and state fields here,
  // they'll be read from a file
  // FIXME: This parameter would be a good candidate to be provided
  // on the command line, rather than requiring recompilation when change.

  // read_checkpoint_filename = "checkpoint_500.bp";

  // -- Set some parameters specific to this case
  g.mi_over_me = 10.;
  g.vA_over_c = .1;
  g.amplitude = .5;
  g.beta_e_par = .1;
  g.beta_i_par = .1;
  g.Ti_perp_over_Ti_par = 1.;
  g.Te_perp_over_Te_par = 1.;

  // calculate derived paramters
  g.B0 = g.vA_over_c;
  g.Te_par = g.beta_e_par * sqr(g.B0) / 2.;
  g.Te_perp = g.Te_perp_over_Te_par * g.Te_par;
  g.Ti_par = g.beta_i_par * sqr(g.B0) / 2.;
  g.Ti_perp = g.Ti_perp_over_Ti_par * g.Ti_par;
  g.mi = 1.;
  g.me = 1. / g.mi_over_me;

  mpi_printf(comm, "d_i = 1., d_e = %g\n", sqrt(g.me));
  mpi_printf(comm, "om_ci = %g, om_ce = %g\n", g.B0, g.B0 / g.me);
  mpi_printf(comm, "\n");
  mpi_printf(comm, "v_i,perp = %g [c] T_i,perp = %g\n", sqrt(2. * g.Ti_perp),
             g.Ti_perp);
  mpi_printf(comm, "v_i,par  = %g [c] T_i,par = %g\n", sqrt(2. * g.Ti_par),
             g.Ti_par);
  mpi_printf(comm, "v_e,perp = %g [c] T_e,perp = %g\n",
             sqrt(2 * g.Te_perp / g.me), g.Te_perp);
  mpi_printf(comm, "v_e,par  = %g [c] T_e,par = %g\n",
             sqrt(2. * g.Te_par / g.me), g.Te_par);
  mpi_printf(comm, "\n");
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
  // -- setup particle kinds
  Grid_t::Kinds kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1., g.me, "e"};
  kinds[KIND_ION] = {1., g.mi, "i"};

  // --- setup domain
  Grid_t::Real3 LL = {5., 5.,
                      100.}; // domain size (normalized units, ie, in d_i)
  Int3 gdims = {8, 8, 200};  // global number of grid points
  Int3 np = {1, 1, 25};      // division into patches

  auto domain = Grid_t::Domain{gdims, LL, {}, np};

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  // --- generic setup
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

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts)
{
  SetupParticles<Mparticles> setup_particles(*grid_ptr);

  double kz = (4. * M_PI / grid_ptr->domain.length[2]);
  double kperp = (2. * M_PI / grid_ptr->domain.length[0]);

  partitionAndSetupParticles(
    setup_particles, balance, grid_ptr, mprts,
    [&](int kind, double crd[3], psc_particle_npt& npt) {
      double x = crd[0], y = crd[1], z = crd[2];

      double envelope1 = exp(-(z - 25.) * (z - 25.) / 40.);
      double envelope2 = exp(-(z - 75.) * (z - 75.) / 40.);

      switch (kind) {
        case KIND_ELECTRON:
          npt.n = 1.;

          npt.T[0] = g.Te_perp;
          npt.T[1] = g.Te_perp;
          npt.T[2] = g.Te_par;

          // Set velocities for first wave:
          npt.p[0] = -g.amplitude * envelope1 * g.B0 * cos(kperp * y + kz * z);
          npt.p[2] = 0.;

          // Set velocities for second wave:
          npt.p[1] = g.amplitude * envelope2 * g.B0 * sin(kperp * x - kz * z);
          break;
        case KIND_ION:
          npt.n = 1.;

          npt.T[0] = g.Ti_perp;
          npt.T[1] = g.Ti_perp;
          npt.T[2] = g.Ti_par;

          // Set velocities for first wave:
          npt.p[0] = -g.amplitude * envelope1 * g.B0 * cos(kperp * y + kz * z);
          npt.p[2] = 0.;

          // Set velocities for second wave:
          npt.p[1] = g.amplitude * envelope2 * g.B0 * sin(kperp * x - kz * z);
          break;
        default: assert(0);
      }
    });
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  const auto& grid = mflds.grid();
  double kz = (4. * M_PI / grid.domain.length[2]);
  double kperp = (2. * M_PI / grid.domain.length[0]);

  setupFields(mflds, [&](int m, double crd[3]) {
    double x = crd[0], y = crd[1], z = crd[2];
    double envelope1 = exp(-(z - 25.) * (z - 25.) / 40.);
    double envelope2 = exp(-(z - 75.) * (z - 75.) / 40.);

    switch (m) {
      case HX: return g.amplitude * envelope1 * cos(kperp * y + kz * z) * g.B0;
      case HY: return g.amplitude * envelope2 * sin(kperp * x - kz * z) * g.B0;
      case HZ: return g.B0;
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
  auto comm = MPI_COMM_WORLD;

  mpi_printf(comm, "*** Setting up...\n");

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
  psc_params.balance_interval = 0;
  Balance balance{.1, false};

  // -- Sort
  psc_params.sort_interval = 100;

  // -- Collision
  int collision_interval = 1;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.continuity.check_interval = 50;
  checks_params.continuity.threshold = 1e-5;
  checks_params.continuity.print_max_always = false;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = 0;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output
  //
  // FIXME, this really is too complicated and not very flexible

  // -- output fields
  OutputFieldsParams outf_params{};
  outf_params.fields.pfield.out_interval = 200;
  OutputFields<MfieldsState, Mparticles, Dim> outf{grid, outf_params};

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
  // setup initial conditions

  if (read_checkpoint_filename.empty()) {
    initializeParticles(balance, grid_ptr, mprts);
    initializeFields(mflds);
  } else {
    read_checkpoint(read_checkpoint_filename, *grid_ptr, mprts, mflds);
  }

  // ----------------------------------------------------------------------
  // hand off to PscIntegrator to run the simulation

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

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
