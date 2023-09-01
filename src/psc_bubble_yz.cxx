
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

using Dim = dim_yz;
#ifdef USE_CUDA
using PscConfig = PscConfig1vbecCuda<Dim>;
#else
using PscConfig = PscConfig1vbecSingle<Dim>;
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
// PscBubbleParams

struct PscBubbleParams
{
  double BB;
  double nnb;
  double nn0;
  double MMach;
  double LLn;
  double LLB;
  double LLz;
  double LLy;
  double TTe;
  double TTi;
  double MMi;
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
PscBubbleParams g;

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
  psc_params.nmax = 1000; // 32000;
  psc_params.stats_every = 100;

  // -- start from checkpoint:
  //
  // Uncomment when wanting to start from a checkpoint, ie.,
  // instead of setting up grid, particles and state fields here,
  // they'll be read from a file
  // FIXME: This parameter would be a good candidate to be provided
  // on the command line, rather than requiring recompilation when change.

  // read_checkpoint_filename = "checkpoint_500.bp";

  // -- Set some parameters specific to this case
  g.BB = .07;
  g.nnb = .1;
  g.nn0 = 1.;
  g.MMach = 3.;
  g.LLn = 200.;
  g.LLB = 200. / 6.;
  g.TTe = .02;
  g.TTi = .02;
  g.MMi = 100.;

  g.LLy = 2. * g.LLn;
  g.LLz = 3. * g.LLn;
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
  auto domain = Grid_t::Domain{{1, 512 * 2, 512 * 3},
                               {g.LLn, g.LLy, g.LLz},
                               {0., -.5 * g.LLy, -.5 * g.LLz},
                               {1, 16 * 2, 16 * 3}};

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1., 1., "e"};
  kinds[KIND_ION] = {1., 100., "i"};

  mpi_printf(MPI_COMM_WORLD, "lambda_D = %g\n", sqrt(g.TTe));

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

  double V0 = g.MMach * sqrt(g.TTe / g.MMi);

  partitionAndSetupParticles(
    setup_particles, balance, grid_ptr, mprts,
    [&](int kind, double crd[3], psc_particle_npt& npt) {
      double r1 = sqrt(sqr(crd[2]) + sqr(crd[1] + .5 * g.LLy));
      double r2 = sqrt(sqr(crd[2]) + sqr(crd[1] - .5 * g.LLy));

      npt.n = g.nnb;
      if (r1 < g.LLn) {
        npt.n += (g.nn0 - g.nnb) * sqr(cos(M_PI / 2. * r1 / g.LLn));
        if (r1 > 0.0) {
          npt.p[2] += V0 * sin(M_PI * r1 / g.LLn) * crd[2] / r1;
          npt.p[1] += V0 * sin(M_PI * r1 / g.LLn) * (crd[1] + .5 * g.LLy) / r1;
        }
      }
      if (r2 < g.LLn) {
        npt.n += (g.nn0 - g.nnb) * sqr(cos(M_PI / 2. * r2 / g.LLn));
        if (r2 > 0.0) {
          npt.p[2] += V0 * sin(M_PI * r2 / g.LLn) * crd[2] / r2;
          npt.p[1] += V0 * sin(M_PI * r2 / g.LLn) * (crd[1] - .5 * g.LLy) / r2;
        }
      }

      switch (kind) {
        case KIND_ELECTRON:
          // electron drift consistent with initial current
          if ((r1 <= g.LLn) && (r1 >= g.LLn - 2. * g.LLB)) {
            npt.p[0] = -g.BB * M_PI / (2. * g.LLB) *
                       cos(M_PI * (g.LLn - r1) / (2. * g.LLB)) / npt.n;
          }
          if ((r2 <= g.LLn) && (r2 >= g.LLn - 2. * g.LLB)) {
            npt.p[0] = -g.BB * M_PI / (2. * g.LLB) *
                       cos(M_PI * (g.LLn - r2) / (2. * g.LLB)) / npt.n;
          }

          npt.T[0] = g.TTe;
          npt.T[1] = g.TTe;
          npt.T[2] = g.TTe;
          break;
        case KIND_ION:
          npt.T[0] = g.TTi;
          npt.T[1] = g.TTi;
          npt.T[2] = g.TTi;
          break;
        default: assert(0);
      }
    });
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  setupFields(mflds, [&](int m, double crd[3]) {
    double z1 = crd[2];
    double y1 = crd[1] + .5 * g.LLy;
    double r1 = sqrt(sqr(z1) + sqr(y1));
    double z2 = crd[2];
    double y2 = crd[1] - .5 * g.LLy;
    double r2 = sqrt(sqr(z2) + sqr(y2));

    double rv = 0.;
    switch (m) {
      case HZ:
        if ((r1 < g.LLn) && (r1 > g.LLn - 2 * g.LLB)) {
          rv += -g.BB * sin(M_PI * (g.LLn - r1) / (2. * g.LLB)) * y1 / r1;
        }
        if ((r2 < g.LLn) && (r2 > g.LLn - 2 * g.LLB)) {
          rv += -g.BB * sin(M_PI * (g.LLn - r2) / (2. * g.LLB)) * y2 / r2;
        }
        return rv;

      case HY:
        if ((r1 < g.LLn) && (r1 > g.LLn - 2 * g.LLB)) {
          rv += g.BB * sin(M_PI * (g.LLn - r1) / (2. * g.LLB)) * z1 / r1;
        }
        if ((r2 < g.LLn) && (r2 > g.LLn - 2 * g.LLB)) {
          rv += g.BB * sin(M_PI * (g.LLn - r2) / (2. * g.LLB)) * z2 / r2;
        }
        return rv;

      case EX:
        if ((r1 < g.LLn) && (r1 > g.LLn - 2 * g.LLB)) {
          rv += g.MMach * sqrt(g.TTe / g.MMi) * g.BB *
                sin(M_PI * (g.LLn - r1) / (2. * g.LLB)) *
                sin(M_PI * r1 / g.LLn);
        }
        if ((r2 < g.LLn) && (r2 > g.LLn - 2 * g.LLB)) {
          rv += g.MMach * sqrt(g.TTe / g.MMi) * g.BB *
                sin(M_PI * (g.LLn - r2) / (2. * g.LLB)) *
                sin(M_PI * r2 / g.LLn);
        }
        return rv;

        // FIXME, JXI isn't really needed anymore (?)
      case JXI:
        if ((r1 < g.LLn) && (r1 > g.LLn - 2 * g.LLB)) {
          rv += g.BB * M_PI / (2. * g.LLB) *
                cos(M_PI * (g.LLn - r1) / (2. * g.LLB));
        }
        if ((r2 < g.LLn) && (r2 > g.LLn - 2 * g.LLB)) {
          rv += g.BB * M_PI / (2. * g.LLB) *
                cos(M_PI * (g.LLn - r2) / (2. * g.LLB));
        }
        return rv;

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

static void run()
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
  psc_params.balance_interval = 0;
  Balance balance{.1};

  // -- Sort
  psc_params.sort_interval = 10;

  // -- Collision
  int collision_interval = 10;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = 0; // 5
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
