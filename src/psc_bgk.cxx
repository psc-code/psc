
#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"

// for parsing
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

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

using PhiField = PscConfig::Mfields;
using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

// ======================================================================
// Parsed

enum
{
  COL_RHO,
  COL_NE,
  COL_V_PHI,
  COL_TE,
  COL_E_RHO,
  COL_PHI,
  n_cols
};

class Parsed
{
private:
  static const int n_rows = 4401;
  static constexpr auto file_path = "../../input/bgk-input-1-phi.txt";

  // the data parsed from file_path
  double data[n_rows][n_cols];
  // step size between first and second entries of rho (and presumably between
  // all consecutive entries)
  double rho_step;

  // gets index of row containing greatest lower bound of given rho
  int get_row(double rho)
  {
    // initial guess; should be precise assuming rho is linearly spaced
    int row = std::max(int(rho / rho_step), n_rows - 1);
    while (rho < data[row][COL_RHO])
      row--;
    return row;
  }

public:
  Parsed()
  {
    // first populate data
    // note the use of normalized version (no carriage returns allowed!!)
    std::ifstream file(file_path);

    // skip first line
    file.ignore(256, '\n');

    // iterate over each line
    int row = 0;
    std::string line;

    while (std::getline(file, line)) {
      assert(row < n_rows);

      // iterate over each entry within a line
      std::istringstream iss(line);
      std::string result;
      int col = 0;

      while (std::getline(iss, result, '\t')) {
        assert(col < n_cols);

        // write entry to data
        data[row][col] = stod(result);
        col++;
      }
      assert(col == n_cols);
      row++;
    }
    assert(row == n_rows);
    file.close();

    // initialize other values
    rho_step = data[1][COL_RHO] - data[0][COL_RHO];
  }

  // calculates and returns an interpolated value (specified by col) at
  // corresponding rho
  double get_interpolated(int col, double rho)
  {
    // ensure we are in bounds
    assert(rho >= data[0][COL_RHO]);
    // FIXME this assert fails because setupFields requests out-of-bounds points
    // assert(rho <= data[n_rows - 1][COL_RHO]);

    int row = get_row(rho);

    // check if we can return an exact value

    if (data[row][COL_RHO] == rho)
      return data[row][col];

    // otherwise, return linear interpolation

    // weights
    double w1 = data[row + 1][COL_RHO] - rho;
    double w2 = rho - data[row][COL_RHO];

    return (w1 * data[row][col] + w2 * data[row + 1][col]) / (w1 + w2);
  }
};

// ======================================================================
// PscBgkParams

struct PscBgkParams
{
  // physical length of region along y and z
  double box_size;

  // strength of transverse magnetic field
  double Hx;

  // ion charge
  double q_i;

  // ion number density
  double n_i;

  // number of grid cells
  int n_grid;

  // whether or not to negate v (if true, we should not see stability)
  bool reverse_v;
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

Parsed parsed;

// Parameters specific to this case. They don't really need to be collected in a
// struct, but maybe it's nice that they are
PscBgkParams g;

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
  psc_params.nmax = 5000; // 32000;
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
  g.box_size = .03;

  g.Hx = .1;

  g.n_i = 1;

  g.q_i = 1;
  // g.q_i = 1.0000111539638505; // from psc-scrap/check_case1.ipynb

  g.n_grid = 32;

  g.reverse_v = true;
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
    {1, g.n_grid, g.n_grid},                  // # grid points
    {1, g.box_size, g.box_size},              // physical lengths
    {0., -.5 * g.box_size, -.5 * g.box_size}, // *offset* for origin
    {1, g.n_grid / 16, g.n_grid / 16}};       // # patches

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1., 1., "e"};
  kinds[KIND_ION] = {1., 1e9, "i"}; // really heavy to keep them fixed

  mpi_printf(MPI_COMM_WORLD, "lambda_D = %g\n",
             sqrt(parsed.get_interpolated(COL_TE, .022)));

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
// write phi

void writePhi(PhiField& phi)
{
  static WriterMRC writer;
  if (!writer) {
    writer.open("phi");
  }
  auto& grid = phi.grid();
  writer.begin_step(grid.timestep(), grid.timestep() * grid.dt);
  writer.write(view_interior(phi.gt(), phi.ibn()), grid, "phi", {"phi"});
  writer.end_step();
}

// ======================================================================
// write grad

void writeGrad(Item_grad<PhiField>& grad)
{
  static WriterMRC writer;
  if (!writer) {
    writer.open("grad");
  }
  auto& grid = grad.grid();
  writer.begin_step(grid.timestep(), grid.timestep() * grid.dt);
  writer.write(grad.gt(), grad.grid(), grad.name(), grad.comp_names());
  writer.end_step();
}

// ======================================================================
// initializeParticles

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts,
                         PhiField& phi)
{
  SetupParticles<Mparticles> setup_particles(*grid_ptr);

  // auto gradPhi = Item_grad<PhiField>(phi).gt();
  // auto divGradPhi = psc::item::div_nc(gradPhi, *grid_ptr);

  auto npt_init = [&](int kind, double crd[3], int p, Int3 idx,
                      psc_particle_npt& npt) {
    double y = crd[1];
    double z = crd[2];
    double rho = sqrt(sqr(y) + sqr(z));

    switch (kind) {
      case KIND_ELECTRON:
        // velocities/momenta
        if (rho != 0) {
          double v_phi = parsed.get_interpolated(COL_V_PHI, rho);
          double sign = g.reverse_v ? -1 : 1;
          npt.p[1] = sign * v_phi * -z / rho;
          npt.p[2] = sign * v_phi * y / rho;
        }
        // otherwise, npt.p[i] = 0

        // number density
        npt.n = parsed.get_interpolated(COL_NE, rho);

        // temperature
        npt.T[0] = parsed.get_interpolated(COL_TE, rho);
        npt.T[1] = npt.T[2] = npt.T[0];
        break;
      case KIND_ION:
        // momenta are 0

        // number density
        npt.n = 1;

        // temperature is 0
        break;
      default: assert(false);
    }
  };

  partitionParticlesGeneralInit(setup_particles, balance, grid_ptr, mprts,
                                npt_init);
  setupParticlesGeneralInit(setup_particles, mprts, npt_init);
}

// ======================================================================
// initializePhi

void initializePhi(PhiField& phi)
{
  setupScalarField(phi, [&](int m, double crd[3]) {
    double rho = sqrt(sqr(crd[1]) + sqr(crd[2]));
    return parsed.get_interpolated(COL_PHI, rho);
  });
}

// ======================================================================
// initializeE from phi

void initializeE(MfieldsState& mflds, PhiField& phi)
{
  auto grad_item = Item_grad<PhiField>(phi);
  // write results so they can be checked
  writeGrad(grad_item);

  auto grad = grad_item.gt();

  view_interior(mflds.storage(), mflds.ibn())
    .view(_all, _all, _all, _s(EX, EX + 3)) = -grad;
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  setupFields(mflds, [&](int m, double crd[3]) {
    // take care of the easy cases
    switch (m) {
      // case EY:
      // case EZ: break;
      case HX: return g.Hx;
      default: return 0.;
    }

    // otherwise, interpolate from file

    double y = crd[1];
    double z = crd[2];
    double rho = sqrt(sqr(y) + sqr(z));

    if (rho == 0)
      return 0.;

    double E_rho = parsed.get_interpolated(COL_E_RHO, rho);

    if (m == EY) {
      return E_rho * y / rho;
    } else if (m == EZ) {
      return E_rho * z / rho;
    }

    return 0.0;
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
  PhiField phi{grid, 1, mflds.ibn()};

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  // -- Balance
  psc_params.balance_interval = 0;
  Balance balance{psc_params.balance_interval, .1};

  // -- Sort
  psc_params.sort_interval = 10;

  // -- Collision
  int collision_interval = 10;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.gauss_every_step = 200;
  // checks_params.gauss_dump_always = true;
  checks_params.gauss_threshold = 1e-4;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = 1; // 5
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output
  //
  // FIXME, this really is too complicated and not very flexible

  // -- output fields
  OutputFieldsParams outf_params{};
  outf_params.fields.pfield_interval = 200;
  outf_params.moments.pfield_interval = 200;
  OutputFields<MfieldsState, Mparticles, Dim> outf{grid, outf_params};

  // -- output particles
  OutputParticlesParams outp_params{};
  outp_params.every_step = 0;
  outp_params.data_dir = ".";
  outp_params.basename = "prt";
  OutputParticles outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // setup initial conditions

  if (read_checkpoint_filename.empty()) {
    initializePhi(phi);
    initializeParticles(balance, grid_ptr, mprts, phi);
    initializeFields(mflds);
    initializeE(mflds, phi);
  } else {
    read_checkpoint(read_checkpoint_filename, *grid_ptr, mprts, mflds);
  }

  // ----------------------------------------------------------------------
  // write phi (so it can be visually checked)
  writePhi(phi);

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
