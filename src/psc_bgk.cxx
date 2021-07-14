
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

#include <iomanip>

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

using Mfields = PscConfig::Mfields;
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

    assert(rho <= data[n_rows - 1][COL_RHO]);

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

  // ion mass
  double m_i;

  // electron charge
  double q_e;

  // electron mass
  double m_e;

  // number of grid cells
  int n_grid;

  // number of patches
  int n_patches;

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

  g.q_e = -1;

  g.m_i = 1e9;

  g.m_e = 1;

  g.n_grid = 128;

  g.n_patches = std::max(g.n_grid / 16, 16);

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
    {1, g.n_patches, g.n_patches}};           // # patches

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {g.q_e, g.m_e, "e"};
  kinds[KIND_ION] = {g.q_i, g.m_i, "i"}; // really heavy to keep them fixed

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

// ======================================================================
// getPadded

template <typename GT>
auto getPadded(const GT& gt, Int3 paddings)
{
  auto shape = gt.shape();
  for (int i = 0; i < 3; i++)
    shape[i] += paddings[i] * 2;

  auto padded =
    gt::full<gt::expr_value_type<GT>, gt::expr_space_type<GT>>(shape, 0);

  view_interior(padded, paddings) = gt;

  return padded;
}

// ======================================================================
// printGT

template <typename GT>
void printGT(const GT& gt)
{
  auto shape = gt.shape();
  for (int p = 0; p < shape[4]; p++) {
    std::cout << "===== " << p << " =====\n";
    for (int j = 0; j < shape[1]; ++j) {
      for (int k = 0; k < shape[2]; ++k) {
        std::cout << std::setprecision(0) << std::scientific
                  << gt(0, j, k, 1, p) << '\t';
      }
      std::cout << "\n";
    }
  }
}

// ======================================================================
// fillGhosts

int getNeighborPatch(int p, int dy, int dz)
{
  // FIXME THIS WILL NOT WORK EXCEPT IN YZ PLANE
  int iy = (p % g.n_patches + dy) % g.n_patches;
  int iz = (p / g.n_patches + dz) % g.n_patches;
  iy += iy < 0 ? g.n_patches : 0;
  iz += iz < 0 ? g.n_patches : 0;
  return iz * g.n_patches + iy;
}

template <typename GT>
void fillGhosts(GT&& gt, Int3 bnd)
{
  auto shape = gt.shape();

  auto rev = _s(_, _, -1);
  auto&& gt_rev = gt.view(rev, rev, rev);

  for (int a = 0; a < 3; ++a) {
    for (int p = 0; p < shape[4]; ++p) {
      gt::gslice sLhs[] = {_all, _all, _all};
      gt::gslice sRhs[] = {_all, _all, _all};

      // selects lower half of ghost points along axis a
      sLhs[a] = _s(_, bnd[a]);
      // selects source for those ghost points
      sRhs[a] = _s(-2 * bnd[a], -bnd[a]);

      int dy = a == 1 ? 1 : 0;
      int dz = a == 2 ? 1 : 0;

      gt.view(sLhs[0], sLhs[1], sLhs[2], _all, p) =
        gt.view(sRhs[0], sRhs[1], sRhs[2], _all, getNeighborPatch(p, dy, dz));
      gt_rev.view(sLhs[0], sLhs[1], sLhs[2], _all, p) = gt_rev.view(
        sRhs[0], sRhs[1], sRhs[2], _all, getNeighborPatch(p, -dy, -dz));
    }
  }
}

template <typename GT>
auto getPaddedWithGhosts(const GT& gt, Int3 bnd)
{
  auto padded = getPadded(gt, bnd);
  fillGhosts(padded, bnd);
  return padded;
}

// ======================================================================
// getCoord

inline double getCoord(double crd)
{
  if (crd < -g.box_size / 2)
    return crd + g.box_size;
  if (crd > g.box_size / 2)
    return crd - g.box_size;
  return crd;
}

// ======================================================================
// initializeParticles

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts,
                         Mfields& phi)
{
  SetupParticles<Mparticles> setup_particles(*grid_ptr);
  setup_particles.centerer = Centering::Centerer(Centering::NC);

  Int3 ibn{0, 2, 2}; // FIXME don't hardcode
  auto gradPhi = getPaddedWithGhosts(Item_grad<Mfields>(phi).gt(), ibn);
  auto divGradPhi = psc::item::div_nc(gradPhi, *grid_ptr);

  writeGT(divGradPhi, *grid_ptr, "divgrad", {"divgrad"});

  auto npt_init = [&](int kind, double crd[3], int p, Int3 idx,
                      psc_particle_npt& npt) {
    double y = getCoord(crd[1]);
    double z = getCoord(crd[2]);
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
        npt.n =
          (-divGradPhi(idx[0], idx[1], idx[2], 0, p) - g.n_i * g.q_i) / g.q_e;

        // temperature
        npt.T[0] = parsed.get_interpolated(COL_TE, rho);
        npt.T[1] = npt.T[2] = npt.T[0];
        break;
      case KIND_ION:
        // momenta are 0

        // number density
        npt.n = g.n_i;

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

void initializePhi(Mfields& phi)
{
  setupScalarField(
    phi, Centering::Centerer(Centering::NC), [&](int m, double crd[3]) {
      double rho = sqrt(sqr(getCoord(crd[1])) + sqr(getCoord(crd[2])));
      return parsed.get_interpolated(COL_PHI, rho);
    });
}

// ======================================================================
// initializeE from phi

void initializeE(MfieldsState& mflds, Mfields& phi)
{
  auto ibn = mflds.ibn();

  auto grad_item = Item_grad<Mfields>(phi);
  auto&& grad = getPadded(grad_item.gt(), ibn);
  fillGhosts(grad, ibn);

  // write results so they can be checked
  writeGT(view_interior(grad, ibn), grad_item.grid(), grad_item.name(),
          grad_item.comp_names());

  mflds.storage().view(_all, _all, _all, _s(EX, EX + 3)) = -grad;
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case HX: return g.Hx;
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
  Mfields phi{grid, 1, mflds.ibn()};

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
  writeGT(view_interior(phi.gt(), phi.ibn()), phi.grid(), "phi", {"phi"});

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
