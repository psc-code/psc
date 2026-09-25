/**
 * This script makes a small amount of data to be used as sample data for
 * pscpy (https://github.com/psc-code/pscpy). As psc output formats change,
 * this script is to be rerun to regenerate pscpy's samples.
 */

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "output_fields.hxx"
#include "psc_config.hxx"
#include "libpsc/psc_output_particles/output_particles_adios2_impl.hxx"

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

// ----------------------------------------------------------------------

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;

// ======================================================================
// Global parameters

PscParams psc_params;

// ======================================================================
// setupGrid

Grid_t* setupGrid()
{
  Int3 gdims = {1, 8, 4};
  Double3 lengths = {1.0, 10.0, 5.0};
  Double3 corner = {0.0, -5.0, -2.5};
  Int3 n_patches = {1, 4, 2};

  auto domain = Grid_t::Domain{gdims, lengths, corner, n_patches};

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 100.0, "i"};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 10;

  double dt = 0.9 * courant_length(domain);
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

  auto init_np = [&](int kind, Double3 pos, int p, Int3 idx,
                     psc_particle_np& np) {
    double t = 1e-3;
    Double3 u = {0.0, 1e-3, 0.0};
    np.n = 1.0;
    np.p =
      setup_particles.createMaxwellian({np.kind, np.n, u, {t, t, t}, np.tag});
  };

  partitionAndSetupParticles(setup_particles, balance, grid_ptr, mprts,
                             init_np);
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case HX: return 0.1;
      default: return 0.;
    }
  });
}

// ======================================================================
// run

static void run(int argc, char** argv)
{
  mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

  // ----------------------------------------------------------------------
  // Set up grid, state fields, particles

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  psc_params.nmax = 2;

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
  checks_params.gauss.check_interval = 1;
  checks_params.gauss.dump_always = true;
  checks_params.continuity.check_interval = 1;
  checks_params.continuity.dump_always = true;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // ----------------------------------------------------------------------
  // Set up output

  // -- output fields
  OutputFields<MfieldsState, Mparticles> out_fields;
  out_fields.pfield.out_interval = 1;

  OutputMoments<MfieldsState, Mparticles, Dim> out_moments;
  out_moments.pfield.out_interval = 1;

  // -- output particles
  OutputParticlesAdios2Params outp_params{};
  outp_params.every_step = 1;
  OutputParticlesAdios2<Mparticles, float> outp{grid, outp_params};

  // ----------------------------------------------------------------------
  // set up initial conditions

  initializeParticles(balance, grid_ptr, mprts);

  // ----------------------------------------------------------------------
  // run the simulation

  auto psc = makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts,
                                          balance, collision, checks);

  psc.add_diagnostic(&out_fields);
  psc.add_diagnostic(&out_moments);
  psc.add_diagnostic(&outp);

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
