

#include <psc.hxx>
#include <psc_bits.h>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"

#ifdef USE_CUDA
#include "cuda_bits.h"
#endif
#define PRINT(c)                                                               \
  {                                                                            \
    std::cout << "JOHN " << c << std::endl;                                    \
  }
// ======================================================================
// Particle kinds
//
// Particle kinds can be used to define different species, or different
// populations of the same species
//
// Here, we only enumerate the types, the actual information gets set up later.
// The last kind (MY_ELECTRON) will be used as "neutralizing kind", ie, in the
// initial setup, the code will add as many electrons as there are ions in a
// cell, at the same position, to ensure the initial plasma is neutral
// (so that Gauss's Law is satisfied).
enum
{
  MY_ELECTRON,
  MY_ION,
  MY_ELECTRON_B,
  // MY_ELECTRON_BG,
  // MY_ION_BG,
  N_MY_KINDS,
};

// ======================================================================
// PscHarrisParams

struct PscHarrisParams
{
  double Lx_di, Ly_di, Lz_di; // Size of box in d_i

  double L_di; // Sheet thickness / ion inertial length

  double BB;
  double Zi;
  double mass_ratio;
  double Ti_Te;
  double Tib_Ti, Teb_Te;
  double nb_n0; // background density
  double ni_n0;

  double bg; // guide field as fraction of B0
  double theta;
  double dby_b0;
  double Lpert_Lz;

  // The following parameters are calculated from the above / and other
  // information

  double d_i;
  double b0; // B0
  // double n0;
  double L; // sheet width in d_e
  // double TTi, TTe;
  double wpe_wce;
  double vae_c;
  double vap_c;
  double v_scale;

  double beta_e_par;
  double beta_i_par;
  double Ti_perp_over_Ti_par;
  double Te_perp_over_Te_par;

  double B0;

  double me;
  double c;
  double eps0;
  double ec;
  double de;

  double amplitude;
  double Te_par;
  double Te_perp;
  double Ti_par;
  double Ti_perp;

  double wpe, wpi, wce, wci;
  double Lx, Ly, Lz; // size of box
  double nx, ny, nz; // number of boxes
  double Lpert;      // wavelength of perturbation
  double dby;        // Perturbation in Bz relative to Bo (Only change here)
  double dbz;        // Set Bx perturbation so that div(B) = 0
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

PscHarrisParams g;

std::string read_checkpoint_filename;

// This is a set of generic PSC params (see include/psc.hxx),
// like number of steps to run, etc, which also should be set by the case
PscParams psc_params;

} // namespace
// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_yz;

#if 0

#ifdef USE_CUDA
using PscConfig = PscConfig1vbecCuda<Dim>;
#else
using PscConfig = PscConfig1vbecSingle<Dim>;
#endif

#else

#include "particle_with_id.h"

using PscConfig =
  PscConfig_<Dim, MparticlesSimple<ParticleWithId<float>>, MfieldsStateSingle,
             MfieldsSingle, PscConfigPushParticles1vbec>;
#endif

using Writer = WriterMRC; // can choose WriterMrc, WriterAdios2

// ----------------------------------------------------------------------

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

// ======================================================================
// setupParameters

void setupParameters()
{
  // -- set some generic PSC parameters
  psc_params.nmax = 200001; // 5001;
  psc_params.cfl = 0.05;
  psc_params.write_checkpoint_every_step = -100;
  psc_params.stats_every = 1;

  // -- start from checkpoint:
  //
  // Uncomment when wanting to start from a checkpoint, ie.,
  // instead of setting up grid, particles and state fields here,
  // they'll be read from a file
  // FIXME: This parameter would be a good candidate to be provided
  // on the command line, rather than requiring recompilation when change.

  // read_checkpoint_filename = "checkpoint_500.bp";

  // -- Set some parameters specific to this case

  // PIC units, only used in this scope
  g.me = 1;   // Mass normalization
  g.ec = 1;   // Charge normalization
  g.c = 1;    // Speed of light
  g.eps0 = 1; // Permittivity of space
  g.de = 1;   // Length normalization (electron inertial length)

  g.Lx_di = 1.;
  g.Ly_di = 5.;
  g.Lz_di = 10.;
  g.L_di = 1.;
  g.Lpert_Lz = 1.;
  g.amplitude = .5;
  g.beta_e_par = 0.005;
  g.beta_i_par = 0.005;

  g.Zi = 1.;
  g.mass_ratio = 500.;
  g.Ti_Te = 5.;
  g.Tib_Ti = 0.333;
  g.Teb_Te = 0.333;
  g.Ti_perp_over_Ti_par = 1.;
  g.Te_perp_over_Te_par = 1.;

  g.vap_c = 0.01;
  g.vae_c = sqrt(g.mass_ratio) * g.vap_c;
  // g.v_scale = 0.1;

  g.nb_n0 = 0.05;
  g.ni_n0 = 0.05;
  g.theta = 0;

  g.B0 = sqrt(g.mass_ratio) * g.vap_c;

  g.wpe_wce = 2.;
  g.Te_par = g.beta_e_par * sqr(g.B0) / 2.;
  g.Te_perp = g.Te_perp_over_Te_par * g.Te_par;
  g.Ti_par = g.beta_i_par * sqr(g.B0) / 2.;
  g.Ti_perp = g.Ti_perp_over_Ti_par * g.Ti_par;
  g.Te_par = 0.5 * g.beta_e_par * g.vap_c * g.vap_c * g.mass_ratio;
  // g.Te_par = 10.;
  g.Ti_par = 0.5 * g.mass_ratio * g.beta_e_par * g.vap_c * g.vap_c;

  g.wci = 1. / (g.mass_ratio * g.wpe_wce); // Ion cyclotron frequency
  g.wce = g.wci * g.mass_ratio;            // Electron cyclotron freqeuncy
  g.wpe = g.wce * g.wpe_wce;               // electron plasma frequency
  g.wpi = g.wpe / sqrt(g.mass_ratio);      // ion plasma frequency

  g.d_i = g.c / g.wpi; // ion inertial length
  g.L = g.L_di * g.d_i;
  g.Lx = g.Lx_di * g.d_i; // size of box in x dimension
  g.Ly = g.Ly_di * g.d_i; // size of box in y dimension
  g.Lz = g.Lz_di * g.d_i; // size of box in z dimension
  g.nx = 1.;              // number of box in x dimension
  g.ny = 640.;            // number of box in y dimension
  g.nz = 960.;            // number of box in z dimension
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
  Grid_t::Kinds kinds(N_MY_KINDS);
  kinds[MY_ION] = {g.Zi, g.mass_ratio, "i"};
  kinds[MY_ELECTRON] = {-1., 1., "e"};
  kinds[MY_ELECTRON_B] = {-1., 1., "eb"};
  // kinds[MY_ION_BG] = {g.Zi, g.mass_ratio * g.Zi, "i_bg"};
  // kinds[MY_ELECTRON_BG] = {-1., 1., "e_bg"};

  mpi_printf(MPI_COMM_WORLD, "d_e = %g, d_i = %g\n", 1., g.d_i);

  Grid_t::Real3 LL = {1., 5., 10.};
  Int3 gdims = {1, 640, 960};
  Int3 np = {1, 40, 60};

  auto domain = Grid_t::Domain{gdims, LL, -.5 * LL, np};

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 100;

  // mprintf("dx %g %g %g\n", domain.dx[0], domain.dx[1], domain.dx[2]);
  // mprintf("1/dx^2
  // %g\n", 1./(sqr(domain.dx[0])+sqr(domain.dx[1])+sqr(domain.dx[2])));
  double dt = psc_params.cfl * courant_length(domain);
  // mprintf("dt %g cfl %g\n", dt, psc_params.cfl);
  // mprintf("courant_length %g\n", courant_length(domain));

  Grid_t::Normalization norm{norm_params};

  // mpi_printf(MPI_COMM_WORLD, "dt = %g\n", dt);

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

void initializeParticles(SetupParticles<Mparticles>& setup_particles,
                         Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts)
{

  // -- set particle initial condition
  partitionAndSetupParticles(
    setup_particles, balance, grid_ptr, mprts,
    [&](int kind, Double3 crd, psc_particle_npt& npt) {
      double x = crd[0], y = crd[1], z = crd[2];
      double a = g.B0 / 3, b = (2 * M_PI / grid_ptr->domain.length[1]),
             c = 2 * (2 * M_PI / grid_ptr->domain.length[2]);
      double constant = 0.1;
      double B0_y = a * sin(b * y) * cos(c * z);
      double B0_z = a * b * cos(b * y) * sin(c * z) / c + g.B0;
      double B_mag = sqrt(B0_y * B0_y + B0_z * B0_z);
      double u = 0.05;
      double uby = u * B0_y / B_mag;
      double ubz = u * B0_z / B_mag;
      double nx = 0.0001;
      double ny = -a * a / 4 * cos(c * z) * cos(c * z) * cos(2 * b * y) *
                  (1 + b * b / (c * c));
      double nz = a * a / 4 * sin(b * y) * sin(b * y) * cos(2 * c * z) *
                  (1 + b * b / (c * c));
      npt.tag = psc::particle::Tag{kind * 10};
      switch (kind) {
        case 0: // electron
          npt.p[0] = 0.;
          npt.p[1] = -uby / 9.;
          npt.p[2] = -ubz / 9.;
          npt.T[0] = g.Te_par;
          npt.T[1] = g.Te_par;
          npt.T[2] = g.Te_par;
          npt.kind = MY_ELECTRON;
          // mprintf("b %g\n", b);
          // mprintf("c %g\n", c);
          // mprintf("ly %g\n", grid_ptr->domain.length[1]);
          // mprintf("lz %g\n", grid_ptr->domain.length[2]);
          // mprintf("ne %g\n", npt.n);
          // mprintf("ne ny %g\n", ny);
          // mprintf("By %g\n", B0_y);
          // mprintf("B general %g\n", (B0_y * B0_y + B0_z * B0_z + g.B0*g.B0));
          // mprintf("T[0]e %g\n", npt.T[0]);
          // mprintf("T[1]e %g\n", npt.T[1]);
          // mprintf("T[2]e %g\n", npt.T[2]);
          // mprintf("T general p %g\n", npt.T);

          break;
        case 1: // ion
          npt.n = 1.;
          npt.T[0] = g.Te_par;
          npt.T[1] = g.Te_par;
          npt.T[2] = g.Te_par;
          npt.p[0] = 0.;
          npt.p[1] = 0.;
          npt.p[2] = 0.;
          npt.kind = MY_ION;
          // mprintf("np %g\n", npt.n);
          // mprintf("np ny %g\n", ny);
          // mprintf("B general %g\n", (B0_y * B0_y + B0_z * B0_z + g.B0*g.B0));
          // mprintf("T[0]p %g\n", npt.T[0]);
          // mprintf("T[1]p %g\n", npt.T[1]);
          // mprintf("T[1]p %g\n", npt.T[2]);
          // mprintf("T general e %g\n", npt.T);
          break;
        case 2: // electron beam
          npt.n = 0.1;

          npt.p[0] = 0.;
          npt.p[1] = uby;
          npt.p[2] = ubz;
          npt.T[0] = g.Te_par;
          npt.T[1] = g.Te_par;
          npt.T[2] = g.Te_par;

          npt.kind = MY_ELECTRON_B;

          break;
        // case 3: // ion bg
        //   npt.n = g.ni_n0;
        //   npt.p[0] = 0.;
        //   npt.T[0] = g.Tib_Ti * g.Ti_perp;
        //   npt.T[1] = g.Tib_Ti * g.Ti_perp;
        //   npt.T[2] = g.Tib_Ti * g.Ti_par;
        //   npt.kind = MY_ION;
        //   break;
        default: assert(0);
      }
    });
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  const auto& grid = mflds.grid();

  // mprintf("L %g\n", L);
  setupFields(mflds, [&](int m, double crd[3]) {
    double x = crd[0], y = crd[1], z = crd[2];
    double a = g.B0 / 3, b = (2 * M_PI / g.Ly_di), c = 2 * (2 * M_PI / g.Lz_di);
    // double a = 3., b = 1.5., c = 3.;
    double B0_y = a * sin(b * y) * cos(c * z);
    double B0_z = a * b * cos(b * y) * sin(c * z) / c + g.B0;

    switch (m) {
        // case HX: return 0.;

        // case HY: return 0.;

      case HX: return 0.;

      case HY: return B0_y;

      case HZ:
        return B0_z;

        // case JXI: return 0.; // FIXME

      default: return 0.;
    }
  });
}

void run()
{
  mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

  // ----------------------------------------------------------------------
  // setup various parameters first

  setupParameters();

  // ----------------------------------------------------------------------
  // Set up grid, state fields, particles

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Mparticles mprts(grid);
  MfieldsState mflds(grid);
  if (!read_checkpoint_filename.empty()) {
    read_checkpoint(read_checkpoint_filename, grid, mprts, mflds);
  }

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  // -- Balance
  psc_params.balance_interval = 2500;
  Balance balance{3};

  // -- Sort
  psc_params.sort_interval = 10;

  // -- Collision
  int collision_interval = 0;
  double collision_nu = 1e-10;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.continuity_every_step = 0;
  checks_params.continuity_dump_always = false;
  checks_params.continuity_threshold = 1e-4;
  checks_params.continuity_verbose = true;

  checks_params.gauss_every_step = -100;
  checks_params.gauss_dump_always = false;
  checks_params.gauss_threshold = 1e-4;
  checks_params.gauss_verbose = true;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = 100;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output
  //
  // FIXME, this really is too complicated and not very flexible

  // -- output fields
  OutputFieldsItemParams outf_item_params{};
  OutputFieldsParams outf_params{};
  outf_item_params.pfield.out_interval = 500; // produce steps
  outf_item_params.tfield.out_interval = -4;
  outf_item_params.tfield.average_every = 50;

  outf_params.fields = outf_item_params;
  outf_params.moments = outf_item_params;
  OutputFields<MfieldsState, Mparticles, Dim, Writer> outf{grid, outf_params};

  // -- output particles
  OutputParticlesParams outp_params{};
  outp_params.every_step = 2000;
  outp_params.data_dir = ".";
  outp_params.basename = "prt";
  OutputParticles outp{grid, outp_params};
  outp_params.lo = {1, 336, 256};
  outp_params.hi = {1, 432, 384};

  int oute_interval = -10;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // Set up objects specific to the Harris case

  SetupParticles<Mparticles> setup_particles(grid);
  setup_particles.fractional_n_particles_per_cell = true;
  setup_particles.neutralizing_population = 4 - 1; // FIXME?!

  // ----------------------------------------------------------------------
  // setup initial conditions

  if (read_checkpoint_filename.empty()) {
    initializeParticles(setup_particles, balance, grid_ptr, mprts);
    initializeFields(mflds);
  }

  // ----------------------------------------------------------------------
  // hand off to PscIntegrator to run the simulation

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

  MEM_STATS();
  psc.integrate();
  MEM_STATS();
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  psc_init(argc, argv);

  run();

  MEM_STATS();

  psc_finalize();
  return 0;
}
