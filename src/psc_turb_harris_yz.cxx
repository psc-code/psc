
#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"
#include "writer_mrc.hxx"

#ifdef USE_CUDA
#include "cuda_bits.h"
#endif

#include "rngpool_iface.h"

//-----------------------------------------------------------------
// To use the complex numbers
//------------------------------
#include <complex>
// using std::complex;
using namespace std;
typedef complex<double> dcomp;

// using namespace std::complex_literals;
//    auto c = 1.0 + 3.0i;
//------------------------------

// extern Grid* vgrid; // FIXME

// To use the random numbers
//-------------------------------
static RngPool* rngpool;
//-------------------------------

// static inline double trunc_granular(double a, double b)
//{
//  return b * (int)(a / b);
//}
//-----------------------------------------------------------------

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
  MY_ELECTRON_UP,
  MY_ION_UP,
  // MY_ELECTRON_BO,
  // MY_ION_BO,
  N_MY_KINDS,
};

enum
{
  PERT_HX,
  PERT_HY,
  PERT_HZ,
  PERT_AZ,
  DIV_B,
  PERT_JX_ext,
  PERT_JY_ext,
  PERT_JZ_ext,
  N_PERT,
};

// ======================================================================
// PscTurbHarrisxzParams

struct PscTurbHarrisxzParams
{

  //----------------------------------
  double BB;
  double Zi;
  double lambda0;
  double background_n;
  double background_Te;
  double background_Ti;
  //----------------------------------

  //----------------------------------
  double Lx_di, Ly_di, Lz_di; // Size of box in di
  double L_di;                // Sheet thickness / ion inertial length
  double Lpert_Lx;            // wavelength of perturbation in terms of Lx

  double taui;                     // simulation wci's to run
  double t_intervali;              // output interval in terms of 1/wci
  double output_particle_interval; // particle output interval in terms of 1/wci
  int ion_sort_interval;
  int electron_sort_interval;
  double overalloc; // Overallocation factor (> 1) for particle arrays

  double wpe_wce;   // electron plasma freq / electron cyclotron freq
  double mi_me;     // Ion mass / electron mass
  double wpedt_max; // Is rhis necessary?
  double Ti_Te;     // Ion temperature / electron temperature
  double nb_n0;     // background plasma density

  double Tbe_Te;  // Ratio of background T_e to Harris T_e
  double Tbi_Ti;  // Ratio of background T_i to Harris T_i
  double Tbi_Tbe; // Ratio of background T_i to background T_e

  double bg; // Guide field
  double theta;
  double cs;
  double sn;

  double db_b0;   // perturbation in Bz relative to B0
  double nppc;    // Average number of macro particle per cell per species
  bool open_bc_x; // Flag to signal we want to do open boundary condition in x
  bool
    driven_bc_z; // Flag to signal we want to do driven boundary condition in z

  Int3 gdims; //
  Int3 np;
  //----------------------------------

  //-------------------------------------
  double ec;
  double me;
  double c;
  double eps0;
  double de;
  double kb;

  double mi;
  double di;
  double wpe;
  double wpi;
  double wce;
  double wci;

  // calculated
  double b0; // B0
  double n0;
  double v_A;
  double rhoi_L;
  double Lx, Ly, Lz; // size of box
  double L;          // Harris sheet thickness
  double Lpert;      // wavelength of perturbation
  // double dbx;        // Perturbation in Bz relative to Bo (Only change here)
  // double dbz;        // Set Bx perturbation so that div(B) = 0
  double dbz; // Perturbation in Bz relative to Bo (Only change here) in the
              // case of yz
  double dby; // Set By perturbation so that div(B) = 0 in the case of yz
  double tanhf;

  double Te;  // Electron temperature main harris
  double Ti;  // Ion temperature main harris
  double Tbe; // Electron temperature backgroung
  double Tbi; // Ion temperature backgroung

  double weight_s; // Charge per macro electron

  double vthe; // Electron thermal velocity
  double vthi; // Ion thermal velocity
  double vdre; // Electron drift velocity
  double vdri; // Ion drift velocity

  double gdri; // gamma of ion drift frame
  double gdre; // gamma of electron drift frame
  double udri; // 4-velocity of ion drift frame
  double udre; // 4-velocity of electron drift frame

  double Ne_back;  // Number of macro electrons in background
  double weight_b; // Charge per macro electron
  double vtheb;    // normalized background e thermal vel.
  double vthib;    // normalized background ion thermal vel.

  int n_global_patches;

  double Ne;        // Total number of macro electrons
  double Ne_sheet;  // Number of macro electrons in Harris sheet
  double Npe_sheet; // N physical e's in sheet
  double Npe_back;  // N physical e's in backgrnd
  double Npe;       //??

  //--------------------------------------------------
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

PscTurbHarrisxzParams g;

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

// There is something not quite right here ask Kai Jeff
//--------------------------------------------------------------------------------
// using Dim = dim_yz;
using Dim = dim_xyz;
//--------------------------------------------------------------------------------

#ifdef USE_CUDA
using PscConfig = PscConfig1vbecCuda<Dim>;
#else
using PscConfig = PscConfig1vbecSingle<Dim>;
#endif

using Writer = WriterMRC; // can choose WriterMRC, WriterAdios2

// ----------------------------------------------------------------------

using MfieldsState = PscConfig::MfieldsState;
using MfieldsAlfven = Mfields<MfieldsState::real_t>;
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
  //-----------------------------------------------
  //-----------------------------------------------
  psc_params.nmax = 2; // 1801;
  psc_params.cfl = 0.75;
  psc_params.write_checkpoint_every_step = -100; // This is not working
  psc_params.stats_every = -1;
  //-----------------------------------------------
  //-----------------------------------------------

  // -- start from checkpoint:
  //
  // Uncomment when wanting to start from a checkpoint, ie.,
  // instead of setting up grid, particles and state fields here,
  // they'll be read from a file
  // FIXME: This parameter would be a good candidate to be provided
  // on the command line, rather than requiring recompilation when change.

  // read_checkpoint_filename = "checkpoint_500.bp";

  //-----------------------------------------------
  // -- Set some parameters specific to this case
  //----------------------------------
  g.BB = 1.;
  g.Zi = 1.;
  g.lambda0 = 20.;

  g.background_n = 1.;
  g.background_Te = .01;
  g.background_Ti = .01;
  //----------------------------------

  //----------------------------------
  // Space dimensions
  //--------------------------------------------------------------------------------
  ///*** // This is ion the case of yz geometry
  // g.Lz_di = 40.;
  // g.Lx_di = 1.;
  // g.Ly_di = 10.;
  // g.gdims = {1, 64, 256};
  g.Lz_di = 10.;
  g.Lx_di = 10.;
  g.Ly_di = 10.;
  g.gdims = {32, 32, 32};
  g.np = {2, 2, 2};
  //***/
  //--------------------------------------------------------------------------------

  g.L_di = 1.0;
  g.Lpert_Lx = 1.;

  // Time dimensions
  g.taui = 10.;
  g.t_intervali = 1.;
  g.output_particle_interval = 100.;
  g.electron_sort_interval = 25;
  g.ion_sort_interval = 25;
  g.overalloc = 2.;

  // Non-dimensional ratios
  g.wpe_wce = 2.5;
  g.mi_me = 1.;
  g.Ti_Te = 1.;
  g.nb_n0 = 0.1;

  g.Tbe_Te = .333; // How to stimate this ratio???? It is now consistent but I
                   // need an extra condition to estimate this ratio.
  g.Tbi_Ti = .333;
  g.Tbi_Tbe = 1.;

  // Background field
  g.bg = 0.;
  g.theta = 0.;
  g.cs = cos(g.theta);
  g.sn = sin(g.theta);

  // Amplitud of the fluctuation
  g.db_b0 = 0.1;

  // Number of macro particles
  g.nppc = 20;

  g.wpedt_max = .36; // what is this for?

  //----------------------------------

  //-----------------------------------
  // use natural PIC units
  g.ec = 1.;   // Charge normalization
  g.me = 1.;   // Mass normalization
  g.c = 1.;    // Speed of light
  g.de = 1.;   // Length normalization (electron inertial length)
  g.eps0 = 1.; // Permittivity of space
  g.wpe = 1.;  // g.wce * g.wpe_wce;             // electron plasma frequency
  g.kb = 1.;   //  k Boltzman

  // derived quantities
  g.wce = g.wpe / g.wpe_wce;     // Electron cyclotron frequency
  g.wci = g.wce / g.mi_me;       // Ion cyclotron frequency
  g.mi = g.me * g.mi_me;         // Ion mass
  g.wpi = g.wpe / sqrt(g.mi_me); // ion plasma frequency

  g.di = g.c / g.wpi; // ion inertial length
  g.L = g.L_di; // Harris sheet thickness in di // Jeff check the thickness. It
                // works best for g.L_di alone
  g.Lx = g.Lx_di; // size of box in x dimension (in di Jeff)
  g.Ly = g.Ly_di; // size of box in y dimension
  g.Lz = g.Lz_di; // size of box in z dimension

  g.b0 = g.me * g.c * g.wce / g.ec; // Asymptotic magnetic field strength
  g.n0 = g.me * g.eps0 * sqr(g.wpe) /
         (sqr(g.ec)); // Peak electron (ion) density this is the cgs correct one
                      // but gives n0 = 0.07

  // g.n0 = 1.;
  // g.b0 = 1.;

  g.Te = g.me * sqr(g.c) /
         (2. * g.kb * sqr(g.wpe_wce) *
          (1. + g.Ti_Te)); // Electron temperature neglecting the background
                           // plasma pressure

  // g.Te = sqr(g.b0) /
  //    (8. * M_PI * g.kb * g.n0 * (1. + g.Ti_Te));    // Electron temperature
  //    neglecting the background  plasma pressure

  // g.Te = sqr(g.b0) /
  //   (8. * M_PI * g.kb * g.n0 * ((1. + g.Ti_Te) - g.nb_n0 * g.Tbe_Te * (1 +
  //   g.Tbi_Tbe) ) );    // Electron temperature INCLUDING the background
  // plasma pressure which is important here due to the way psc inject particles

  g.Ti = g.Te * g.Ti_Te; // Ion temperature

  g.Tbe = g.Te * g.Tbe_Te;
  g.Tbi = g.Ti * g.Tbi_Ti;

  g.v_A = g.c * (g.wci /
                 g.wpi); //                      / sqrt(g.nb_n0); // based on nb
  // Include the relativistic alfven speed correction.
  g.rhoi_L = sqrt(g.Ti_Te / (1. + g.Ti_Te)) / g.L_di;

  g.vthe = sqrt(g.Te / g.me); // Electron thermal velocity
  g.vthi = sqrt(g.Ti / g.mi); // Ion thermal velocity
  g.vtheb =
    sqrt(g.Tbe_Te * g.Te / g.me); // normalized background e thermal vel.
  g.vthib =
    sqrt(g.Tbi_Ti * g.Ti / g.mi); // normalized background ion thermal vel.

  g.vdri =
    g.c * g.b0 /
    (8 * M_PI * g.L * g.ec * g.n0 * (1 + 1 / g.Ti_Te)); // Ion drift velocity
  g.vdre = -g.vdri / (g.Ti_Te); // electron drift velocity

  g.n_global_patches = g.np[0] * g.np[1] * g.np[2];

  g.Npe_sheet = 2 * g.n0 * g.Lx * g.Ly * g.L *
                tanh(0.5 * g.Lz / g.L); // N physical e's in sheet
  g.Npe_back =
    g.nb_n0 * g.n0 * g.Ly * g.Lz * g.Lx; // N physical e's in backgrnd
  g.Npe = g.Npe_sheet + g.Npe_back;

  g.Ne = g.nppc * g.gdims[0] * g.gdims[1] *
         g.gdims[2]; // total macro electrons in box
  g.Ne_sheet = g.Ne * g.Npe_sheet / g.Npe;
  g.Ne_back = g.Ne * g.Npe_back / g.Npe;

  // g.Ne_sheet = trunc_granular(g.Ne_sheet, g.n_global_patches); // Make it
  // divisible by # subdomains g.Ne_back = trunc_granular( g.Ne_back,
  // g.n_global_patches); // Make it divisible by # subdomains
  g.Ne = g.Ne_sheet + g.Ne_back;
  // g.weight_s = g.ec * g.Npe_sheet / g.Ne_sheet; // Charge per macro electron
  // g.weight_b = g.ec * g.Npe_back / g.Ne_back;   // Charge per macro electron

  g.gdri = 1. / sqrt(1. - sqr(g.vdri) / sqr(g.c)); // gamma of ion drift frame
  g.gdre =
    1. / sqrt(1. - sqr(g.vdre) / sqr(g.c)); // gamma of electron drift frame

  g.udri = g.vdri * g.gdri; // 4-velocity of ion drift frame
  g.udre = g.vdre * g.gdre; // 4-velocity of electron drift frame
  g.tanhf = tanh(0.5 * g.Lz / g.L);
  g.Lpert = g.Lpert_Lx * g.Lx; // wavelength of perturbation

  //--------------------------------------------------------------------------------
  // This is in the case of yz geometry
  g.dbz =
    g.db_b0 * g.b0 * M_PI * g.L /
    g.Ly; // Set Bz perturbation so that div(B) = 0 in the case of yz geometry
  g.dby = -g.db_b0 * g.b0 * 2 * M_PI * g.L /
          g.Lz; // Perturbation in By relative to Bo (Only change here) in the
                // case of yz geometry
  //--------------------------------------------------------------------------------

  //-----------------------------------
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

  //---------------------------------------------------------------
  mpi_printf(MPI_COMM_WORLD,
             "***********************************************\n");
  mpi_printf(MPI_COMM_WORLD, "* Topology: %d x %d x %d\n", g.np[0], g.np[1],
             g.np[2]);
  mpi_printf(MPI_COMM_WORLD, "tanhf    = %g\n", g.tanhf);
  mpi_printf(MPI_COMM_WORLD, "L_di     = %g\n", g.L_di);
  mpi_printf(MPI_COMM_WORLD, "rhoi/L   = %g\n", g.rhoi_L);
  mpi_printf(MPI_COMM_WORLD, "Ti/Te    = %g\n", g.Ti_Te);
  mpi_printf(MPI_COMM_WORLD, "Ti    = %g\n", g.Ti);
  mpi_printf(MPI_COMM_WORLD, "Te    = %g\n", g.Te);
  mpi_printf(MPI_COMM_WORLD, "nb/n0    = %g\n", g.nb_n0);
  mpi_printf(MPI_COMM_WORLD, "wpe/wce  = %g\n", g.wpe_wce);
  mpi_printf(MPI_COMM_WORLD, "mi/me    = %g\n", g.mi_me);
  mpi_printf(MPI_COMM_WORLD, "theta    = %g\n", g.theta);
  mpi_printf(MPI_COMM_WORLD, "Lpert/Lx = %g\n", g.Lpert_Lx);
  mpi_printf(MPI_COMM_WORLD, "dbz/b0   = %g\n", g.db_b0);
  mpi_printf(MPI_COMM_WORLD, "taui     = %g\n", g.taui);
  mpi_printf(MPI_COMM_WORLD, "t_intervali = %g\n", g.t_intervali);
  mpi_printf(MPI_COMM_WORLD, "num_step = %d\n", psc_params.nmax);
  mpi_printf(MPI_COMM_WORLD, "n0 = %g\n", g.n0);
  mpi_printf(MPI_COMM_WORLD, "Lx/di = %g\n", g.Lx);
  mpi_printf(MPI_COMM_WORLD, "Lx/de = %g\n", g.Lx * g.de / g.di);
  mpi_printf(MPI_COMM_WORLD, "Ly/di = %g\n", g.Ly);
  mpi_printf(MPI_COMM_WORLD, "Ly/de = %g\n", g.Ly * g.de / g.di);
  mpi_printf(MPI_COMM_WORLD, "Lz/di = %g\n", g.Lz);
  mpi_printf(MPI_COMM_WORLD, "Lz/de = %g\n", g.Lz * g.de / g.di);
  mpi_printf(MPI_COMM_WORLD, "nx = %d\n", g.gdims[0]);
  mpi_printf(MPI_COMM_WORLD, "ny = %d\n", g.gdims[1]);
  mpi_printf(MPI_COMM_WORLD, "nz = %d\n", g.gdims[2]);
  mpi_printf(MPI_COMM_WORLD, "n_global_patches = %d\n", g.n_global_patches);
  mpi_printf(MPI_COMM_WORLD, "nppc = %g\n", g.nppc);
  mpi_printf(MPI_COMM_WORLD, "b0 = %g\n", g.b0);
  mpi_printf(MPI_COMM_WORLD, "v_A (based on nb) = %g\n", g.v_A);
  mpi_printf(MPI_COMM_WORLD, "di = %g\n", g.di);
  mpi_printf(MPI_COMM_WORLD, "Ne = %g\n", g.Ne);
  mpi_printf(MPI_COMM_WORLD, "Ne_sheet = %g\n", g.Ne_sheet);
  mpi_printf(MPI_COMM_WORLD, "Ne_back = %g\n", g.Ne_back);
  mpi_printf(MPI_COMM_WORLD, "total # of particles = %g\n", 2 * g.Ne);
  // mpi_printf(MPI_COMM_WORLD, "dt*wpe = %g\n", g.wpe * grid.dt);
  // mpi_printf(MPI_COMM_WORLD, "dt*wce = %g\n", g.wce * grid.dt);
  // mpi_printf(MPI_COMM_WORLD, "dt*wci = %g\n", g.wci * grid.dt);
  mpi_printf(MPI_COMM_WORLD, "dx/de = %g\n", g.Lx / (g.de * g.gdims[0]));
  mpi_printf(MPI_COMM_WORLD, "dy/de = %g\n", g.Ly / (g.de * g.gdims[1]));
  mpi_printf(MPI_COMM_WORLD, "dz/de = %g\n", g.Lz / (g.de * g.gdims[2]));
  mpi_printf(MPI_COMM_WORLD, "dx/rhoi = %g\n",
             (g.Lx / g.gdims[0]) / (g.vthi / g.wci));
  mpi_printf(MPI_COMM_WORLD, "dx/rhoe = %g\n",
             (g.Lx / g.gdims[0]) / (g.vthe / g.wce));
  mpi_printf(MPI_COMM_WORLD, "L/debye = %g\n", g.L / (g.vthe / g.wpe));
  mpi_printf(MPI_COMM_WORLD, "dx/debye = %g\n",
             (g.Lx / g.gdims[0]) / (g.vthe / g.wpe));
  mpi_printf(MPI_COMM_WORLD, "n0 = %g\n", g.n0);
  mpi_printf(MPI_COMM_WORLD, "vthi/c = %g\n", g.vthi / g.c);
  mpi_printf(MPI_COMM_WORLD, "vthe/c = %g\n", g.vthe / g.c);
  mpi_printf(MPI_COMM_WORLD, "vdri/c = %g\n", g.vdri / g.c);
  mpi_printf(MPI_COMM_WORLD, "vdre/c = %g\n", g.vdre / g.c);
  mpi_printf(MPI_COMM_WORLD, "gdri = %g\n", g.gdri);
  mpi_printf(MPI_COMM_WORLD, "gdre = %g\n", g.gdre);
  //-------------------------------------------------------------

  // --- setup domain
  Grid_t::Real3 LL = {
    g.Lx_di, g.Ly_di,
    g.Lz_di}; // domain size (in d_i) !!!! This is important (jeff)
  // Grid_t::Real3 LL = {3. * 80., 1., 80.}; // domain size (in d_e)
  // Int3 gdims = {3 * 80, 1, 80};           // global number of grid points
  // Int3 np = {3*5, 1, 2};                // division into patches

  Grid_t::Domain domain{g.gdims, LL, -.5 * LL, g.np};
  // There was an issue with the conducting and reflective boundary conditions.
  // This returns continuity diff messages.
  // Both and each of them generate the discontinuity error

  //--------------------------------------------------------------------------------
  // This is in the case of yz geometry
  psc::grid::BC bc{
    {BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL,
     BND_FLD_PERIODIC}, // this is in the case of yz geometry
    {BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_PERIODIC},
    {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC},
    {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC}};
  //--------------------------------------------------------------------------------

  // -- setup particle kinds
  // last population ("i") is neutralizing
  Grid_t::Kinds kinds(N_MY_KINDS);
  kinds[MY_ION_UP] = {g.Zi, g.mi_me * g.Zi, "i_UP"};
  kinds[MY_ELECTRON_UP] = {-1., 1., "e_UP"};
  // kinds[MY_ION_BO] = {g.Zi, g.mi_me * g.Zi, "i_BO"};
  // kinds[MY_ELECTRON_BO] = {-1., 1., "e_BO"};

  g.di = sqrt(kinds[MY_ION_UP].m / kinds[MY_ION_UP].q);

  mpi_printf(MPI_COMM_WORLD, "de = %g, di = %g\n", 1., g.di);
  mpi_printf(MPI_COMM_WORLD, "lambda_De (background) = %g\n",
             sqrt(g.background_Te));

  // -- setup normalization
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = g.nppc;

  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  mpi_printf(MPI_COMM_WORLD, "dt = %g\n", dt);

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
// Langevin antenna

struct Langevin
{
  Langevin();

  void step(double dt);

  int Nk;
  double dB_bar;
  double omega_0;
  double gamma_0;
  double g_rate;
  double dBn;

  std::array<dcomp, 8> bn_k;
  std::array<Double3, 8> k;
  std::array<double, 8> k_per;

  Rng* rng_;
  int rank_;
};

Langevin::Langevin()
{
  Nk = 8;

  double L_per = sqrt(sqr(g.Lx) + sqr(g.Ly));
  dB_bar = 0.5 * g.b0 * L_per / g.Lz; // 0.5 * g.b0 * g.db_b0 * L_per/g.Lz ;

  dBn = dB_bar;
  // Checking the iterations, assuming B_bar remains constant,
  // then dBn is always dB0.

  omega_0 = 0.9 * (2. * M_PI * g.v_A /
                   g.Lz); // These are the values according to Daniel Groselj
  gamma_0 = 0.6 * omega_0;
  g_rate = 0.1; // This is chosen so omega_0 << g_rate << wpe;

  double k_x = 2. * M_PI / g.Lx;
  double k_y = 2. * M_PI / g.Ly;
  double k_z = 2. * M_PI / g.Lz;

  k[0] = {1. * k_x, 0. * k_y, 1. * k_z};
  k[1] = {1. * k_x, 0. * k_y, -1. * k_z};
  k[2] = {0. * k_x, 1. * k_y, 1. * k_z};
  k[3] = {0. * k_x, 1. * k_y, -1. * k_z};
  k[4] = {-1. * k_x, 0. * k_y, 1. * k_z};
  k[5] = {-1. * k_x, 0. * k_y, -1. * k_z};
  k[6] = {0. * k_x, -1. * k_y, 1. * k_z};
  k[7] = {0. * k_x, -1. * k_y, -1. * k_z};

  for (int n = 0; n < Nk; n++) {
    k_per[n] = std::sqrt(sqr(k[n][0]) + sqr(k[n][1]));
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  rngpool = RngPool_create();
  RngPool_seed(rngpool, rank_);
  rng_ = RngPool_get(rngpool, 0);
  if (rank_ == 0) {
    double rph_a = 0.;
    double rph_b = 2. * M_PI;

    //-------------------------------------------
    // const dcomp i(0.0,1.0);
    //-----------------------------------------------------------

    for (int n = 0; n < Nk; n++) {
      // I think this numbers need to change at each time
      double rand_ph = Rng_uniform(rng_, rph_a, rph_b);
      bn_k[n] = std::polar(dBn, rand_ph);
    }
  }
  MPI_Bcast(bn_k.data(), Nk, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}

void Langevin::step(double dt)
{
  const double ua = -0.5;
  const double ub = 0.5;

  // This is the time step that has to be calculated
  double delta_t_n = dt;

  //  double Cnp1 = 1.; // Since dBn = dB0 always, Cnp1=1 always.
  double Cnp1 = 1. + delta_t_n * g_rate * (dB_bar - dBn) / dB_bar;
  double dBnp1 = Cnp1 * dBn;

  dcomp unk[8];
  if (rank_ == 0) {
    for (int n = 0; n < Nk; n++) {
      unk[n] = Rng_uniform(rng_, ua, ub) + 1.i * Rng_uniform(rng_, ua, ub);
    }
  }
  MPI_Bcast(unk, Nk, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

  // This is an iterative formula

  dcomp bnp1_k[8];
  for (int n = 0; n < Nk; n++) {
    bnp1_k[n] =
      Cnp1 * bn_k[n] * std::exp(-(gamma_0 + omega_0 * 1.i) * delta_t_n) +
      dBnp1 * sqrt(12. * gamma_0 * delta_t_n) * unk[n];
  }

  // update n -> n + 1
  dBn = dBnp1;
  for (int n = 0; n < Nk; n++) {
    bn_k[n] = bnp1_k[n];
  }
};

struct Langevin* lng;

// ======================================================================
// initializeAlfven

void initializeAlfven(MfieldsAlfven& mflds, Langevin& lng)
{
  const auto& grid = mflds.grid();
  double ky = 2. * M_PI / grid.domain.length[1];

  //--------------------------------------------------------------------------------
  mpi_printf(grid.comm(), "**** Setting up Alfven fields...\n");

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto& patch = grid.patches[p];
    auto F = make_Fields3d<Dim>(mflds[p]);

    int n_ghosts = std::max(
      {mflds.ibn()[0], mflds.ibn()[1], mflds.ibn()[2]}); // FIXME, not pretty

    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      Int3 index{jx, jy, jz};
      auto crd_ec_z = Centering::getPos(patch, index, Centering::EC, 2);
      //--------------------------------------------------------------------------------
      // Calculate the vector potential
      dcomp Bext_x = 0., Bext_y = 0., Az = 0.;
      for (int n = 0; n < lng.Nk; n++) {
        Az += (lng.bn_k[n] *
               exp(1i * (lng.k[n][0] * crd_ec_z[0] + lng.k[n][1] * crd_ec_z[1] +
                         lng.k[n][2] * crd_ec_z[2]))) /
              lng.k_per[n];
      }
      F(PERT_AZ, jx, jy, jz) = Az.real();
    });
    //--------------------------------------------------------------------------------
    // Calculate the magnetic field components from the vector potential
    grid.Foreach_3d(1, 0, [&](int jx, int jy, int jz) {
      Int3 index{jx, jy, jz};

      dcomp Bext_x = 0., Bext_y = 0., Az = 0.;

      // Bx_{i, j+1/2, k+1/2} = (Az(i, j+1, k+1/2) - Az(i, j, k+1/2)) / dy
      // Bx_{i, j, k} =  (Az(i, j+1, k) - Az(i, j, k)) / dy

      // Bx_{i, j+1/2, k+1/2}
      // By_{i+1/2, j, k+1/2}
      // Bz_{i+1/2, j+1/2, k}
      F(PERT_HX, jx, jy, jz) =
        (F(PERT_AZ, jx, jy + 1, jz) - F(PERT_AZ, jx, jy, jz)) /
        (patch.y_nc(1) - patch.y_nc(0));
      F(PERT_HY, jx, jy, jz) =
        -(F(PERT_AZ, jx + 1, jy, jz) - F(PERT_AZ, jx, jy, jz)) /
        (patch.x_nc(1) - patch.x_nc(0));
      F(PERT_HZ, jx, jy, jz) = 0.;

      F(DIV_B, jx, jy, jz) =
        (F(PERT_HX, jx + 1, jy, jz) - F(PERT_HX, jx, jy, jz)) /
          (patch.x_nc(1) - patch.x_nc(0)) +
        (F(PERT_HY, jx, jy + 1, jz) - F(PERT_HY, jx, jy, jz)) /
          (patch.y_nc(1) - patch.y_nc(0)) +
        (F(PERT_HZ, jx, jy, jz + 1) - F(PERT_HZ, jx, jy, jz)) /
          (patch.z_nc(1) - patch.z_nc(0));
    });
    //--------------------------------------------------------------------------------
    // Calculate the external current componentfrom the magnetuc field
    grid.Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
      Int3 index{jx, jy, jz};

      double Jext_x = 0., Jext_y = 0., Jext_z = 0.;

      // Bx_{i, j+1/2, k+1/2}
      // By_{i+1/2, j, k+1/2}
      // Bz_{i+1/2, j+1/2, k}

      // This is where the current compomemts live
      //---------------------------------------------------------------------------
      // Jx_{i+1/2, j, k} = (Bz(i+1/2, j+1/2, k+1) - Bz(i+1/2, j-1/2, k+1)) / dy
      //                  - (By(i, j+1/2, k+1/2) - By(i, j+1/2, k-1/2)) / dz

      // Jy_{i+1, j+1/2, k} = -(Bz(i+1/2, j+1/2, k+1) - Bz(i-1/2, j+1/2, k+1)) /
      // dx
      //                  + (Bx(i+1, j+1/2, k+1/2) - Bx(i+1, j+1/2, k-1/2)) / dz

      // Jz_{i, j, k+1/2} = (By(i+1/2, j, k+1/2) - By(i-1/2, j, k+1/2)) / dx
      //                  - (Bx(i, j+1/2, k+1/2) - Bx(i, j-1/2, k+1/2)) / dy
      //--------------------------------------------------------------------------------
      // This is how it is implemented assuming that the right coordinates take
      // care of the 1/2's
      //--------------------------------------------------------------------------------
      // Jx_{i, j, k} =  (Bz(i, j, k) - Bz(i, j-1, k)) / dy
      //                  - (By(i, j, k) - By(i, j, k-1)) / dz

      // Jy_{i, j, k} =  -(Bz(i, j, k) - Bz(i-1, j, k)) / dx
      //                  + (Bx(i, j, k) - Bx(i, j, k-1)) / dz

      // Jz_{i, j, k} =  (By(i, j, k) - By(i-1, j, k)) / dx
      //                  - (Bx(i, j, k) - Bx(i, j-1, k)) / dy
      //--------------------------------------------------------------------------------

      F(PERT_JX_ext, jx, jy, jz) =
        (F(PERT_HZ, jx, jy, jz) - F(PERT_HZ, jx, jy - 1, jz)) /
          (patch.y_nc(1) - patch.y_nc(0)) -
        (F(PERT_HY, jx, jy, jz) - F(PERT_HY, jx, jy, jz - 1)) /
          (patch.z_nc(1) - patch.z_nc(0));

      F(PERT_JY_ext, jx, jy, jz) =
        (F(PERT_HX, jx, jy, jz) - F(PERT_HX, jx, jy, jz - 1)) /
          (patch.z_nc(1) - patch.z_nc(0)) -
        (F(PERT_HZ, jx, jy, jz) - F(PERT_HZ, jx - 1, jy, jz)) /
          (patch.x_nc(1) - patch.x_nc(0));

      F(PERT_JZ_ext, jx, jy, jz) =
        (F(PERT_HY, jx, jy, jz) - F(PERT_HY, jx - 1, jy, jz)) /
          (patch.x_nc(1) - patch.x_nc(0)) -
        (F(PERT_HX, jx, jy, jz) - F(PERT_HX, jx, jy - 1, jz)) /
          (patch.y_nc(1) - patch.y_nc(0));
    });
  }
  //--------------------------------------------------------------------------------

#if 0

  //------------------------------------------------------------------------------------------------------------
#endif
}

///***
// ======================================================================
// initializeParticles

void initializeParticles(SetupParticles<Mparticles>& setup_particles,
                         Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts,
                         MfieldsAlfven& mflds_alfven)
{
  //***
  // -- set particle initial condition
  partitionAndSetupParticlesGeneral(
    setup_particles, balance, grid_ptr, mprts,
    [&](int kind, Double3 crd, int p, Int3 idx, psc_particle_npt& npt) {
      double x = crd[0], y = crd[1], z = crd[2];
      switch (kind) {

        //--------------------------------------------------------------------------------
        /*** // This is a block for the yz configuration using a single
        poppulation double psi; if (y<=g.L && y>=0.) psi=1.; else if (y<0. &&
        y>=-g.L) psi=1.; else psi=0.; case 0: //Ion drifting up npt.n = g.n0 *
        (g.nb_n0 + (1 / sqr(cosh(y / g.L))) ) ; npt.T[0] = g.Ti * psi + g.Tbi;
          npt.T[1] = g.Ti * psi + g.Tbi;
          npt.T[2] = g.Ti * psi + g.Tbi;
          npt.p[0] = g.udri * psi;
          npt.p[1] = 0.;
          npt.p[2] = 0.;
          npt.kind = MY_ION_UP;
          break;
        case 1: //Electron drifting up
          npt.n =  g.n0 * (g.nb_n0 + (1 / sqr(cosh(y / g.L))) ) ;
          npt.T[0] = g.Te * psi + g.Tbe;
          npt.T[1] = g.Te * psi + g.Tbe;
          npt.T[2] = g.Te * psi + g.Tbe;
          npt.p[0] = g.udre * psi ;
          npt.p[1] = 0.;
          npt.p[2] = 0.;
          npt.kind = MY_ELECTRON_UP;
          break;
        ***/
        //--------------------------------------------------------------------------------

        //--------------------------------------------------------------------------------
        //*** // This is a block for the yz configuration using two
        // popullations
        case 0: // Ion drifting up
          npt.n = g.n0 * (g.nb_n0 + (1 / sqr(cosh(y / g.L))));
          npt.T[0] = g.Ti;
          npt.T[1] = g.Ti;
          npt.T[2] = g.Ti;
          npt.p[0] = g.udri;
          npt.p[1] = 0.;
          npt.p[2] = 0.;
          npt.kind = MY_ION_UP;
          break;
        case 1: // Electron drifting up
          npt.n = g.n0 * (g.nb_n0 + (1 / sqr(cosh(y / g.L))));
          npt.T[0] = g.Te;
          npt.T[1] = g.Te;
          npt.T[2] = g.Te;
          npt.p[0] = g.udre;
          npt.p[1] = 0.;
          npt.p[2] = 0.;
          npt.kind = MY_ELECTRON_UP;
          break;
        case 2: // Ion background up
          npt.n = g.n0 * g.nb_n0;
          npt.T[0] = g.Tbi;
          npt.T[1] = g.Tbi;
          npt.T[2] = g.Tbi;
          npt.p[0] = 0.;
          npt.p[1] = 0.;
          npt.p[2] = 0.;
          npt.kind = MY_ION_UP;
          break;
        case 3: // Electron background up
          npt.n = g.n0 * g.nb_n0;
          npt.T[0] = g.Tbe;
          npt.T[1] = g.Tbe;
          npt.T[2] = g.Tbe;
          npt.p[0] = 0.;
          npt.p[1] = 0.;
          npt.p[2] = 0.;
          npt.kind = MY_ELECTRON_UP;
          break;
          /***/
          //--------------------------------------------------------------------------------

        default: assert(0);
      }
    });
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds, MfieldsAlfven& mflds_alfven)
{
  setupFieldsGeneral(
    mflds, [&](int m, Int3 idx, int p, double crd[3]) -> MfieldsState::real_t {
      double x = crd[0], y = crd[1], z = crd[2];
      switch (m) {
        //--------------------------------------------------------------------------------
        // case HY: return mflds_alfven(PERT_HY, idx[0], idx[1], idx[2], p);
        // default: return 0.;
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        ///*** This is the magnetic field in the case yz geometry
        // case HX:
        //   return mflds_alfven(PERT_HX, idx[0], idx[1], idx[2], p); // 0.
        // case HY:
        //   return mflds_alfven(PERT_HY, idx[0], idx[1], idx[2], p);
        //   // return 0. + g.dby * sin(2. * M_PI * (z - 0.5 * g.Lz) / g.Lz) *
        //   //               cos(M_PI * y / g.Ly); // + dB_azT //In  the case
        //   of
        //   //               yz geometry
        // case HZ:
        //   //return mflds_alfven(PERT_HZ, idx[0], idx[1], idx[2], p);
        //   // return g.b0 * tanh(y / g.L) +
        //   //        g.dbz * cos(2. * M_PI * (z - 0.5 * g.Lz) / g.Lz) *
        //   //          sin(M_PI * y / g.Ly);
        //   return mflds_alfven(DIV_B, idx[0], idx[1], idx[2], p); // To check
        //   the divB
        case HX: return mflds_alfven(PERT_JX_ext, idx[0], idx[1], idx[2], p);
        case HY: return mflds_alfven(PERT_JY_ext, idx[0], idx[1], idx[2], p);
        case HZ: return mflds_alfven(PERT_JZ_ext, idx[0], idx[1], idx[2], p);
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
  psc_params.balance_interval = 180;
  Balance balance{3};

  // -- Sort
  psc_params.sort_interval = 10;

  // -- Collision
  int collision_interval = 0;
  double collision_nu = 1e-10;
  //    3.76 * std::pow(g.target_Te, 2.) / g.Zi / g.lambda0;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.continuity_every_step = -2;
  checks_params.continuity_dump_always = false;
  checks_params.continuity_threshold = 1e-4;
  checks_params.continuity_verbose = true;

  checks_params.gauss_every_step = -2;
  checks_params.gauss_dump_always = false;
  checks_params.gauss_threshold = 1e-4;
  checks_params.gauss_verbose = true;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = 10;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output
  //
  // FIXME, this really is too complicated and not very flexible

  // -- output fields
  OutputFieldsItemParams outf_item_params{};
  OutputFieldsParams outf_params{};
  outf_item_params.pfield_interval = 45;
  outf_item_params.tfield_interval = -10;

  outf_params.fields = outf_item_params;
  outf_params.moments = outf_item_params;
  OutputFields<MfieldsState, Mparticles, Dim, Writer> outf{grid, outf_params};

  // -- output particles
  OutputParticlesParams outp_params{};
  outp_params.every_step = -100;
  outp_params.data_dir = ".";
  outp_params.basename = "prt";
  OutputParticles outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // Set up objects specific to the TurbHarrisxz case

  SetupParticles<Mparticles> setup_particles(grid);
  setup_particles.fractional_n_particles_per_cell = true;
  setup_particles.neutralizing_population = MY_ION_UP; // It has to be the ions

  // ----------------------------------------------------------------------
  // setup initial conditions

  lng = new Langevin();
  for (int n = 0; n < 10; n++) {
    lng->step(.1);
  }

  if (read_checkpoint_filename
        .empty()) { // This is the block which is returning the
    MfieldsAlfven mflds_alfven(grid, N_PERT, grid.ibn);
    initializeAlfven(mflds_alfven, *lng);
    initializeParticles(setup_particles, balance, grid_ptr, mprts,
                        mflds_alfven);
    initializeFields(mflds, mflds_alfven);
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
