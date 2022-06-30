
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
  MY_ELECTRON,
  MY_ION,
  N_MY_KINDS,
};

enum
{
  PERT_HX,
  PERT_HY,
  PERT_HZ,
  PERT_VX,
  PERT_VY,
  PERT_VZ,
  N_PERT,
};

// ======================================================================
// PscTurbHarrisxzParams

struct PscTurbHarrisxzParams
{
  double BB;
  double Zi;
  double mass_ratio;
  double lambda0;

  double background_n;
  double background_Te;
  double background_Ti;

  // The following parameters are calculated from the above / and other
  // information

  double d_i;
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

using Dim = dim_xyz;

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
  psc_params.nmax = 2001;
  psc_params.cfl = 0.75;
  psc_params.write_checkpoint_every_step = -100; //This is not working
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
  g.BB = 1.;
  g.Zi = 1.;
  g.mass_ratio = 100.;
  g.lambda0 = 20.;

  double vA_over_c_ = .1; //Why 0.1?? 
  double amplitude_ = .5;     
  double beta_e_par_ = 1.; //Ques_jeff what beta 0.1, 1??
  double beta_i_par_ = 1.;
  double Ti_perp_over_Ti_par_ = 1.;
  double Te_perp_over_Te_par_ = 1.;

  double B0 = vA_over_c_;     
  double debye_length_ = 1.* vA_over_c_ *sqrt(beta_i_par_ / 2. );
  double Te_par = beta_e_par_ * sqr(B0) / 2.;
  double Te_perp = Te_perp_over_Te_par_ * Te_par;
  double Ti_par = beta_i_par_ * sqr(B0) / 2.;
  double Ti_perp = Ti_perp_over_Ti_par_ * Ti_par;

  g.background_n = 1.;
  g.background_Te = Ti_par;
  g.background_Ti = Te_par;
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
  // --- setup domain
  //Grid_t::Real3 LL = {1., 80., 3. * 80.}; // domain size (in d_e)
  //Int3 gdims = {1, 80, 3 * 80};           // global number of grid points
  //Int3 np = {1, 2, 3 * 5};                // division into patches
  Grid_t::Real3 LL = {2. * M_PI, 2. * M_PI, 2. * M_PI}; // domain size (in d_e) 
  Int3 gdims = {32, 32, 32};           // global number of grid points
  Int3 np = {2, 2, 2};                // division into patches

  Grid_t::Domain domain{gdims, LL, -.5 * LL, np};

  psc::grid::BC bc{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  // -- setup particle kinds
  // last population ("i") is neutralizing
  Grid_t::Kinds kinds(N_MY_KINDS);
  kinds[MY_ION] = {g.Zi, g.mass_ratio * g.Zi, "i"};
  kinds[MY_ELECTRON] = {-1., 1., "e"};

  g.d_i = sqrt(kinds[MY_ION].m / kinds[MY_ION].q);

  mpi_printf(MPI_COMM_WORLD, "d_e = %g, d_i = %g\n", 1., g.d_i);
  mpi_printf(MPI_COMM_WORLD, "lambda_De (background) = %g\n",
   sqrt(g.background_Te));

  // -- setup normalization
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 20;

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
// initializeAlfven

void initializeAlfven(MfieldsAlfven& mflds)
{  
  const auto& grid = mflds.grid();
  std::array<Double3, 8> k;
  std::array<double, 8> mn_per;
  std::array<double, 8> mn_par;
  std::array<double, 8> k_a_per;
  std::array<double, 8> k_a_par;
  std::array<double, 8> phi;
  std::array<double, 8> phase;
  std::array<double, 8> Amp;
  std::array<double, 8> dB_ax;
  std::array<double, 8> dB_ay;
  std::array<double, 8> dv_ax;
  std::array<double, 8> dv_ay;
  //-----------------------------------------------------------------------
  double vA_over_c_ = .1; 
  double B0 = vA_over_c_; 
  double beta_e_par_ = 1.; 
  double beta_i_par_ = 1.;

  //double x = crd[0], y = crd[1], z = crd[2];    
  double Lx = grid.domain.length[0];  
  double Ly = grid.domain.length[1];       
  double Lz = grid.domain.length[2];

  double k_x = 2. * M_PI / Lx;
  double k_y = 2. * M_PI / Ly;
  double k_z = 2. * M_PI / Lz;
  //-------------------------------------------------------------------
  double rho_i = sqrt(beta_i_par_)*1.; //check how to set this as di !!
  double sp = 1./3.;  //Spectral index for the alfven wave (1/3 for AW and 2/3 for KAW)
  int p = 1;// fraction of B0 so dB=(1/p)B0
  double crit_fact = 0.1;   //critical balance /normalization coefficient       
  double C1; // normalization factor
  
  int m_per=2; //modes in the perpendicular directions
  int m_par=1; //modes in the parallel direction acording to critical balance    
  int Nk = 8;

  // This part didn't work. there is a conflict that I don't understand yet
  //---------------------------------------------------------------------
  //Rng* rng_;
  //int rank_;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  //rngpool = RngPool_create();
  //RngPool_seed(rngpool, rank_);
  //rng_ = RngPool_get(rngpool, 0);
  //for (int n = 0; n < Nk; n++) {
  //        phase[n] = Rng_uniform(rng_, rph_a, rph_b); // random phase
  //        }  
  //---------------------------------------------------------------------
 
  //if (rank_ == 0) {
  double rph_a = 0.;
  double rph_b = 2. * M_PI;
  //----------------------------------------------------------------------
  double dB_axT, dB_ayT, dv_axT, dv_ayT; 

  for (int n = 0; n < Nk; n++) {
    mn_per[n] = 1. * m_per;
    mn_par[n] = 1. * m_par;
  }

  k[0] = {1. * k_x * mn_per[0], 0. * k_y * mn_per[0], 1. * k_z * mn_par[0]};
  k[1] = {0. * k_x * mn_per[1], 1. * k_y * mn_per[1], -1. * k_z * mn_par[1]};
  k[2] = {-1. * k_x * mn_per[2], 0. * k_y * mn_per[2], 1. * k_z * mn_par[2]};
  k[3] = {0. * k_x * mn_per[3], -1. * k_y * mn_per[3], -1. * k_z * mn_par[3]};
  k[4] = {1. * k_x * mn_per[4], 1. * k_y * mn_per[4], 1. * k_z * mn_par[4]};
  k[5] = {-1. * k_x * mn_per[5], 1. * k_y * mn_per[5], -1. * k_z * mn_par[5]};
  k[6] = {-1. * k_x * mn_per[6], -1. * k_y * mn_per[6], 1. * k_z * mn_par[6]};
  k[7] = {1. * k_x * mn_per[7], -1. * k_y * mn_per[7], -1. * k_z * mn_par[7]};

  phase={0.987*2.*M_PI, 0.666*2.*M_PI, 0.025*2.*M_PI, 0.954*2.*M_PI,
         0.781*2.*M_PI, 0.846*2.*M_PI, 0.192*2.*M_PI, 0.778*2.*M_PI}; // random phases   

  //-----------------------------------------------------------------------
  mpi_printf(grid.comm(), "**** Setting up Alfven fields...\n");

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto& patch = grid.patches[p];
    auto F = make_Fields3d<dim_xyz>(mflds[p]);

    int n_ghosts = std::max(
      {mflds.ibn()[0], mflds.ibn()[1], mflds.ibn()[2]}); // FIXME, not pretty

    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      Int3 index{jx, jy, jz};
      auto crd_fc = Centering::getPos(patch, index, Centering::FC);
      auto crd_cc = Centering::getPos(patch, index, Centering::CC); 

      for (int n = 0; n < Nk; n++) {
        k_a_per[n] = sqrt(sqr(k[n][0]) + sqr(k[n][1])); // k_per
        k_a_par[n] = k[n][2]; //k_par

        phi[n] = n * M_PI/2.; // Polarization angle
        //phase[n] = Rng_uniform(rng_, rph_a, rph_b); // random phase

        Amp[n] = B0 * crit_fact *pow((k_a_per[n]), -sp);//According to the critical balance Amp[i]=B0*crit_fact * K_a_perp[i]^(-sp);
        dB_ax[n] = -Amp[n] * cos ( k[n][0] * crd_fc[0] + k[n][1] * crd_fc[1] + k[n][2] * crd_fc[2] + phase[n] ) * sin(phi[n]) ;
        dB_ay[n] =  Amp[n] * cos ( k[n][0] * crd_fc[0] + k[n][1] * crd_fc[1] + k[n][2] * crd_fc[2] + phase[n] ) * cos(phi[n]) ;
        
        dv_ax[n] = -Amp[n] * cos ( k[n][0] * crd_cc[0] + k[n][1] * crd_cc[1] + k[n][2] * crd_cc[2] + phase[n] ) * sin(phi[n]) ;
        dv_ay[n] =  Amp[n] * cos ( k[n][0] * crd_cc[0] + k[n][1] * crd_cc[1] + k[n][2] * crd_cc[2] + phase[n] ) * cos(phi[n]) ;        
      }

      dB_axT = dB_ax[0] + dB_ax[1] + dB_ax[2] + dB_ax[3] + dB_ax[4] + dB_ax[5] + dB_ax[6] + dB_ax[7];
      dB_ayT = dB_ay[0] + dB_ay[1] + dB_ay[2] + dB_ay[3] + dB_ay[4] + dB_ay[5] + dB_ay[6] + dB_ay[7];
      
      dv_axT = -dv_ax[0] + dv_ax[1] - dv_ax[2] + dv_ax[3] - dv_ax[4] + dv_ax[5] - dv_ax[6] + dv_ax[7];
      dv_ayT = -dv_ay[0] + dv_ay[1] - dv_ay[2] + dv_ay[3] - dv_ay[4] + dv_ay[5] - dv_ay[6] + dv_ay[7];

      C1 = B0/sqrt(sqr(Amp[0])+sqr(Amp[1])+sqr(Amp[2])+sqr(Amp[3])+sqr(Amp[4])+sqr(Amp[5])+sqr(Amp[6])+sqr(Amp[7])); //    

      F(PERT_HX, jx, jy, jz) = C1 * dB_axT;//g.BB + .1 * sin(ky * crd_fc[1]);
      F(PERT_HY, jx, jy, jz) = C1 * dB_ayT;//g.BB + .1 * sin(ky * crd_fc[1]);
      F(PERT_HZ, jx, jy, jz) = B0;//g.BB + .1 * sin(ky * crd_fc[1]);
      
      F(PERT_VX, jx, jy, jz) = C1 * dv_axT;//-.1 * sin(ky * crd_cc[1]);
      F(PERT_VY, jx, jy, jz) = C1 * dv_ayT;//-.1 * sin(ky * crd_cc[1]);
      F(PERT_VZ, jx, jy, jz) = 0.;//g.BB + .1 * sin(ky * crd_fc[1]);
    });
  }
  //}
}

// ======================================================================
// initializeParticles

    void initializeParticles(SetupParticles<Mparticles>& setup_particles,
     Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts,
     MfieldsAlfven& mflds_alfven)
    {
  // -- set particle initial condition
      partitionAndSetupParticlesGeneral(
        setup_particles, balance, grid_ptr, mprts,
        [&](int kind, Double3 crd, int p, Int3 idx, psc_particle_npt& npt) {
          switch (kind) {
            case MY_ION:
            npt.n = g.background_n;
            npt.T[0] = g.background_Ti;
            npt.T[1] = g.background_Ti;
            npt.T[2] = g.background_Ti;
            npt.p[0] = mflds_alfven(PERT_VX, idx[0], idx[1], idx[2], p);
            npt.p[1] = mflds_alfven(PERT_VY, idx[0], idx[1], idx[2], p);
            npt.p[2] = mflds_alfven(PERT_VZ, idx[0], idx[1], idx[2], p);
            break;
            case MY_ELECTRON:
            npt.n = g.background_n;
            npt.T[0] = g.background_Te;
            npt.T[1] = g.background_Te;
            npt.T[2] = g.background_Te;
            npt.p[0] = mflds_alfven(PERT_VX, idx[0], idx[1], idx[2], p);
            npt.p[1] = mflds_alfven(PERT_VY, idx[0], idx[1], idx[2], p);
            npt.p[2] = mflds_alfven(PERT_VZ, idx[0], idx[1], idx[2], p);
            break;
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
          switch (m) {
            case HX: return mflds_alfven(PERT_HX, idx[0], idx[1], idx[2], p);
            case HY: return mflds_alfven(PERT_HY, idx[0], idx[1], idx[2], p); 
            case HZ: return mflds_alfven(PERT_HZ, idx[0], idx[1], idx[2], p);                   
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
      psc_params.balance_interval = 200;
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
      checks_params.continuity_every_step = 10;
      checks_params.continuity_dump_always = false;
      checks_params.continuity_threshold = 1e-4;
      checks_params.continuity_verbose = true;

      checks_params.gauss_every_step = 10;
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
      outf_item_params.pfield_interval = 50;
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
  //setup_particles.neutralizing_population = MY_ION;

  // ----------------------------------------------------------------------
  // setup initial conditions

      if (read_checkpoint_filename.empty()) {
        MfieldsAlfven mflds_alfven(grid, N_PERT, grid.ibn);
        initializeAlfven(mflds_alfven);
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
