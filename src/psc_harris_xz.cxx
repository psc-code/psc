
#define VPIC 1

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "../libpsc/vpic/fields_item_vpic.hxx"
#include "../libpsc/vpic/setup_fields_vpic.hxx"
#include "OutputFieldsDefault.h"
#include "psc_config.hxx"

#include "rngpool_iface.h"

extern Grid* vgrid; // FIXME

static RngPool* rngpool; // FIXME, should be member (of struct psc, really)

// FIXME, helper should go somewhere...

static inline double trunc_granular(double a, double b)
{
  return b * (int)(a / b);
}

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

#ifdef VPIC
#ifdef DO_VPIC
using PscConfig = PscConfigVpicWrap;
#else
using PscConfig = PscConfigVpicPsc;
#endif
#else
using PscConfig = PscConfig1vbecSingle<dim_xz>;
#endif

// ----------------------------------------------------------------------

using MfieldsState = typename PscConfig::MfieldsState;
#ifdef VPIC
using MaterialList = typename MfieldsState::MaterialList;
#endif
using Mparticles = typename PscConfig::Mparticles;
using Balance = typename PscConfig::Balance;
using Collision = typename PscConfig::Collision;
using Checks = typename PscConfig::Checks;
using Marder = typename PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

// FIXME!
MfieldsC evalMfields(const MfieldsState& _exp)
{
  auto& exp = const_cast<MfieldsState&>(_exp);
  MfieldsC mflds{exp.grid(), exp.n_comps(), exp.ibn()};

  for (int p = 0; p < mflds.n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(mflds[p]);
    auto _exp = make_Fields3d<dim_xyz>(exp[p]);
    for (int m = 0; m < exp.n_comps(); m++) {
      mflds.Foreach_3d(0, 0, [&](int i, int j, int k) {
        flds(m, i, j, k) = _exp(m, i, j, k);
      });
    }
  }
  return mflds;
}

// ======================================================================
// PscHarrisParams

struct PscHarrisParams
{
  double L_di;   // Sheet thickness / ion inertial length
  double Ti_Te;  // Ion temperature / electron temperature
  double nb_n0;  // background plasma density
  double Tbe_Te; // Ratio of background T_e to Harris T_e
  double Tbi_Ti; // Ratio of background T_i to Harris T_i

  double bg; // Guide field
  double theta;

  double Lpert_Lx; // wavelength of perturbation in terms of Lx
  double dbz_b0;   // perturbation in Bz relative to B0
  double nppc;     // Average number of macro particle per cell per species
  bool open_bc_x;  // Flag to signal we want to do open boundary condition in x
  bool
    driven_bc_z; // Flag to signal we want to do driven boundary condition in z

  // FIXME, not really harris-specific
  double wpedt_max;

  double wpe_wce; // electron plasma freq / electron cyclotron freq
  double mi_me;   // Ion mass / electron mass

  double Lx_di, Ly_di, Lz_di; // Size of box in d_i

  int ion_sort_interval;
  int electron_sort_interval;

  double taui;                     // simulation wci's to run
  double t_intervali;              // output interval in terms of 1/wci
  double output_particle_interval; // particle output interval in terms of 1/wci

  double overalloc; // Overallocation factor (> 1) for particle arrays

  Int3 gdims;
  Int3 np;
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
// setupHarrisParams()

void setupHarrisParams()
{
  g.wpedt_max = .36;
  g.wpe_wce = 2.;
  g.mi_me = 25.;

  g.Lx_di = 40.;
  g.Ly_di = 1.;
  g.Lz_di = 10.;

  g.electron_sort_interval = 25;
  g.ion_sort_interval = 25;

  g.taui = 40.;
  g.t_intervali = 1.;
  g.output_particle_interval = 10.;

  g.overalloc = 2.;

  g.gdims = {512, 1, 128};
  g.np = {4, 1, 1};

  g.L_di = .5;
  g.Ti_Te = 5.;
  g.nb_n0 = .05;
  g.Tbe_Te = .333;
  g.Tbi_Ti = .333;

  g.bg = 0.;
  g.theta = 0.;

  g.Lpert_Lx = 1.;
  g.dbz_b0 = .03;
  g.nppc = 100;
  g.open_bc_x = false;
  g.driven_bc_z = false;
}

// ======================================================================
// globals_physics
//
// FIXME rename / merge?

struct globals_physics
{
  double ec;
  double me;
  double c;
  double eps0;
  double de;

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
  double dbx;        // Perturbation in Bz relative to Bo (Only change here)
  double dbz;        // Set Bx perturbation so that div(B) = 0
  double tanhf;

  double Ne;       // Total number of macro electrons
  double Ne_sheet; // Number of macro electrons in Harris sheet
  double weight_s; // Charge per macro electron
  double vthe;     // Electron thermal velocity
  double vthi;     // Ion thermal velocity
  double vdre;     // Electron drift velocity
  double vdri;     // Ion drift velocity
  double gdri;     // gamma of ion drift frame
  double gdre;     // gamma of electron drift frame
  double udri;     // 4-velocity of ion drift frame
  double udre;     // 4-velocity of electron drift frame

  double Ne_back;  // Number of macro electrons in background
  double weight_b; // Charge per macro electron
  double vtheb;    // normalized background e thermal vel.
  double vthib;    // normalized background ion thermal vel.

  int n_global_patches;

  // ----------------------------------------------------------------------
  // ctor

  // FIXME, do we want to keep this?
  globals_physics() {}

  globals_physics(const PscHarrisParams& p)
  {
    assert(p.np[2] <= 2); // For load balance, keep "1" or "2" for Harris sheet

    // FIXME, the general normalization stuff should be shared somehow

    // use natural PIC units
    ec = 1;   // Charge normalization
    me = 1;   // Mass normalization
    c = 1;    // Speed of light
    de = 1;   // Length normalization (electron inertial length)
    eps0 = 1; // Permittivity of space

    // derived quantities
    mi = me * p.mi_me; // Ion mass

    double Te =
      me * sqr(c) /
      (2. * eps0 * sqr(p.wpe_wce) * (1. + p.Ti_Te)); // Electron temperature
    double Ti = Te * p.Ti_Te;                        // Ion temperature
    vthe = sqrt(Te / me);             // Electron thermal velocity
    vthi = sqrt(Ti / mi);             // Ion thermal velocity
    vtheb = sqrt(p.Tbe_Te * Te / me); // normalized background e thermal vel.
    vthib = sqrt(p.Tbi_Ti * Ti / mi); // normalized background ion thermal vel.
    wci = 1. / (p.mi_me * p.wpe_wce); // Ion cyclotron frequency
    wce = wci * p.mi_me;              // Electron cyclotron freqeuncy
    wpe = wce * p.wpe_wce;            // electron plasma frequency
    wpi = wpe / sqrt(p.mi_me);        // ion plasma frequency
    di = c / wpi;                     // ion inertial length
    L = p.L_di * di;                  // Harris sheet thickness
    rhoi_L = sqrt(p.Ti_Te / (1. + p.Ti_Te)) / p.L_di;
    v_A = (wci / wpi) / sqrt(p.nb_n0); // based on nb

    Lx = p.Lx_di * di; // size of box in x dimension
    Ly = p.Ly_di * di; // size of box in y dimension
    Lz = p.Lz_di * di; // size of box in z dimension

    b0 = me * c * wce / ec; // Asymptotic magnetic field strength
    n0 = me * eps0 * wpe * wpe / (ec * ec); // Peak electron (ion) density
    vdri = 2 * c * Ti / (ec * b0 * L);      // Ion drift velocity
    vdre = -vdri / (p.Ti_Te);               // electron drift velocity

    n_global_patches = p.np[0] * p.np[1] * p.np[2];
    double Npe_sheet =
      2 * n0 * Lx * Ly * L * tanh(0.5 * Lz / L);   // N physical e's in sheet
    double Npe_back = p.nb_n0 * n0 * Ly * Lz * Lx; // N physical e's in backgrnd
    double Npe = Npe_sheet + Npe_back;
    Ne = p.nppc * p.gdims[0] * p.gdims[1] *
         p.gdims[2]; // total macro electrons in box
    Ne_sheet = Ne * Npe_sheet / Npe;
    Ne_back = Ne * Npe_back / Npe;
    Ne_sheet = trunc_granular(
      Ne_sheet, n_global_patches); // Make it divisible by # subdomains
    Ne_back = trunc_granular(
      Ne_back, n_global_patches); // Make it divisible by # subdomains
    Ne = Ne_sheet + Ne_back;
    weight_s = ec * Npe_sheet / Ne_sheet; // Charge per macro electron
    weight_b = ec * Npe_back / Ne_back;   // Charge per macro electron

    gdri = 1. / sqrt(1. - sqr(vdri) / sqr(c)); // gamma of ion drift frame
    gdre = 1. / sqrt(1. - sqr(vdre) / sqr(c)); // gamma of electron drift frame
    udri = vdri * gdri;                        // 4-velocity of ion drift frame
    udre = vdre * gdre; // 4-velocity of electron drift frame
    tanhf = tanh(0.5 * Lz / L);
    Lpert = p.Lpert_Lx * Lx; // wavelength of perturbation
    dbz = p.dbz_b0 * b0; // Perturbation in Bz relative to Bo (Only change here)
    dbx = -dbz * Lpert / (2. * Lz); // Set Bx perturbation so that div(B) = 0
  }
};

globals_physics phys;

// ======================================================================
// setupGrid
//
// This helper function is responsible for setting up the "Grid",
// which is really more than just the domain and its decomposition, it
// also encompasses PC normalization parameters, information about the
// particle kinds, etc.

Grid_t* setupGrid()
{
  auto comm = MPI_COMM_WORLD;

  // --- set up domain

  auto domain = Grid_t::Domain{g.gdims,
                               {phys.Lx, phys.Ly, phys.Lz},
                               {0., -.5 * phys.Ly, -.5 * phys.Lz},
                               g.np};

  mpi_printf(comm, "Conducting fields on Z-boundaries\n");
  mpi_printf(comm, "Reflect particles on Z-boundaries\n");
  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_REFLECTING},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_REFLECTING}};
  if (g.open_bc_x) {
    mpi_printf(comm, "Absorbing fields on X-boundaries\n");
    bc.fld_lo[0] = BND_FLD_ABSORBING;
    bc.fld_hi[0] = BND_FLD_ABSORBING;
    mpi_printf(comm, "Absorb particles on X-boundaries\n");
    bc.prt_lo[1] = BND_PRT_ABSORBING;
    bc.prt_hi[1] = BND_PRT_ABSORBING;
  }

  if (g.driven_bc_z) {
    mpi_printf(comm, "Absorb particles on Z-boundaries\n");
    bc.prt_lo[2] = BND_PRT_ABSORBING;
    bc.prt_hi[2] = BND_PRT_ABSORBING;
  }

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-phys.ec, phys.me, "e"};
  kinds[KIND_ION] = {phys.ec, phys.mi, "i"};

  // determine the time step
  double dg = courant_length(domain);
  double dt = psc_params.cfl * dg / phys.c; // courant limited time step
  if (phys.wpe * dt > g.wpedt_max) {
    dt =
      g.wpedt_max / phys.wpe; // override timestep if plasma frequency limited
  }

  assert(phys.c == 1. && phys.eps0 == 1.);
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;
  auto norm = Grid_t::Normalization{norm_params};

#ifdef VPIC
  Int3 ibn = {1, 1, 1};
#else
  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();
#endif

  auto grid_ptr = new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
  vpic_define_grid(*grid_ptr);
  return grid_ptr;
}

// ----------------------------------------------------------------------
// setupMaterials

void setupMaterials(MaterialList& material_list)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  mpi_printf(comm, "Setting up materials.\n");

  // -- set up MaterialList
  vpic_define_material(material_list, "vacuum", 1., 1., 0., 0.);
#if 0
  struct material *resistive =
    vpic_define_material(material_list, "resistive", 1., 1., 1., 0.);
#endif

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.
}

// ----------------------------------------------------------------------
// vpic_setup_species
//
// FIXME, half-redundant to the PSC species setup

void vpic_setup_species(Mparticles& mprts)
{
  mpi_printf(mprts.grid().comm(), "Setting up species.\n");
  double nmax = g.overalloc * phys.Ne / phys.n_global_patches;
  double nmovers = .1 * nmax;
  double sort_method = 1; // 0=in place and 1=out of place

  mprts.define_species("electron", -phys.ec, phys.me, nmax, nmovers,
                       g.electron_sort_interval, sort_method);
  mprts.define_species("ion", phys.ec, phys.mi, nmax, nmovers,
                       g.ion_sort_interval, sort_method);
}

// ----------------------------------------------------------------------
// setup_log

void setup_log(const Grid_t& grid)
{
  MPI_Comm comm = grid.comm();

  mpi_printf(comm, "***********************************************\n");
  mpi_printf(comm, "* Topology: %d x %d x %d\n", g.np[0], g.np[1], g.np[2]);
  mpi_printf(comm, "tanhf    = %g\n", phys.tanhf);
  mpi_printf(comm, "L_di     = %g\n", g.L_di);
  mpi_printf(comm, "rhoi/L   = %g\n", phys.rhoi_L);
  mpi_printf(comm, "Ti/Te    = %g\n", g.Ti_Te);
  mpi_printf(comm, "nb/n0    = %g\n", g.nb_n0);
  mpi_printf(comm, "wpe/wce  = %g\n", g.wpe_wce);
  mpi_printf(comm, "mi/me    = %g\n", g.mi_me);
  mpi_printf(comm, "theta    = %g\n", g.theta);
  mpi_printf(comm, "Lpert/Lx = %g\n", g.Lpert_Lx);
  mpi_printf(comm, "dbz/b0   = %g\n", g.dbz_b0);
  mpi_printf(comm, "taui     = %g\n", g.taui);
  mpi_printf(comm, "t_intervali = %g\n", g.t_intervali);
  mpi_printf(comm, "num_step = %d\n", psc_params.nmax);
  mpi_printf(comm, "Lx/di = %g\n", phys.Lx / phys.di);
  mpi_printf(comm, "Lx/de = %g\n", phys.Lx / phys.de);
  mpi_printf(comm, "Ly/di = %g\n", phys.Ly / phys.di);
  mpi_printf(comm, "Ly/de = %g\n", phys.Ly / phys.de);
  mpi_printf(comm, "Lz/di = %g\n", phys.Lz / phys.di);
  mpi_printf(comm, "Lz/de = %g\n", phys.Lz / phys.de);
  mpi_printf(comm, "nx = %d\n", g.gdims[0]);
  mpi_printf(comm, "ny = %d\n", g.gdims[1]);
  mpi_printf(comm, "nz = %d\n", g.gdims[2]);
  mpi_printf(comm, "n_global_patches = %d\n", phys.n_global_patches);
  mpi_printf(comm, "nppc = %g\n", g.nppc);
  mpi_printf(comm, "b0 = %g\n", phys.b0);
  mpi_printf(comm, "v_A (based on nb) = %g\n", phys.v_A);
  mpi_printf(comm, "di = %g\n", phys.di);
  mpi_printf(comm, "Ne = %g\n", phys.Ne);
  mpi_printf(comm, "Ne_sheet = %g\n", phys.Ne_sheet);
  mpi_printf(comm, "Ne_back = %g\n", phys.Ne_back);
  mpi_printf(comm, "total # of particles = %g\n", 2 * phys.Ne);
  mpi_printf(comm, "dt*wpe = %g\n", phys.wpe * grid.dt);
  mpi_printf(comm, "dt*wce = %g\n", phys.wce * grid.dt);
  mpi_printf(comm, "dt*wci = %g\n", phys.wci * grid.dt);
  mpi_printf(comm, "dx/de = %g\n", phys.Lx / (phys.de * g.gdims[0]));
  mpi_printf(comm, "dy/de = %g\n", phys.Ly / (phys.de * g.gdims[1]));
  mpi_printf(comm, "dz/de = %g\n", phys.Lz / (phys.de * g.gdims[2]));
  mpi_printf(comm, "dx/rhoi = %g\n",
             (phys.Lx / g.gdims[0]) / (phys.vthi / phys.wci));
  mpi_printf(comm, "dx/rhoe = %g\n",
             (phys.Lx / g.gdims[0]) / (phys.vthe / phys.wce));
  mpi_printf(comm, "L/debye = %g\n", phys.L / (phys.vthe / phys.wpe));
  mpi_printf(comm, "dx/debye = %g\n",
             (phys.Lx / g.gdims[0]) / (phys.vthe / phys.wpe));
  mpi_printf(comm, "n0 = %g\n", phys.n0);
  mpi_printf(comm, "vthi/c = %g\n", phys.vthi / phys.c);
  mpi_printf(comm, "vthe/c = %g\n", phys.vthe / phys.c);
  mpi_printf(comm, "vdri/c = %g\n", phys.vdri / phys.c);
  mpi_printf(comm, "vdre/c = %g\n", phys.vdre / phys.c);
  mpi_printf(comm, "Open BC in x?   = %d\n", g.open_bc_x);
  mpi_printf(comm, "Driven BC in z? = %d\n", g.driven_bc_z);
}

// ======================================================================
// Diagnostics

class Diagnostics
{
public:
  Diagnostics(Grid_t& grid,
              OutputFields<MfieldsState, Mparticles, dim_xz>& outf,
              OutputParticles& outp, DiagEnergies& oute)
    : outf_{outf},
      outp_{outp},
      oute_{oute},
      outf_state_(),
      outf_hydro_(grid),
      mflds_acc_state_(grid, outf_state_.n_comps(), {1, 1, 1}),
      mflds_acc_hydro_(grid, outf_hydro_.n_comps(), {1, 1, 1})
  {
    io_pfd_.open("pfd");
    io_tfd_.open("tfd");
  }

  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    vpic_run_diagnostics(mprts, mflds);

    const auto& grid = mprts.grid();
    MPI_Comm comm = grid.comm();

    int timestep = grid.timestep();
    if (outf_.fields.pfield.out_interval > 0 &&
        timestep % outf_.fields.pfield.out_interval == 0) {
      mpi_printf(comm, "***** Writing PFD output\n");
      io_pfd_.begin_step(grid);

      {
        auto result = outf_state_(mflds);
        io_pfd_.write(adapt(evalMfields(result.mflds)), grid, result.name,
                      result.comp_names);
      }

      {
        auto result = outf_hydro_(mprts, *hydro, *interpolator);
        io_pfd_.write(adapt(result.mflds), grid, result.name,
                      result.comp_names);
      }

      io_pfd_.end_step();
    }

    if (outf_.fields.tfield.out_interval > 0) {
      auto result_state = outf_state_(mflds);

      for (int p = 0; p < mflds_acc_state_.n_patches(); p++) {
        auto flds_acc_state = make_Fields3d<dim_xyz>(mflds_acc_state_[p]);
        auto flds_state = make_Fields3d<dim_xyz>(result_state.mflds[p]);
        for (int m = 0; m < mflds_acc_state_.n_comps(); m++) {
          mflds_acc_state_.grid().Foreach_3d(0, 0, [&](int i, int j, int k) {
            flds_acc_state(m, i, j, k) += flds_state(m, i, j, k);
          });
        }
      }

      auto result_hydro = outf_hydro_(mprts, *hydro, *interpolator);
      for (int p = 0; p < mflds_acc_hydro_.n_patches(); p++) {
        auto flds_acc_hydro = make_Fields3d<dim_xyz>(mflds_acc_hydro_[p]);
        auto flds_hydro = make_Fields3d<dim_xyz>(result_hydro.mflds[p]);
        for (int m = 0; m < mflds_acc_hydro_.n_comps(); m++) {
          mflds_acc_hydro_.grid().Foreach_3d(0, 0, [&](int i, int j, int k) {
            flds_acc_hydro(m, i, j, k) += flds_hydro(m, i, j, k);
          });
        }
      }

      n_accum_++;

      if (timestep % outf_.fields.tfield.out_interval == 0) {
        mpi_printf(comm, "***** Writing TFD output\n");
        io_tfd_.begin_step(grid);
        mflds_acc_state_.storage() =
          (1. / n_accum_) * mflds_acc_state_.storage();
        io_tfd_.write(adapt(mflds_acc_state_), grid, result_state.name,
                      result_state.comp_names);
        mflds_acc_state_.storage().view() = 0.;

        mflds_acc_hydro_.storage() =
          (1. / n_accum_) * mflds_acc_hydro_.storage();
        io_tfd_.write(adapt(mflds_acc_hydro_), grid, result_hydro.name,
                      result_hydro.comp_names);
        mflds_acc_hydro_.storage().view() = 0.;

        io_tfd_.end_step();

        n_accum_ = 0;
      }
    }

    psc_stats_start(st_time_output);
    outp_(mprts);
    psc_stats_stop(st_time_output);

#ifndef VPIC
    oute_(mprts, mflds);
#endif
  }

private:
  WriterMRC io_pfd_;
  WriterMRC io_tfd_;
  OutputFields<MfieldsState, Mparticles, dim_xz>& outf_;
  OutputFieldsVpic<MfieldsState> outf_state_;
  OutputHydro outf_hydro_;
  OutputParticles& outp_;
  DiagEnergies& oute_;
  MfieldsSingle mflds_acc_state_;
  MfieldsSingle mflds_acc_hydro_;
  int n_accum_ = 0;
};

// ----------------------------------------------------------------------
// setup_particles
//
// set particles x^{n+1/2}, p^{n+1/2}

void setup_particles(Mparticles& mprts,
                     std::vector<uint>& nr_particles_by_patch, bool count_only)
{
  const auto& grid = mprts.grid();
  MPI_Comm comm = grid.comm();

  double cs = cos(g.theta), sn = sin(g.theta);
  double Ne_sheet = phys.Ne_sheet, vthe = phys.vthe, vthi = phys.vthi;
  int n_global_patches = phys.n_global_patches;
  double weight_s = phys.weight_s;
  double tanhf = phys.tanhf, L = phys.L;
  double gdre = phys.gdre, udre = phys.udre, gdri = phys.gdri, udri = phys.udri;
  double Ne_back = phys.Ne_back, vtheb = phys.vtheb, vthib = phys.vthib;
  double weight_b = phys.weight_b;

  if (count_only) {
    for (int p = 0; p < grid.n_patches(); p++) {
      nr_particles_by_patch[p] =
        2 * (Ne_sheet / n_global_patches + Ne_back / n_global_patches);
    }
    return;
  }

  // LOAD PARTICLES

  mpi_printf(comm, "Loading particles\n");

  // Do a fast load of the particles

  rngpool =
    RngPool_create(); // FIXME, should be part of ctor (of struct psc, really)

  int rank;
  MPI_Comm_rank(comm, &rank);
  RngPool_seed(rngpool, rank);
  Rng* rng = RngPool_get(rngpool, 0);

  assert(grid.n_patches() > 0);
  const Grid_t::Patch& patch = grid.patches[0];
  double xmin = patch.xb[0], xmax = patch.xe[0];
  double ymin = patch.xb[1], ymax = patch.xe[1];
  double zmin = patch.xb[2], zmax = patch.xe[2];

  // Load Harris population

  {
    auto inj = mprts.injector();
    auto injector = inj[0];

    mpi_printf(comm, "-> Main Harris Sheet\n");

    for (int64_t n = 0; n < Ne_sheet / n_global_patches; n++) {
      double x, y, z, ux, uy, uz, d0;

      do {
        z = L * atanh(Rng_uniform(rng, -1., 1.) * tanhf);
      } while (z <= zmin || z >= zmax);
      x = Rng_uniform(rng, xmin, xmax);
      y = Rng_uniform(rng, ymin, ymax);

      // inject_particles() will return an error for particles not on this
      // node and will not inject particle locally

      ux = Rng_normal(rng, 0, vthe);
      uy = Rng_normal(rng, 0, vthe);
      uz = Rng_normal(rng, 0, vthe);
      d0 = gdre * uy + sqrt(ux * ux + uy * uy + uz * uz + 1) * udre;
      uy = d0 * cs - ux * sn;
      ux = d0 * sn + ux * cs;

      injector.reweight(psc::particle::Inject{
        {x, y, z}, {ux, uy, uz}, weight_s, KIND_ELECTRON});

      ux = Rng_normal(rng, 0, vthi);
      uy = Rng_normal(rng, 0, vthi);
      uz = Rng_normal(rng, 0, vthi);
      d0 = gdri * uy + sqrt(ux * ux + uy * uy + uz * uz + 1) * udri;
      uy = d0 * cs - ux * sn;
      ux = d0 * sn + ux * cs;

      injector.reweight(
        psc::particle::Inject{{x, y, z}, {ux, uy, uz}, weight_s, KIND_ION});
    }

    mpi_printf(comm, "-> Background Population\n");

    for (int64_t n = 0; n < Ne_back / n_global_patches; n++) {
      Double3 pos{Rng_uniform(rng, xmin, xmax), Rng_uniform(rng, ymin, ymax),
                  Rng_uniform(rng, zmin, zmax)};

      Double3 u{Rng_normal(rng, 0, vtheb), Rng_normal(rng, 0, vtheb),
                Rng_normal(rng, 0, vtheb)};
      injector.reweight(psc::particle::Inject{pos, u, weight_b, KIND_ELECTRON});

      u = {Rng_normal(rng, 0, vthib), Rng_normal(rng, 0, vthib),
           Rng_normal(rng, 0, vthib)};
      injector.reweight(psc::particle::Inject{pos, u, weight_b, KIND_ION});
    }
  }
  mpi_printf(comm, "Finished loading particles\n");
}

// ======================================================================
// initializeParticles

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts)
{
  auto comm = grid_ptr->comm();

  mpi_printf(comm, "**** Partitioning...\n");
  std::vector<uint> n_prts_by_patch(mprts.n_patches());
  setup_particles(mprts, n_prts_by_patch, true);

  balance.initial(grid_ptr, n_prts_by_patch);
  mprts.reset(*grid_ptr);

  mpi_printf(comm, "**** Setting up particles...\n");
  mprts.reserve_all(n_prts_by_patch);
  setup_particles(mprts, n_prts_by_patch, false);
}

// ======================================================================
// initializeFields

void initializeFields(MfieldsState& mflds)
{
  double b0 = phys.b0, dbx = phys.dbx, dbz = phys.dbz;
  double L = phys.L, Lx = phys.Lx, Lz = phys.Lz, Lpert = phys.Lpert;
  double cs = cos(g.theta), sn = sin(g.theta);

  setupFields(mflds, [&](int m, double crd[3]) {
    double x = crd[0], z = crd[2];

    switch (m) {
      case HX:
        return cs * b0 * tanh(z / L) +
               dbx * cos(2. * M_PI * (x - .5 * Lx) / Lpert) *
                 sin(M_PI * z / Lz);

      case HY: return -sn * b0 * tanh(z / L) + b0 * g.bg;

      case HZ:
        return dbz * cos(M_PI * z / Lz) *
               sin(2.0 * M_PI * (x - 0.5 * Lx) / Lpert);

      case JYI: return 0.; // FIXME

      default: return 0.;
    }
  });
}

// ======================================================================
// run

void run()
{
  auto comm = MPI_COMM_WORLD;

  mpi_printf(comm, "*** Setting up simulation\n");

  setupHarrisParams();
  phys = globals_physics{g};

  psc_params.cfl = 0.99;
  psc_params.stats_every = 100;

  // ----------------------------------------------------------------------
  // Set up grid, state fields, particles

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  psc_params.nmax =
    int(g.taui / (phys.wci * grid.dt)); // number of steps from taui

  // --- create Simulation
#if 0
  // set high level VPIC simulation parameters
  // FIXME, will be unneeded eventually
  setParams(psc_params.nmax, psc_params.stats_every,
	    psc_params.stats_every / 2, psc_params.stats_every / 2,
	    psc_params.stats_every / 2);
#endif

  // --- setup field data structure
#ifdef VPIC
  // --- setup materials
  MaterialList material_list;
  setupMaterials(material_list);

  // FIXME, mv assert into MfieldsState ctor
  assert(!material_list.empty());

  double damp = 0.;
  MfieldsState mflds{grid, vgrid, material_list, damp};
  vpic_define_fields(grid);
#else
  MfieldsState mflds{grid};
#endif

  mpi_printf(comm, "*** Finalizing Field Advance\n");
#if 0
  assert(grid.nr_patches() > 0);
  Simulation_set_region_resistive_harris(sub->sim, &sub->prm, phys, psc_->patch[0].dx,
					 0., resistive);
#endif

  /// --- setup particle data structure
#ifdef VPIC
  Mparticles mprts{grid, vgrid};
  vpic_setup_species(mprts);
#else
  Mparticles mprts{grid};
#endif

  // -- Balance
  psc_params.balance_interval = 0;
  Balance balance{};

  // -- Sort
  // FIXME: the "vpic" sort actually keeps track of per-species sorting
  // intervals internally, so it needs to be called every step
#ifdef VPIC
  psc_params.sort_interval = 1;
#endif

  // -- Collision
  int collision_interval = 0;
  double collision_nu = .1; // FIXME, != 0 needed to avoid crash
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  Checks checks{grid, comm, checks_params};

  // -- Marder correction
  // FIXME, these are ignored for vpic (?)
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  // FIXME, how do we make sure that we don't forget to set this?
  // (maybe make it part of the Marder object, or provide a base class
  // interface define_marder() that takes the object and the interval
#ifdef VPIC
  psc_params.marder_interval = 1;
#else
  psc_params.marder_interval = 0;
#endif
#if 0
  // FIXME, marder "vpic" manages its own cleaning intervals
  psc_marder_set_param_int(psc_->marder, "every_step", 1);
  psc_marder_set_param_int(psc_->marder, "clean_div_e_interval", 50);
  psc_marder_set_param_int(psc_->marder, "clean_div_b_interval", 50);
  psc_marder_set_param_int(psc_->marder, "sync_shared_interval", 50);
  psc_marder_set_param_int(psc_->marder, "num_div_e_round", 2);
  psc_marder_set_param_int(psc_->marder, "num_div_b_round", 2);
#endif

  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // -- output fields
  OutputFieldsParams outf_params;
  double output_field_interval = .1;
  outf_params.fields.pfield.out_interval = 100;
  //    int((output_field_interval / (phys.wci * grid.dt)));
  outf_params.fields.tfield.out_interval = -1;
  //    int((output_field_interval / (phys.wci * grid.dt)));
  OutputFields<MfieldsState, Mparticles, dim_xz> outf{grid, outf_params};

  OutputParticlesParams outp_params{};
  outp_params.every_step =
    int((g.output_particle_interval / (phys.wci * grid.dt)));
  outp_params.data_dir = ".";
  outp_params.basename = "prt";
  outp_params.lo = {192, 0, 48};
  outp_params.hi = {320, 0, 80};
  OutputParticles outp{grid, outp_params};

  int oute_interval = 100;
  DiagEnergies oute{grid.comm(), oute_interval};

  Diagnostics diagnostics{grid, outf, outp, oute};

  // ---

  int interval = int(g.t_intervali / (phys.wci * grid.dt));
  vpic_create_diagnostics(interval);
  vpic_setup_diagnostics();
  setup_log(grid);

  mpi_printf(comm, "*** Finished with user-specified initialization ***\n");

  // ----------------------------------------------------------------------

  initializeParticles(balance, grid_ptr, mprts);
  initializeFields(mflds);

  // ----------------------------------------------------------------------
  // hand off to PscIntegrator to run the simulation

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

#if 0
  // FIXME, checkpoint reading should be moved to before the integrator
  if (!read_checkpoint_filename.empty()) {
    mpi_printf(MPI_COMM_WORLD, "**** Reading checkpoint...\n");
    psc.read_checkpoint(read_checkpoint_filename);
  }
#endif

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
