
#include <psc_config.h>

#include <psc.h>
#include <psc.hxx>

#ifdef USE_VPIC
#include "../libpsc/vpic/vpic_iface.h"
#endif

#include <psc_method.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_collision.h>
#include <psc_balance.h>
#include <psc_marder.h>
#include <psc_sort.h>
#include <psc_collision.h>
#include <psc_bnd_particles.h>
#include <psc_bnd.h>
#include <psc_bnd_fields.h>
#include <psc_checks.h>
#include <psc_output_fields_collection_private.h>
#include <psc_output_fields_private.h>
#include <psc_output_particles.h>

#include <push_particles.hxx>
#include <push_fields.hxx>
#include <sort.hxx>
#include <collision.hxx>
#include <bnd_particles.hxx>
#include <checks.hxx>
#include <marder.hxx>

#include "psc_fields_single.h"
#include "balance.hxx"
#include "fields3d.hxx"
#include "setup_particles.hxx"
#include "setup_fields.hxx"

#include <psc_particles_single.h>
#include <psc_particles_vpic.h>

#include "rngpool_iface.h"

#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct PscConfig
{
  using Mparticles_t = MparticlesVpic;
};

static RngPool *rngpool; // FIXME, should be member (of struct psc, really)

// FIXME, helper should go somewhere...

static inline double trunc_granular( double a, double b )
{
  return b * (int)(a/b);
}

// ----------------------------------------------------------------------
// courant length
//
// FIXME, the dt calculating should be consolidated with what regular PSC does

static inline double
courant_length(double length[3], int gdims[3])
{
  double inv_sum = 0.;
  for (int d = 0; d < 3; d++) {
    if (gdims[d] > 1) {
      inv_sum += sqr(gdims[d] / length[d]);
    }
  }
  return sqrt(1. / inv_sum);
}
 
// ======================================================================
// PscHarrisParams

struct PscHarrisParams
{
  double L_di;                     // Sheet thickness / ion inertial length
  double Ti_Te;                    // Ion temperature / electron temperature
  double nb_n0;                    // background plasma density
  double Tbe_Te;                   // Ratio of background T_e to Harris T_e
  double Tbi_Ti;                   // Ratio of background T_i to Harris T_i

  double bg;                       // Guide field
  double theta;
  
  double Lpert_Lx;                 // wavelength of perturbation in terms of Lx
  double dbz_b0;                   // perturbation in Bz relative to B0
  double nppc;                     // Average number of macro particle per cell per species
  bool open_bc_x;                  // Flag to signal we want to do open boundary condition in x
  bool driven_bc_z;                // Flag to signal we want to do driven boundary condition in z

  // FIXME, not really harris-specific
  double wpedt_max;

  double wpe_wce;                  // electron plasma freq / electron cyclotron freq
  double mi_me;                    // Ion mass / electron mass
  
  double Lx_di, Ly_di, Lz_di;      // Size of box in d_i

  int ion_sort_interval;
  int electron_sort_interval;

  double taui;                     // simulation wci's to run
  double t_intervali;              // output interval in terms of 1/wci
  double output_field_interval;    // field output interval in terms of 1/wci
  double output_particle_interval; // particle output interval in terms of 1/wci

  double overalloc;                // Overallocation factor (> 1) for particle arrays

  Int3 gdims;
  Int3 np;
};

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
  
  // general
  double dg;       // courant length
  double dt;       // timestep
  
  // calculated
  double b0;       // B0
  double n0;
  double v_A;
  double rhoi_L;
  double Lx, Ly, Lz; // size of box
  double L;        // Harris sheet thickness
  double Lpert;    // wavelength of perturbation
  double dbx;      // Perturbation in Bz relative to Bo (Only change here)
  double dbz;      // Set Bx perturbation so that div(B) = 0
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

  globals_physics(const PscHarrisParams& p)
  {
    assert(p.np[2] <= 2); // For load balance, keep "1" or "2" for Harris sheet

    // FIXME, the general normalization stuff should be shared somehow

    // use natural PIC units
    ec   = 1;         // Charge normalization
    me   = 1;         // Mass normalization
    c    = 1;         // Speed of light
    de   = 1;         // Length normalization (electron inertial length)
    eps0 = 1;         // Permittivity of space

    //derived quantities
    mi = me * p.mi_me;       // Ion mass

    double Te = me * sqr(c) / (2. * eps0 * sqr(p.wpe_wce) * (1. + p.Ti_Te)); // Electron temperature
    double Ti = Te * p.Ti_Te;           // Ion temperature
    vthe = sqrt(Te / me);               // Electron thermal velocity
    vthi = sqrt(Ti / mi);               // Ion thermal velocity
    vtheb = sqrt(p.Tbe_Te * Te / me);   // normalized background e thermal vel.
    vthib = sqrt(p.Tbi_Ti * Ti / mi);   // normalized background ion thermal vel.
    wci  = 1. / (p.mi_me * p.wpe_wce);  // Ion cyclotron frequency
    wce  = wci * p.mi_me;               // Electron cyclotron freqeuncy
    wpe  = wce * p.wpe_wce;             // electron plasma frequency
    wpi  = wpe / sqrt(p.mi_me);         // ion plasma frequency
    di   = c / wpi;                     // ion inertial length
    L    = p.L_di * di;                 // Harris sheet thickness
    rhoi_L = sqrt(p.Ti_Te / (1. + p.Ti_Te)) / p.L_di;
    v_A = (wci / wpi) / sqrt(p.nb_n0);  // based on nb

    Lx    = p.Lx_di * di;               // size of box in x dimension
    Ly    = p.Ly_di * di;               // size of box in y dimension
    Lz    = p.Lz_di * di;               // size of box in z dimension

    b0 = me*c*wce/ec;                   // Asymptotic magnetic field strength
    n0 = me*eps0*wpe*wpe/(ec*ec);       // Peak electron (ion) density
    vdri = 2*c*Ti/(ec*b0*L);            // Ion drift velocity
    vdre = -vdri/(p.Ti_Te);             // electron drift velocity

    n_global_patches = p.np[0] * p.np[1] * p.np[2];
    double Npe_sheet = 2*n0*Lx*Ly*L*tanh(0.5*Lz/L);         // N physical e's in sheet
    double Npe_back  = p.nb_n0 * n0 * Ly*Lz*Lx;             // N physical e's in backgrnd
    double Npe       = Npe_sheet + Npe_back;
    Ne         = p.nppc * p.gdims[0] * p.gdims[1] * p.gdims[2];   // total macro electrons in box
    Ne_sheet   = Ne * Npe_sheet / Npe;
    Ne_back    = Ne * Npe_back / Npe;
    Ne_sheet   = trunc_granular(Ne_sheet,n_global_patches); // Make it divisible by # subdomains
    Ne_back    = trunc_granular(Ne_back, n_global_patches); // Make it divisible by # subdomains
    Ne         = Ne_sheet + Ne_back;
    weight_s   = ec*Npe_sheet/Ne_sheet;                     // Charge per macro electron
    weight_b   = ec*Npe_back/Ne_back;                       // Charge per macro electron

    gdri  = 1./sqrt(1.-sqr(vdri)/sqr(c));  // gamma of ion drift frame
    gdre  = 1./sqrt(1.-sqr(vdre)/sqr(c));  // gamma of electron drift frame
    udri  = vdri * gdri;                   // 4-velocity of ion drift frame
    udre  = vdre * gdre;                   // 4-velocity of electron drift frame
    tanhf = tanh(0.5*Lz/L);
    Lpert = p.Lpert_Lx * Lx;               // wavelength of perturbation
    dbz   = p.dbz_b0 * b0;                 // Perturbation in Bz relative to Bo (Only change here)
    dbx   = -dbz * Lpert / (2.*Lz);        // Set Bx perturbation so that div(B) = 0
  }
};

// ======================================================================
// PscHarris

struct PscHarris : Psc<PscConfig>, PscHarrisParams
{
  // ----------------------------------------------------------------------
  // PscHarris ctor
  
  PscHarris(const PscHarrisParams& params, psc *psc)
    : Psc{psc},
      PscHarrisParams(params),
      phys_{params},
      mprts_{createMparticles(setup_grid())}
  {
    mprts__ = &mprts_;
    // partition and initial balancing
    std::vector<uint> n_prts_by_patch_old(psc_->n_patches());
    setup_particles(n_prts_by_patch_old, true);
    psc_balance_setup(psc_->balance);
    auto balance = PscBalanceBase{psc_->balance};
    auto n_prts_by_patch_new = balance.initial(psc_, n_prts_by_patch_old);
    mprts_.reset(psc_->grid());

    MPI_Comm comm = psc_comm(psc_);
    // create and set up base mflds
    psc_->flds = PscMfieldsCreate(comm, psc_->grid(),
				  psc_->n_state_fields, psc_->ibn, psc_->prm.fields_base).mflds();
    
    mprts_.reserve_all(n_prts_by_patch_new.data());
    setup_particles(n_prts_by_patch_new, false);

    // FIXME MfieldsSingle
    SetupFields<MfieldsSingle>::set(*PscMfieldsBase(psc->flds).sub(), [&](int m, double xx[3]) {
	return init_field(xx, m);
      });

    if (strcmp(psc_method_type(psc_->method), "vpic") != 0) {
    } else {
      Simulation_diagnostics_setup(sim_);
    }

    psc_setup_member_objs(psc_);
    initialize_stats();
    setup_log();
    
    if (output_field_interval > 0) {
      struct psc_output_fields *out;
      mrc_obj_for_each_child(out, psc_->output_fields_collection, struct psc_output_fields) {
	psc_output_fields_set_param_int(out, "pfield_step",
					(int) (output_field_interval / (phys_.wci*phys_.dt)));
      }
    }
  
    if (output_particle_interval > 0) {
      psc_output_particles_set_param_int(psc_->output_particles, "every_step",
					 (int) (output_particle_interval / (phys_.wci*phys_.dt)));
    }
  
    mpi_printf(comm, "*** Finished with user-specified initialization ***\n");
  }

  // ----------------------------------------------------------------------
  // PscHarris dtor
  
  ~PscHarris()
  {
    Simulation_delete(sim_);
  }

  // ----------------------------------------------------------------------
  // setup_grid

  const Grid_t& setup_grid()
  {
    auto domain = Grid_t::Domain{gdims,
				 {phys_.Lx, phys_.Ly, phys_.Lz},
				 {0., -.5 * phys_.Ly, -.5 * phys_.Lz}, np};

    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_REFLECTING },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_REFLECTING }};
 
    auto kinds = Grid_t::Kinds{{-phys_.ec, phys_.me, "e"},
			       { phys_.ec, phys_.mi, "i"}};

    // Determine the time step
    phys_.dg = courant_length(domain.length, domain.gdims);
    phys_.dt = psc_->prm.cfl * phys_.dg / phys_.c; // courant limited time step
    if (phys_.wpe * phys_.dt > wpedt_max) {
      phys_.dt = wpedt_max / phys_.wpe;  // override timestep if plasma frequency limited
    }
    psc_->dt = phys_.dt;
    
    psc_->prm.nmax = (int) (taui / (phys_.wci*phys_.dt)); // number of steps from taui

    psc_setup_coeff(psc_);
    psc_setup_domain(psc_, domain, grid_bc, kinds);
    if (strcmp(psc_method_type(psc_->method), "vpic") != 0) {
    } else {
      sim_ = Simulation_create();
      psc_method_set_param_ptr(psc_->method, "sim", sim_);
      // set high level VPIC simulation parameters
      // FIXME, will be unneeded eventually
      Simulation_set_params(sim_, psc_->prm.nmax, psc_->prm.stats_every,
			    psc_->prm.stats_every / 2, psc_->prm.stats_every / 2,
			    psc_->prm.stats_every / 2);
      setup_domain(domain);
      setup_fields();
      setup_species();

      int interval = (int) (t_intervali / (phys_.wci*phys_.dt));
      Simulation_diagnostics_init(sim_, interval);

      psc_->n_state_fields = VPIC_MFIELDS_N_COMP;
      psc_->ibn[0] = psc_->ibn[1] = psc_->ibn[2] = 1;
    }

    return psc_->grid();
  }

  // ----------------------------------------------------------------------
  // createMparticles

  Mparticles_t createMparticles(const Grid_t& grid)
  {
    mpi_printf(psc_comm(psc_), "**** Creating particle data structure...\n");
    return Mparticles_t{grid};
  }
  
  // ----------------------------------------------------------------------
  // setup_domain

  void setup_domain(const Grid_t::Domain& domain)
  {
    MPI_Comm comm = psc_comm(psc_);

    // Setup basic grid parameters
    double dx[3], xl[3], xh[3];
    for (int d = 0; d < 3; d++) {
      dx[d] = domain.length[d] / domain.gdims[d];
      xl[d] = domain.corner[d];
      xh[d] = xl[d] + domain.length[d];
    }
    Simulation_setup_grid(sim_, dx, phys_.dt, phys_.c, phys_.eps0);

    // Define the grid
    Simulation_define_periodic_grid(sim_, xl, xh, domain.gdims, domain.np);

    int p = 0;
    bool left = psc_at_boundary_lo(psc_, p, 0);
    bool right = psc_at_boundary_hi(psc_, p, 0);

    bool bottom = psc_at_boundary_lo(psc_, p, 2);
    bool top = psc_at_boundary_hi(psc_, p, 2);

    // ***** Set Field Boundary Conditions *****
    if (open_bc_x) {
      mpi_printf(comm, "Absorbing fields on X-boundaries\n");
      if (left ) Simulation_set_domain_field_bc(sim_, BOUNDARY(-1,0,0), BND_FLD_ABSORBING);
      if (right) Simulation_set_domain_field_bc(sim_, BOUNDARY( 1,0,0), BND_FLD_ABSORBING);
    }
  
    mpi_printf(comm, "Conducting fields on Z-boundaries\n");
    if (bottom) Simulation_set_domain_field_bc(sim_, BOUNDARY(0,0,-1), BND_FLD_CONDUCTING_WALL);
    if (top   ) Simulation_set_domain_field_bc(sim_, BOUNDARY(0,0, 1), BND_FLD_CONDUCTING_WALL);

    // ***** Set Particle Boundary Conditions *****
    if (driven_bc_z) {
      mpi_printf(comm, "Absorb particles on Z-boundaries\n");
      if (bottom) Simulation_set_domain_particle_bc(sim_, BOUNDARY(0,0,-1), BND_PRT_ABSORBING);
      if (top   ) Simulation_set_domain_particle_bc(sim_, BOUNDARY(0,0, 1), BND_PRT_ABSORBING);
    } else {
      mpi_printf(comm, "Reflect particles on Z-boundaries\n");
      if (bottom) Simulation_set_domain_particle_bc(sim_, BOUNDARY(0,0,-1), BND_PRT_REFLECTING);
      if (top   ) Simulation_set_domain_particle_bc(sim_, BOUNDARY(0,0, 1), BND_PRT_REFLECTING);
    }
    if (open_bc_x) {
      mpi_printf(comm, "Absorb particles on X-boundaries\n");
      if (left)  Simulation_set_domain_particle_bc(sim_, BOUNDARY(-1,0,0), BND_PRT_ABSORBING);
      if (right) Simulation_set_domain_particle_bc(sim_, BOUNDARY( 1,0,0), BND_PRT_ABSORBING);
    }
  }

  // ----------------------------------------------------------------------
  // setup_fields

  void setup_fields()
  {
    MPI_Comm comm = psc_comm(psc_);

    mpi_printf(comm, "Setting up materials.\n");

    Simulation_define_material(sim_, "vacuum", 1., 1., 0., 0.);
#if 0
    struct material *resistive =
      Simulation_define_material(sub->sim, "resistive", 1., 1., 1., 0.);
#endif
    Simulation_define_field_array(sim_, 0.);

    // Note: define_material defaults to isotropic materials with mu=1,sigma=0
    // Tensor electronic, magnetic and conductive materials are supported
    // though. See "shapes" for how to define them and assign them to regions.
    // Also, space is initially filled with the first material defined.

    //////////////////////////////////////////////////////////////////////////////
    // Finalize Field Advance

    mpi_printf(comm, "Finalizing Field Advance\n");

#if 0
    assert(psc_->nr_patches > 0);
    struct globals_physics *phys = &sub->phys;
    Simulation_set_region_resistive_harris(sub->sim, &sub->prm, phys, psc_->patch[0].dx,
					   0., resistive);
#endif
  }

  // ----------------------------------------------------------------------
  // setup_species

  void setup_species()
  {
    MPI_Comm comm = psc_comm(psc_);
    
    mpi_printf(comm, "Setting up species.\n");
    double nmax = overalloc * phys_.Ne / phys_.n_global_patches;
    double nmovers = .1 * nmax;
    double sort_method = 1;   // 0=in place and 1=out of place
    
    Simulation_define_species(sim_, "electron", -phys_.ec, phys_.me, nmax, nmovers,
			    electron_sort_interval, sort_method);
    Simulation_define_species(sim_, "ion", phys_.ec, phys_.mi, nmax, nmovers,
			      ion_sort_interval, sort_method);
  }

  // ----------------------------------------------------------------------
  // setup_log

  void setup_log()
  {
    MPI_Comm comm = psc_comm(psc_);

    mpi_printf(comm, "***********************************************\n");
    mpi_printf(comm, "* Topology: %d x %d x %d\n", np[0], np[1], np[2]);
    mpi_printf(comm, "tanhf    = %g\n", phys_.tanhf);
    mpi_printf(comm, "L_di     = %g\n", L_di);
    mpi_printf(comm, "rhoi/L   = %g\n", phys_.rhoi_L);
    mpi_printf(comm, "Ti/Te    = %g\n", Ti_Te) ;
    mpi_printf(comm, "nb/n0    = %g\n", nb_n0) ;
    mpi_printf(comm, "wpe/wce  = %g\n", wpe_wce);
    mpi_printf(comm, "mi/me    = %g\n", mi_me);
    mpi_printf(comm, "theta    = %g\n", theta);
    mpi_printf(comm, "Lpert/Lx = %g\n", Lpert_Lx);
    mpi_printf(comm, "dbz/b0   = %g\n", dbz_b0);
    mpi_printf(comm, "taui     = %g\n", taui);
    mpi_printf(comm, "t_intervali = %g\n", t_intervali);
    mpi_printf(comm, "num_step = %d\n", psc_->prm.nmax);
    mpi_printf(comm, "Lx/di = %g\n", phys_.Lx/phys_.di);
    mpi_printf(comm, "Lx/de = %g\n", phys_.Lx/phys_.de);
    mpi_printf(comm, "Ly/di = %g\n", phys_.Ly/phys_.di);
    mpi_printf(comm, "Ly/de = %g\n", phys_.Ly/phys_.de);
    mpi_printf(comm, "Lz/di = %g\n", phys_.Lz/phys_.di);
    mpi_printf(comm, "Lz/de = %g\n", phys_.Lz/phys_.de);
    mpi_printf(comm, "nx = %d\n", gdims[0]);
    mpi_printf(comm, "ny = %d\n", gdims[1]);
    mpi_printf(comm, "nz = %d\n", gdims[2]);
    mpi_printf(comm, "courant = %g\n", phys_.c*phys_.dt/phys_.dg);
    mpi_printf(comm, "n_global_patches = %d\n", phys_.n_global_patches);
    mpi_printf(comm, "nppc = %g\n", nppc);
    mpi_printf(comm, "b0 = %g\n", phys_.b0);
    mpi_printf(comm, "v_A (based on nb) = %g\n", phys_.v_A);
    mpi_printf(comm, "di = %g\n", phys_.di);
    mpi_printf(comm, "Ne = %g\n", phys_.Ne);
    mpi_printf(comm, "Ne_sheet = %g\n", phys_.Ne_sheet);
    mpi_printf(comm, "Ne_back = %g\n", phys_.Ne_back);
    mpi_printf(comm, "total # of particles = %g\n", 2*phys_.Ne);
    mpi_printf(comm, "dt*wpe = %g\n", phys_.wpe*phys_.dt);
    mpi_printf(comm, "dt*wce = %g\n", phys_.wce*phys_.dt);
    mpi_printf(comm, "dt*wci = %g\n", phys_.wci*phys_.dt);
    mpi_printf(comm, "dx/de = %g\n", phys_.Lx/(phys_.de*gdims[0]));
    mpi_printf(comm, "dy/de = %g\n", phys_.Ly/(phys_.de*gdims[1]));
    mpi_printf(comm, "dz/de = %g\n", phys_.Lz/(phys_.de*gdims[2]));
    mpi_printf(comm, "dx/rhoi = %g\n", (phys_.Lx/gdims[0])/(phys_.vthi/phys_.wci));
    mpi_printf(comm, "dx/rhoe = %g\n", (phys_.Lx/gdims[0])/(phys_.vthe/phys_.wce));
    mpi_printf(comm, "L/debye = %g\n", phys_.L/(phys_.vthe/phys_.wpe));
    mpi_printf(comm, "dx/debye = %g\n", (phys_.Lx/gdims[0])/(phys_.vthe/phys_.wpe));
    mpi_printf(comm, "n0 = %g\n", phys_.n0);
    mpi_printf(comm, "vthi/c = %g\n", phys_.vthi/phys_.c);
    mpi_printf(comm, "vthe/c = %g\n", phys_.vthe/phys_.c);
    mpi_printf(comm, "vdri/c = %g\n", phys_.vdri/phys_.c);
    mpi_printf(comm, "vdre/c = %g\n", phys_.vdre/phys_.c);
    mpi_printf(comm, "Open BC in x?   = %d\n", open_bc_x);
    mpi_printf(comm, "Driven BC in z? = %d\n", driven_bc_z);
  }

  // ----------------------------------------------------------------------
  // setup_particles
  //
  // set particles x^{n+1/2}, p^{n+1/2}

  void setup_particles(std::vector<uint>& nr_particles_by_patch, bool count_only)
  {
    MPI_Comm comm = psc_comm(psc_);

    double cs = cos(theta), sn = sin(theta);
    double Ne_sheet = phys_.Ne_sheet, vthe = phys_.vthe, vthi = phys_.vthi;
    int n_global_patches = phys_.n_global_patches;
    double weight_s = phys_.weight_s;
    double tanhf = phys_.tanhf, L = phys_.L;
    double gdre = phys_.gdre, udre = phys_.udre, gdri = phys_.gdri, udri = phys_.udri;
    double Ne_back = phys_.Ne_back, vtheb = phys_.vtheb, vthib = phys_.vthib;
    double weight_b = phys_.weight_b;

    if (count_only) {
      for (int p = 0; p < psc_->n_patches(); p++) {
	nr_particles_by_patch[p] = 2 * (Ne_sheet / n_global_patches + Ne_back / n_global_patches);
      }
      return;
    }
  
    // LOAD PARTICLES

    mpi_printf(comm, "Loading particles\n");

    // Do a fast load of the particles

    rngpool = RngPool_create(); // FIXME, should be part of ctor (of struct psc, really)

    int rank;
    MPI_Comm_rank(comm, &rank);
    RngPool_seed(rngpool, rank);
    Rng *rng = RngPool_get(rngpool, 0);

    assert(psc_->n_patches() > 0);
    const Grid_t& grid = psc_->grid();
    const Grid_t::Patch& patch = grid.patches[0];
    double xmin = patch.xb[0], xmax = patch.xe[0];
    double ymin = patch.xb[1], ymax = patch.xe[1];
    double zmin = patch.xb[2], zmax = patch.xe[2];

    // Load Harris population

    mpi_printf(comm, "-> Main Harris Sheet\n");

    for (int64_t n = 0; n < Ne_sheet / n_global_patches; n++) {
      double x, y, z, ux, uy, uz, d0;

      particle_inject_t prt;

      do {
	z = L*atanh(Rng_uniform(rng, -1., 1.)*tanhf);
      } while (z <= zmin || z >= zmax);
      x = Rng_uniform(rng, xmin, xmax);
      y = Rng_uniform(rng, ymin, ymax);

      // inject_particles() will return an error for particles no on this
      // node and will not inject particle locally

      ux = Rng_normal(rng, 0, vthe);
      uy = Rng_normal(rng, 0, vthe);
      uz = Rng_normal(rng, 0, vthe);
      d0 = gdre*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udre;
      uy = d0*cs - ux*sn;
      ux = d0*sn + ux*cs;

      prt.x[0] = x; prt.x[1] = y; prt.x[2] = z;
      prt.u[0] = ux; prt.u[1] = uy; prt.u[2] = uz;
      prt.w = weight_s;
      prt.kind = KIND_ELECTRON;
      mprts_.inject_reweight(0, prt);

      ux = Rng_normal(rng, 0, vthi);
      uy = Rng_normal(rng, 0, vthi);
      uz = Rng_normal(rng, 0, vthi);
      d0 = gdri*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udri;
      uy = d0*cs - ux*sn;
      ux = d0*sn + ux*cs;

      prt.u[0] = ux; prt.u[1] = uy; prt.u[2] = uz;
      prt.kind = KIND_ION;
      mprts_.inject_reweight(0, prt);
    }

    mpi_printf(comm, "-> Background Population\n");

    for (int64_t n = 0; n < Ne_back / n_global_patches; n++) {
      particle_inject_t prt;
      double x = Rng_uniform(rng, xmin, xmax);
      double y = Rng_uniform(rng, ymin, ymax);
      double z = Rng_uniform(rng, zmin, zmax);

      prt.x[0] = x; prt.x[1] = y; prt.x[2] = z;
      prt.u[0] = Rng_normal(rng, 0, vtheb);
      prt.u[1] = Rng_normal(rng, 0, vtheb);
      prt.u[2] = Rng_normal(rng, 0, vtheb);
      prt.w = weight_b;
      prt.kind = KIND_ELECTRON;
      mprts_.inject_reweight(0, prt);
    
      prt.u[0] = Rng_normal(rng, 0, vthib);
      prt.u[1] = Rng_normal(rng, 0, vthib);
      prt.u[2] = Rng_normal(rng, 0, vthib);
      prt.kind = KIND_ION;
      mprts_.inject_reweight(0, prt);
    }

    mpi_printf(comm, "Finished loading particles\n");
  }

  // ----------------------------------------------------------------------
  // init_field

  double init_field(double crd[3], int m)
  {
    double b0 = phys_.b0, dbx = phys_.dbx, dbz = phys_.dbz;
    double L = phys_.L, Lx = phys_.Lx, Lz = phys_.Lz, Lpert = phys_.Lpert;
    double x = crd[0], z = crd[2];
    
    double cs = cos(theta), sn = sin(theta);
    
    switch (m) {
    case HX:
      return cs*b0*tanh(z/L)+dbx*cos(2.*M_PI*(x-.5*Lx)/Lpert)*sin(M_PI*z/Lz);
      
    case HY:
    return -sn*b0*tanh(z/L) + b0*bg;
    
    case HZ:
      return dbz*cos(M_PI*z/Lz)*sin(2.0*M_PI*(x-0.5*Lx)/Lpert);
      
    case JYI:
      return 0.; // FIXME
      
    default:
      return 0.;
    }
  }

  // ----------------------------------------------------------------------
  // step

  void step()
  {
    if (!pr_time_step_no_comm) {
      pr_time_step_no_comm = prof_register("time step w/o comm", 1., 0, 0);
    }

#if 0
    mpi_printf(psc_comm(psc), "**** Step %d / %d, Time %g\n", psc->timestep + 1,
	       psc->prm.nmax, psc->timestep * psc->dt);
#endif

    // x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}

    PscMfieldsBase mflds(psc_->flds);
    MfieldsBase& mflds_{*mflds.sub()};
    PscPushParticlesBase pushp(psc_->push_particles);
    PscPushFieldsBase pushf(psc_->push_fields);
    PscSortBase sort(psc_->sort);
    PscCollisionBase collision(psc_->collision);
    PscBndParticlesBase bndp(psc_->bnd_particles);

    auto balance = PscBalanceBase{psc_->balance};
    balance(psc_, mprts_);

    prof_start(pr_time_step_no_comm);
    prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

    sort(mprts_);
    collision(mprts_);
  
    //psc_bnd_particles_open_calc_moments(psc_->bnd_particles, psc_->particles);

    PscChecksBase{psc_->checks}.continuity_before_particle_push(psc_, mprts_);

    // particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    pushp(mprts_, mflds_);
    // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}
    
    // field propagation B^{n+1/2} -> B^{n+1}
    pushf.advance_H(mflds_, .5);
    // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    bndp(mprts_);
  
    // field propagation E^{n+1/2} -> E^{n+3/2}
    pushf.advance_b2(mflds_);
    // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}

    // field propagation B^{n+1} -> B^{n+3/2}
    pushf.advance_a(mflds_);
    // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}

    PscChecksBase{psc_->checks}.continuity_after_particle_push(psc_, mprts_);

    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not it's natural time,
    // but div B should be == 0 at any time...)
    PscMarderBase{psc_->marder}(mflds, mprts_);
    
    PscChecksBase{psc_->checks}.gauss(psc_, mprts_);

    psc_push_particles_prep(psc_->push_particles, mprts_, mflds_);
  }

protected:
  globals_physics phys_;
  Mparticles_t mprts_;
  Simulation* sim_;
};

// ======================================================================
// PscHarrisBuilder

struct PscHarrisBuilder
{
  PscHarris* makePsc();
};

// ----------------------------------------------------------------------
// PscHarrisBuilder::makePsc

PscHarris* PscHarrisBuilder::makePsc()
{
  auto psc_ = psc_create(MPI_COMM_WORLD);
  MPI_Comm comm = psc_comm(psc_);
  
  mpi_printf(comm, "*** Setting up...\n");

  PscHarrisParams params;

  params.wpedt_max = .36;
  params.wpe_wce = 2.;
  params.mi_me = 25.;

  params.Lx_di = 40.;
  params.Ly_di = 1.;
  params.Lz_di = 10.;

  params.electron_sort_interval = 25;
  params.ion_sort_interval = 25;
  
  params.taui = 40.;
  params.t_intervali = 1.;
  params.output_field_interval = 1.;
  params.output_particle_interval = 10.;

  params.overalloc = 2.;

  params.gdims = {512, 1, 128};
  params.np = { 4, 1, 1 };
  
  params.L_di = .5;
  params.Ti_Te = 5.;
  params.nb_n0 = .05;
  params.Tbe_Te = .333;
  params.Tbi_Ti = .333;

  params.bg = 0.;
  params.theta = 0.;

  params.Lpert_Lx = 1.;
  params.dbz_b0 = .03;
  params.nppc = 10;
  params.open_bc_x = false;
  params.driven_bc_z = false;
  
  psc_default_dimensionless(psc_);

  psc_->prm.nicell = 1;
  psc_->prm.cfl = 0.99;

  psc_->prm.stats_every = 100;

  psc_method_set_type(psc_->method, "vpic");

  psc_sort_set_type(psc_->sort, "vpic");
  // FIXME: the "vpic" sort actually keeps track of per-species sorting intervals
  // internally
  psc_sort_set_param_int(psc_->sort, "every", 1);

  psc_collision_set_type(psc_->collision, "vpic");

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc_->push_particles, "vpic");
  psc_bnd_particles_set_type(psc_->bnd_particles, "vpic");

  psc_push_fields_set_type(psc_->push_fields, "vpic");
  psc_bnd_set_type(psc_->bnd, "vpic");
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc_->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "vpic");

  psc_marder_set_type(psc_->marder, "vpic");
  // FIXME, marder "vpic" manages its own cleaning intervals
  psc_marder_set_param_int(psc_->marder, "every_step", 1);
  psc_marder_set_param_int(psc_->marder, "clean_div_e_interval", 50);
  psc_marder_set_param_int(psc_->marder, "clean_div_b_interval", 50);
  psc_marder_set_param_int(psc_->marder, "sync_shared_interval", 50);
  psc_marder_set_param_int(psc_->marder, "num_div_e_round", 2);
  psc_marder_set_param_int(psc_->marder, "num_div_b_round", 2);

  psc_set_from_options(psc_);

  return new PscHarris{params, psc_};
}

// ======================================================================
// main

int
main(int argc, char **argv)
{
#ifdef USE_VPIC
  vpic_base_init(&argc, &argv);
#else
  MPI_Init(&argc, &argv);
#endif
  libmrc_params_init(argc, argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);

  auto builder = PscHarrisBuilder{};
  auto psc = builder.makePsc();

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
