
#include <psc_config.h>

#include <psc.h>

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
};

// ======================================================================
// PscHarris

struct PscHarris : PscHarrisParams
{
  using Mparticles_t = MparticlesSingle; // FIXME!!!
  
  PscHarris(const PscHarrisParams& params, psc *psc)
    : PscHarrisParams(params),
      psc_{psc}
  {
    psc_harris* sub = psc_harris(psc_);
    globals_physics* phys = &sub->phys;
    MPI_Comm comm = psc_comm(psc_);

    setup_ic();

    // Determine the time step
    phys->dg = courant_length(psc_->domain_.length, psc_->domain_.gdims);
    phys->dt = psc_->prm.cfl * phys->dg / phys->c; // courant limited time step
    if (phys->wpe * phys->dt > wpedt_max) {
      phys->dt = wpedt_max / phys->wpe;  // override timestep if plasma frequency limited
    }
    psc_->dt = phys->dt;
    
    psc_->prm.nmax = (int) (taui / (phys->wci*phys->dt)); // number of steps from taui

    bool split;
    psc_method_get_param_bool(psc_->method, "split", &split);
    if (strcmp(psc_method_type(psc_->method), "vpic") != 0 || !split) {
      psc_method_do_setup(psc_->method, psc_);
    } else {
      sub->sim = Simulation_create();
      psc_method_set_param_ptr(psc_->method, "sim", sub->sim);
      // set high level VPIC simulation parameters
      // FIXME, will be unneeded eventually
      Simulation_set_params(sub->sim, psc_->prm.nmax, psc_->prm.stats_every,
			    psc_->prm.stats_every / 2, psc_->prm.stats_every / 2,
			    psc_->prm.stats_every / 2);
      setup_domain();
      setup_fields();
      setup_species();

      int interval = (int) (t_intervali / (phys->wci*phys->dt));
      Simulation_diagnostics_init(sub->sim, interval);

      psc_->n_state_fields = VPIC_MFIELDS_N_COMP;
      psc_->ibn[0] = psc_->ibn[1] = psc_->ibn[2] = 1;
    }

    // partition and initial balancing
    auto n_prts_by_patch_old = psc_method_setup_partition(psc_->method, psc_);
    psc_balance_setup(psc_->balance);
    auto balance = PscBalanceBase{psc_->balance};
    auto n_prts_by_patch_new = balance.initial(psc_, n_prts_by_patch_old);

    // create and initialize base particle data structure x^{n+1/2}, p^{n+1/2}
    psc_->particles = PscMparticlesCreate(comm, psc_->grid(),
					  psc_->prm.particles_base).mprts();
    
    // create and set up base mflds
    psc_->flds = PscMfieldsCreate(comm, psc_->grid(),
				  psc_->n_state_fields, psc_->ibn, psc_->prm.fields_base).mflds();
    
    auto mprts_base = PscMparticlesBase{psc_->particles};
    mprts_base->reserve_all(n_prts_by_patch_new.data());
    psc_ops(psc)->setup_particles(psc, n_prts_by_patch_new, false);

    // FIXME MfieldsSingle
    SetupFields<MfieldsSingle>::set(*PscMfieldsBase(psc->flds).sub(), [&](int m, double xx[3]) {
	return init_field(xx, m);
      });

    if (strcmp(psc_method_type(psc_->method), "vpic") != 0 || !split) {
    } else {
      Simulation_diagnostics_setup(sub->sim);
    }

    psc_setup_member_objs(psc_);
    setup_log();
    
    if (output_field_interval > 0) {
      struct psc_output_fields *out;
      mrc_obj_for_each_child(out, psc_->output_fields_collection, struct psc_output_fields) {
	psc_output_fields_set_param_int(out, "pfield_step",
					(int) (output_field_interval / (phys->wci*phys->dt)));
      }
    }
  
    if (output_particle_interval > 0) {
      psc_output_particles_set_param_int(psc_->output_particles, "every_step",
					 (int) (output_particle_interval / (phys->wci*phys->dt)));
    }
  
    mpi_printf(comm, "*** Finished with user-specified initialization ***\n");
  }

  // ----------------------------------------------------------------------
  // setup_ic

  void setup_ic()
  {
    struct psc_harris *sub = psc_harris(psc_);
    struct globals_physics *phys = &sub->phys;
    struct vpic_harris_params *prm = &sub->prm;

    Int3 gdims = {512, 1, 128};
    Int3 np = {4, 1, 1};
    sub->n_global_patches = np[0] * np[1] * np[2];
  
    assert(np[2] <= 2); // For load balance, keep "1" or "2" for Harris sheet

    // FIXME, the general normalization stuff should be shared somehow

    // use natural PIC units
    phys->ec   = 1;         // Charge normalization
    phys->me   = 1;         // Mass normalization
    phys->c    = 1;         // Speed of light
    phys->de   = 1;         // Length normalization (electron inertial length)
    phys->eps0 = 1;         // Permittivity of space

    double c = phys->c;
    //derived quantities
    phys->mi = phys->me*mi_me;       // Ion mass
    double Te = phys->me*c*c/(2*phys->eps0*sqr(wpe_wce)*(1.+Ti_Te)); // Electron temperature
    double Ti = Te*Ti_Te;       // Ion temperature
    phys->vthe = sqrt(Te/phys->me);         // Electron thermal velocity
    phys->vthi = sqrt(Ti/phys->mi);         // Ion thermal velocity
    phys->vtheb = sqrt(Tbe_Te*Te/phys->me);  // normalized background e thermal vel.
    phys->vthib = sqrt(Tbi_Ti*Ti/phys->mi);  // normalized background ion thermal vel.
    phys->wci  = 1.0/(mi_me*wpe_wce);  // Ion cyclotron frequency
    phys->wce  = phys->wci*mi_me;            // Electron cyclotron freqeuncy
    phys->wpe  = phys->wce*wpe_wce;          // electron plasma frequency
    phys->wpi  = phys->wpe/sqrt(mi_me);      // ion plasma frequency
    phys->di   = c/phys->wpi;                      // ion inertial length
    phys->L    = L_di*phys->di;              // Harris sheet thickness
    phys->rhoi_L = sqrt(Ti_Te/(1.+Ti_Te))/L_di;
    phys->v_A = (phys->wci/phys->wpi)/sqrt(nb_n0); // based on nb

    phys->Lx    = Lx_di*phys->di; // size of box in x dimension
    phys->Ly    = Ly_di*phys->di; // size of box in y dimension
    phys->Lz    = Lz_di*phys->di; // size of box in z dimension

    phys->b0 = phys->me*c*phys->wce/phys->ec; // Asymptotic magnetic field strength
    phys->n0 = phys->me*phys->eps0*phys->wpe*phys->wpe/(phys->ec*phys->ec);  // Peak electron (ion) density
    phys->vdri = 2*c*Ti/(phys->ec*phys->b0*phys->L);   // Ion drift velocity
    phys->vdre = -phys->vdri/(Ti_Te);      // electron drift velocity

    double Lx = phys->Lx, Ly = phys->Ly, Lz = phys->Lz, L = phys->L;
    double Npe_sheet = 2*phys->n0*Lx*Ly*L*tanh(0.5*Lz/L); // N physical e's in sheet
    double Npe_back  = nb_n0*phys->n0 * Ly*Lz*Lx;          // N physical e's in backgrnd
    double Npe       = Npe_sheet + Npe_back;
    phys->Ne         = prm->nppc * gdims[0] * gdims[1] * gdims[2];  // total macro electrons in box
    phys->Ne_sheet   = phys->Ne*Npe_sheet/Npe;
    phys->Ne_back    = phys->Ne*Npe_back/Npe;
    phys->Ne_sheet   = trunc_granular(phys->Ne_sheet,sub->n_global_patches); // Make it divisible by # subdomains
    phys->Ne_back    = trunc_granular(phys->Ne_back, sub->n_global_patches); // Make it divisible by # subdomains
    phys->Ne         = phys->Ne_sheet + phys->Ne_back;
    phys->weight_s   = phys->ec*Npe_sheet/phys->Ne_sheet;  // Charge per macro electron
    phys->weight_b   = phys->ec*Npe_back/phys->Ne_back;  // Charge per macro electron

    phys->gdri  = 1/sqrt(1-phys->vdri*phys->vdri/(c*c));  // gamma of ion drift frame
    phys->gdre  = 1/sqrt(1-phys->vdre*phys->vdre/(c*c)); // gamma of electron drift frame
    phys->udri  = phys->vdri*phys->gdri;                 // 4-velocity of ion drift frame
    phys->udre  = phys->vdre*phys->gdre;                 // 4-velocity of electron drift frame
    phys->tanhf = tanh(0.5*Lz/L);
    phys->Lpert = prm->Lpert_Lx*Lx; // wavelength of perturbation
    phys->dbz   = prm->dbz_b0*phys->b0; // Perturbation in Bz relative to Bo (Only change here)
    phys->dbx   = -phys->dbz*phys->Lpert/(2.0*Lz); // Set Bx perturbation so that div(B) = 0

    psc_->domain_ = Grid_t::Domain{gdims,
				  {phys->Lx, phys->Ly, phys->Lz},
				  {0., -.5 * phys->Ly, -.5 * phys->Lz}, np};
  }

  // ----------------------------------------------------------------------
  // setup_domain

  void setup_domain()
  {
    struct psc_harris *sub = psc_harris(psc_);
    struct globals_physics *phys = &sub->phys;
    MPI_Comm comm = psc_comm(psc_);

    psc_setup_coeff(psc_); // FIXME -- in the middle of things here, will be done again later
    psc_setup_domain(psc_, psc_->domain_, psc_->bc_, psc_->kinds_);
  
    // Setup basic grid parameters
    double dx[3], xl[3], xh[3];
    for (int d = 0; d < 3; d++) {
      dx[d] = psc_->domain_.length[d] / psc_->domain_.gdims[d];
      xl[d] = psc_->domain_.corner[d];
      xh[d] = xl[d] + psc_->domain_.length[d];
    }
    Simulation_setup_grid(sub->sim, dx, phys->dt, phys->c, phys->eps0);

    // Define the grid
    Simulation_define_periodic_grid(sub->sim, xl, xh, psc_->domain_.gdims, psc_->domain_.np);

    int p = 0;
    bool left = psc_at_boundary_lo(psc_, p, 0);
    bool right = psc_at_boundary_hi(psc_, p, 0);

    bool bottom = psc_at_boundary_lo(psc_, p, 2);
    bool top = psc_at_boundary_hi(psc_, p, 2);

    // ***** Set Field Boundary Conditions *****
    if (sub->prm.open_bc_x) {
      mpi_printf(comm, "Absorbing fields on X-boundaries\n");
      if (left ) Simulation_set_domain_field_bc(sub->sim, BOUNDARY(-1,0,0), BND_FLD_ABSORBING);
      if (right) Simulation_set_domain_field_bc(sub->sim, BOUNDARY( 1,0,0), BND_FLD_ABSORBING);
    }
  
    mpi_printf(comm, "Conducting fields on Z-boundaries\n");
    if (bottom) Simulation_set_domain_field_bc(sub->sim, BOUNDARY(0,0,-1), BND_FLD_CONDUCTING_WALL);
    if (top   ) Simulation_set_domain_field_bc(sub->sim, BOUNDARY(0,0, 1), BND_FLD_CONDUCTING_WALL);

    // ***** Set Particle Boundary Conditions *****
    if (sub->prm.driven_bc_z) {
      mpi_printf(comm, "Absorb particles on Z-boundaries\n");
      if (bottom) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0,-1), BND_PRT_ABSORBING);
      if (top   ) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0, 1), BND_PRT_ABSORBING);
    } else {
      mpi_printf(comm, "Reflect particles on Z-boundaries\n");
      if (bottom) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0,-1), BND_PRT_REFLECTING);
      if (top   ) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0, 1), BND_PRT_REFLECTING);
    }
    if (sub->prm.open_bc_x) {
      mpi_printf(comm, "Absorb particles on X-boundaries\n");
      if (left)  Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(-1,0,0), BND_PRT_ABSORBING);
      if (right) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY( 1,0,0), BND_PRT_ABSORBING);
    }
  }

  // ----------------------------------------------------------------------
  // setup_fields

  void setup_fields()
  {
    struct psc_harris *sub = psc_harris(psc_);
    MPI_Comm comm = psc_comm(psc_);

    mpi_printf(comm, "Setting up materials.\n");

    Simulation_define_material(sub->sim, "vacuum", 1., 1., 0., 0.);
#if 0
    struct material *resistive =
      Simulation_define_material(sub->sim, "resistive", 1., 1., 1., 0.);
#endif
    Simulation_define_field_array(sub->sim, 0.);

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
    psc_harris* sub = psc_harris(psc_);
    globals_physics* phys = &sub->phys;
    MPI_Comm comm = psc_comm(psc_);
    
    mpi_printf(comm, "Setting up species.\n");
    double nmax = sub->prm.overalloc * phys->Ne / sub->n_global_patches;
    double nmovers = .1 * nmax;
    double sort_method = 1;   // 0=in place and 1=out of place
    
    psc_set_kinds(psc_, {{-phys->ec, phys->me, "e"}, {phys->ec, phys->mi, "i"}});
    
    Simulation_define_species(sub->sim, "electron", -phys->ec, phys->me, nmax, nmovers,
			    electron_sort_interval, sort_method);
    Simulation_define_species(sub->sim, "ion", phys->ec, phys->mi, nmax, nmovers,
			      ion_sort_interval, sort_method);
  }

  // ----------------------------------------------------------------------
  // setup_log

  void setup_log()
  {
    struct psc_harris *sub = psc_harris(psc_);
    struct globals_physics *phys = &sub->phys;
    MPI_Comm comm = psc_comm(psc_);

    mpi_printf(comm, "***********************************************\n");
    mpi_printf(comm, "* Topology: %d x %d x %d\n",
	       psc_->domain_.np[0], psc_->domain_.np[1], psc_->domain_.np[2]);
    mpi_printf(comm, "tanhf    = %g\n", phys->tanhf);
    mpi_printf(comm, "L_di     = %g\n", L_di);
    mpi_printf(comm, "rhoi/L   = %g\n", phys->rhoi_L);
    mpi_printf(comm, "Ti/Te    = %g\n", Ti_Te) ;
    mpi_printf(comm, "nb/n0    = %g\n", nb_n0) ;
    mpi_printf(comm, "wpe/wce  = %g\n", wpe_wce);
    mpi_printf(comm, "mi/me    = %g\n", mi_me);
    mpi_printf(comm, "theta    = %g\n", sub->prm.theta);
    mpi_printf(comm, "Lpert/Lx = %g\n", sub->prm.Lpert_Lx);
    mpi_printf(comm, "dbz/b0   = %g\n", sub->prm.dbz_b0);
    mpi_printf(comm, "taui     = %g\n", taui);
    mpi_printf(comm, "t_intervali = %g\n", t_intervali);
    mpi_printf(comm, "num_step = %d\n", psc_->prm.nmax);
    mpi_printf(comm, "Lx/di = %g\n", phys->Lx/phys->di);
    mpi_printf(comm, "Lx/de = %g\n", phys->Lx/phys->de);
    mpi_printf(comm, "Ly/di = %g\n", phys->Ly/phys->di);
    mpi_printf(comm, "Ly/de = %g\n", phys->Ly/phys->de);
    mpi_printf(comm, "Lz/di = %g\n", phys->Lz/phys->di);
    mpi_printf(comm, "Lz/de = %g\n", phys->Lz/phys->de);
    mpi_printf(comm, "nx = %d\n", psc_->domain_.gdims[0]);
    mpi_printf(comm, "ny = %d\n", psc_->domain_.gdims[1]);
    mpi_printf(comm, "nz = %d\n", psc_->domain_.gdims[2]);
    mpi_printf(comm, "courant = %g\n", phys->c*phys->dt/phys->dg);
    mpi_printf(comm, "n_global_patches = %d\n", sub->n_global_patches);
    mpi_printf(comm, "nppc = %g\n", sub->prm.nppc);
    mpi_printf(comm, "b0 = %g\n", phys->b0);
    mpi_printf(comm, "v_A (based on nb) = %g\n", phys->v_A);
    mpi_printf(comm, "di = %g\n", phys->di);
    mpi_printf(comm, "Ne = %g\n", phys->Ne);
    mpi_printf(comm, "Ne_sheet = %g\n", phys->Ne_sheet);
    mpi_printf(comm, "Ne_back = %g\n", phys->Ne_back);
    mpi_printf(comm, "total # of particles = %g\n", 2*phys->Ne);
    mpi_printf(comm, "dt*wpe = %g\n", phys->wpe*phys->dt);
    mpi_printf(comm, "dt*wce = %g\n", phys->wce*phys->dt);
    mpi_printf(comm, "dt*wci = %g\n", phys->wci*phys->dt);
    mpi_printf(comm, "dx/de = %g\n", phys->Lx/(phys->de*psc_->domain_.gdims[0]));
    mpi_printf(comm, "dy/de = %g\n", phys->Ly/(phys->de*psc_->domain_.gdims[1]));
    mpi_printf(comm, "dz/de = %g\n", phys->Lz/(phys->de*psc_->domain_.gdims[2]));
    mpi_printf(comm, "dx/rhoi = %g\n", (phys->Lx/psc_->domain_.gdims[0])/(phys->vthi/phys->wci));
    mpi_printf(comm, "dx/rhoe = %g\n", (phys->Lx/psc_->domain_.gdims[0])/(phys->vthe/phys->wce));
    mpi_printf(comm, "L/debye = %g\n", phys->L/(phys->vthe/phys->wpe));
    mpi_printf(comm, "dx/debye = %g\n", (phys->Lx/psc_->domain_.gdims[0])/(phys->vthe/phys->wpe));
    mpi_printf(comm, "n0 = %g\n", phys->n0);
    mpi_printf(comm, "vthi/c = %g\n", phys->vthi/phys->c);
    mpi_printf(comm, "vthe/c = %g\n", phys->vthe/phys->c);
    mpi_printf(comm, "vdri/c = %g\n", phys->vdri/phys->c);
    mpi_printf(comm, "vdre/c = %g\n", phys->vdre/phys->c);
    mpi_printf(comm, "Open BC in x?   = %d\n", sub->prm.open_bc_x);
    mpi_printf(comm, "Driven BC in z? = %d\n", sub->prm.driven_bc_z);
  }

  // ----------------------------------------------------------------------
  // init_field

  double init_field(double crd[3], int m)
  {
    struct psc_harris *sub = psc_harris(psc_);
    struct globals_physics *phys = &sub->phys;
    double theta = sub->prm.theta;
    double b0 = phys->b0, bg = sub->prm.bg, dbx = phys->dbx, dbz = phys->dbz;
    double L = phys->L, Lx = phys->Lx, Lz = phys->Lz, Lpert = phys->Lpert;
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
  // integrate

  void integrate()
  {
    psc_integrate(psc_);
  }

private:
  psc* psc_;
};

using Mparticles_t = MparticlesSingle;

static RngPool *rngpool; // FIXME, should be member (of struct psc, really)

// ----------------------------------------------------------------------

#define VAR(x) (void *)offsetof(struct psc_harris, x)
static struct param psc_harris_descr[] = {
  { "nppc"                   , VAR(prm.nppc)                 , PARAM_DOUBLE(100.),
    .help = "average number of macro particle per cell per species" },

  { "bg"                    , VAR(prm.bg)                    , PARAM_DOUBLE(0.),
    .help = "Guide field" },
  { "theta"                 , VAR(prm.theta)                 , PARAM_DOUBLE(0.),
    .help = "Theta" },
  { "Lpert_Lx"              , VAR(prm.Lpert_Lx)              , PARAM_DOUBLE(1.),
    .help = "wavelength of perturbation in terms of Lx" },
  { "dbz_b0"                , VAR(prm.dbz_b0)                , PARAM_DOUBLE(.03),
    .help = "perturbation in Bz relative to B0" },
  { "open_bc_x"             , VAR(prm.open_bc_x)             , PARAM_BOOL(false),
    .help = "use open b.c. at x boundaries" },
  { "driven_bc_z"           , VAR(prm.driven_bc_z)           , PARAM_BOOL(false),
    .help = "use driven b.c. at z boundaries" },

  { "overalloc"             , VAR(prm.overalloc)             , PARAM_DOUBLE(2.),
    .help = "over-allocate particle arrays by this factor" },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_harris_setup_particles
//
// set particles x^{n+1/2}, p^{n+1/2}

static void
psc_harris_setup_particles(struct psc *psc, std::vector<uint>& nr_particles_by_patch, bool count_only)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  MPI_Comm comm = psc_comm(psc);

  double cs = cos(sub->prm.theta), sn = sin(sub->prm.theta);
  double Ne_sheet = phys->Ne_sheet, vthe = phys->vthe, vthi = phys->vthi;
  double weight_s = phys->weight_s;
  double tanhf = phys->tanhf, L = phys->L;
  double gdre = phys->gdre, udre = phys->udre, gdri = phys->gdri, udri = phys->udri;
  double Ne_back = phys->Ne_back, vtheb = phys->vtheb, vthib = phys->vthib;
  double weight_b = phys->weight_b;
  int n_global_patches = sub->n_global_patches;

  if (count_only) {
    for (int p = 0; p < psc->n_patches(); p++) {
      nr_particles_by_patch[p] = 2 * (Ne_sheet / n_global_patches + Ne_back / n_global_patches);
    }
    return;
  }
  
  PscMparticlesBase mprts(psc->particles);
  
  // LOAD PARTICLES

  mpi_printf(comm, "Loading particles\n");

  // Do a fast load of the particles

  rngpool = RngPool_create(); // FIXME, should be part of ctor (of struct psc, really)

  int rank;
  MPI_Comm_rank(comm, &rank);
  RngPool_seed(rngpool, rank);
  Rng *rng = RngPool_get(rngpool, 0);

  assert(psc->n_patches() > 0);
  const Grid_t& grid = psc->grid();
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
    mprts->inject_reweight(0, prt);

    ux = Rng_normal(rng, 0, vthi);
    uy = Rng_normal(rng, 0, vthi);
    uz = Rng_normal(rng, 0, vthi);
    d0 = gdri*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udri;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

    prt.u[0] = ux; prt.u[1] = uy; prt.u[2] = uz;
    prt.kind = KIND_ION;
    mprts->inject_reweight(0, prt);
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
    mprts->inject_reweight(0, prt);
    
    prt.u[0] = Rng_normal(rng, 0, vthib);
    prt.u[1] = Rng_normal(rng, 0, vthib);
    prt.u[2] = Rng_normal(rng, 0, vthib);
    prt.kind = KIND_ION;
    mprts->inject_reweight(0, prt);
  }

  mpi_printf(comm, "Finished loading particles\n");
}

// ----------------------------------------------------------------------
// psc_harris_read

static void
psc_harris_read(struct psc *psc, struct mrc_io *io)
{
  psc_read_super(psc, io);
}

// ----------------------------------------------------------------------
// psc_harris_destroy

static void
psc_harris_destroy(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);

  Simulation_delete(sub->sim);
}

// ======================================================================
// psc_harris_ops

struct psc_ops_harris : psc_ops {
  psc_ops_harris() {
    name             = "harris";
    size             = sizeof(struct psc_harris);
    param_descr      = psc_harris_descr;
    destroy          = psc_harris_destroy;
    read             = psc_harris_read;
    setup_particles  = psc_harris_setup_particles;
  }
} psc_harris_ops;

// ======================================================================
// PscHarrisBuilder

struct PscHarrisBuilder
{
  PscHarrisBuilder()
    : psc_(psc_create(MPI_COMM_WORLD))
  {}

  PscHarris* makePscHarris();

  PscHarrisParams params;
  psc* psc_;
};

// ----------------------------------------------------------------------
// PscHarrisBuilder::makePscHarris

PscHarris* PscHarrisBuilder::makePscHarris()
{
  MPI_Comm comm = psc_comm(psc_);
  
  mpi_printf(comm, "*** Setting up...\n");

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


  params.L_di = .5;
  params.Ti_Te = 5.;
  params.nb_n0 = .05;
  params.Tbe_Te = .333;
  params.Tbi_Ti = .333;

  psc_default_dimensionless(psc_);

  psc_->prm.nicell = 1;
  psc_->prm.cfl = 0.99;

  psc_->prm.stats_every = 100;

  auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL },
			{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL },
			{ BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_REFLECTING },
			{ BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_REFLECTING }};
  psc_->bc_ = grid_bc;
 
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
  MPI_Init(argc, argv);
#endif
  libmrc_params_init(argc, argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);

  mrc_class_register_subclass(&mrc_class_psc, &psc_harris_ops);

  auto sim = PscHarrisBuilder{};
  auto harris = sim.makePscHarris();

  psc_view(sim.psc_);
  psc_mparticles_view(sim.psc_->particles);
  psc_mfields_view(sim.psc_->flds);
  
  harris->integrate();
  
  delete harris;

  psc_destroy(sim.psc_);
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
