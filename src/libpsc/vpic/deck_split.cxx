
#include "vpic_iface.h"

#include <cassert>

//////////////////////////////////////////////////////
//
//   Harris Sheet Reconnection - Open Boundary Model
//
//////////////////////////////////////////////////////

// ----------------------------------------------------------------------
// params

struct params : vpic_params, vpic_harris_params {
  // Harris
  double L_di;     // Sheet thickness / ion inertial length
  double Ti_Te;    // Ion temperature / electron temperature
  double nb_n0;    // background plasma density
  double Tbe_Te;   // Ratio of background T_e to Harris T_e
  double Tbi_Ti;   // Ratio of background T_i to Harris T_i
  double bg;       // Guide field
  double theta;
  double Lpert_Lx; // wavelength of perturbation in terms of Lx
  double dbz_b0;   // perturbation in Bz relative to B0
  double nppc;     // Average number of macro particle per cell per species
  int open_bc_x;   // Flag to signal we want to do open boundary condition in x
  int driven_bc_z; // Flag to signal we want to do driven boundary condition in z
};

// ----------------------------------------------------------------------
// globals_diag

struct globals_diag {
  int interval;
  int energies_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  double quota_sec;          // Run quota in seconds

  // state
  int rtoggle;               // enables save of last 2 restart dumps for safety
  // Output variables
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;

  // Vadim: modified restart machinary
  int write_end_restart; // global flag for all to write restart files
};

// ======================================================================

struct globals_physics {
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
};

begin_globals {
  struct params prm;
  struct globals_diag diag;
  struct globals_physics phys;
};

// ======================================================================

static void user_init_params(params *prm);
static void user_init_harris(vpic_simulation* simulation, globals_physics *phys,
			     params *prm, int nproc);
static void user_init_diagnostics(globals_diag *diag, params *prm, globals_physics *phys);
static void user_init_grid(vpic_simulation *simulation, params *prm, globals_physics *phys);
static void user_init_log(vpic_simulation *simulation, params *prm,
			  globals_physics *phys, globals_diag *diag);
static void user_load_fields(vpic_simulation *simulation, params *prm, globals_physics *phys);
static void user_load_particles(vpic_simulation *simulation, params *prm, globals_physics *phys,
				species_t *electron, species_t *ion);
static void user_setup_diagnostics(vpic_simulation *simulation, globals_diag *diag,
				   species_t *electron, species_t *ion);

// ======================================================================
// initialization

begin_initialization {
  params *prm = &global->prm;
  globals_physics *phys = &global->phys;
  globals_diag *diag = &global->diag;

  // FIXME, what a hacky way of passing arguments...
  vpic_params *vpic_prm = reinterpret_cast<vpic_params *>(cmdline_argument[0]);
  vpic_harris_params *vpic_harris_prm = reinterpret_cast<vpic_harris_params *>(cmdline_argument[1]);
  *static_cast<vpic_params *>(prm) = *vpic_prm;
  *static_cast<vpic_harris_params *>(prm) = * vpic_harris_prm;

  user_init_params(prm);

  user_init_harris(this, phys, prm, nproc());

  // Determine the time step
  phys->dg = courant_length(phys->Lx,phys->Ly,phys->Lz,prm->gdims[0],prm->gdims[1],prm->gdims[2]);  // courant length
  phys->dt = prm->cfl_req*phys->dg/phys->c;                       // courant limited time step
  if (phys->wpe * phys->dt > prm->wpedt_max)
    phys->dt = prm->wpedt_max / phys->wpe;  // override timestep if plasma frequency limited

  ///////////////////////////////////////////////
  // Setup high level simulation parameters

  num_step             = int(prm->taui / (phys->wci*phys->dt));
  status_interval      = prm->status_interval;
  sync_shared_interval = status_interval/2;
  clean_div_e_interval = status_interval/2;
  clean_div_b_interval = status_interval/2;

  //////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  user_init_grid(this, prm, phys);
  
  //////////////////////////////////////////////////////////////////////////////
  // Setup materials

  sim_log("Setting up materials. ");

  define_material( "vacuum", 1 );
  material_t * resistive = define_material( "resistive",1,1,1 );

  define_field_array(NULL); // second argument is damp, default to 0

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.


  //////////////////////////////////////////////////////////////////////////////
  // Finalize Field Advance

  sim_log("Finalizing Field Advance");

  double hx = phys->Lx/prm->gdims[0];
  double hz = phys->Lz/prm->gdims[1];

  // Define resistive layer surrounding boundary --> set thickness=0
  // to eliminate this feature
  double thickness = 0;
#define resistive_layer ((prm->open_bc_x && x < hx*thickness) ||	\
			 (prm->open_bc_x && x > phys->Lx-hx*thickness)	\
                         || z <-phys->Lz/2+hz*thickness  || z > phys->Lz/2-hz*thickness )

  if (thickness > 0) {
    sim_log("Setting resistive layer of thickness "<< thickness);
    set_region_material(resistive_layer, resistive, resistive);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Setup the species

  sim_log("Setting up species. ");
  double nmax = 2.0*phys->Ne/nproc();
  double nmovers = 0.1*nmax;
  double sort_method = 1;   // 0=in place and 1=out of place
  species_t *electron = define_species("electron", -phys->ec, phys->me, nmax, nmovers,
				       prm->electron_sort_interval, sort_method);
  species_t *ion = define_species("ion", phys->ec, phys->mi, nmax, nmovers,
				  prm->ion_sort_interval, sort_method);

  ////////////////////////////////////////////////////////////////////////

  user_init_diagnostics(diag, prm, phys);

  user_init_log(this, prm, phys, diag);

  user_load_fields(this, prm, phys);

  user_load_particles(this, prm, phys, electron, ion);

  user_setup_diagnostics(this, diag, electron, ion);

  sim_log("*** Finished with user-specified initialization ***");

  // Upon completion of the initialization, the following occurs:
  // - The synchronization error (tang E, norm B) is computed between domains
  //   and tang E / norm B are synchronized by averaging where discrepancies
  //   are encountered.
  // - The initial divergence error of the magnetic field is computed and
  //   one pass of cleaning is done (for good measure)
  // - The bound charge density necessary to give the simulation an initially
  //   clean divergence e is computed.
  // - The particle momentum is uncentered from u_0 to u_{-1/2}
  // - The user diagnostics are called on the initial state
  // - The physics loop is started
  //
  // The physics loop consists of:
  // - Advance particles from x_0,u_{-1/2} to x_1,u_{1/2}
  // - User particle injection at x_{1-age}, u_{1/2} (use inject_particles)
  // - User current injection (adjust field(x,y,z).jfx, jfy, jfz)
  // - Advance B from B_0 to B_{1/2}
  // - Advance E from E_0 to E_1
  // - User field injection to E_1 (adjust field(x,y,z).ex,ey,ez,cbx,cby,cbz)
  // - Advance B from B_{1/2} to B_1
  // - (periodically) Divergence clean electric field
  // - (periodically) Divergence clean magnetic field
  // - (periodically) Synchronize shared tang e and norm b
  // - Increment the time step
  // - Call user diagnostics
  // - (periodically) Print a status message

} //begin_initialization

// ----------------------------------------------------------------------
// user_init_params

static void user_init_params(params *prm)
{
  // Physics parameters
  prm->L_di    = 0.5;     // Sheet thickness / ion inertial length
  prm->Ti_Te   = 5.0;     // Ion temperature / electron temperature
  prm->nb_n0   = 0.228;   // background plasma density
  prm->Tbe_Te  = 0.7598;  // Ratio of background T_e to Harris T_e
  prm->Tbi_Ti  = 0.3039;  // Ratio of background T_i to Harris T_i
  prm->bg = 0.0;
  prm->theta = 0;         // B0 = Bx

  prm->Lpert_Lx = 1.;     // wavelength of perturbation in terms of Lx
  prm->dbz_b0 = 0.03;     // perturbation in Bz relative to B0

  // Parameters for Open BC model
  prm->open_bc_x = 0;
  prm->driven_bc_z = 0;

  // Numerical parameters

  prm->nppc  = 100; // Average number of macro particle per cell per species

  assert(prm->np[2] <= 2); // For load balance, keep "1" or "2" for Harris sheet

  prm->quota   = 11.0;   // run quota in hours
  prm->quota_check_interval = 100;

  prm->restart_interval = 8000;
}

// ----------------------------------------------------------------------
// user_init_harris

static void user_init_harris(vpic_simulation* simulation, globals_physics *phys,
			     params *prm, int nproc)
{
  // use natural PIC units
  phys->ec   = 1;         // Charge normalization
  phys->me   = 1;         // Mass normalization
  phys->c    = 1;         // Speed of light
  phys->de   = 1;         // Length normalization (electron inertial length)
  phys->eps0 = 1;         // Permittivity of space

  double c = phys->c;
  //derived quantities
  phys->mi = phys->me*prm->mi_me;       // Ion mass
  double Te = phys->me*c*c/(2*phys->eps0*prm->wpe_wce*prm->wpe_wce*(1+prm->Ti_Te)); // Electron temperature
  double Ti = Te*prm->Ti_Te;       // Ion temperature
  phys->vthe = sqrt(Te/phys->me);         // Electron thermal velocity
  phys->vthi = sqrt(Ti/phys->mi);         // Ion thermal velocity
  phys->vtheb = sqrt(prm->Tbe_Te*Te/phys->me);  // normalized background e thermal vel.
  phys->vthib = sqrt(prm->Tbi_Ti*Ti/phys->mi);  // normalized background ion thermal vel.
  phys->wci  = 1.0/(prm->mi_me*prm->wpe_wce);  // Ion cyclotron frequency
  phys->wce  = phys->wci*prm->mi_me;            // Electron cyclotron freqeuncy
  phys->wpe  = phys->wce*prm->wpe_wce;          // electron plasma frequency
  phys->wpi  = phys->wpe/sqrt(prm->mi_me);      // ion plasma frequency
  phys->di   = c/phys->wpi;                      // ion inertial length
  phys->L    = prm->L_di*phys->di;              // Harris sheet thickness
  phys->rhoi_L = sqrt(prm->Ti_Te/(1.0+prm->Ti_Te))/prm->L_di;
  phys->v_A = (phys->wci/phys->wpi)/sqrt(prm->nb_n0); // based on nb

  phys->Lx    = prm->Lx_di*phys->di; // size of box in x dimension
  phys->Ly    = prm->Ly_di*phys->di; // size of box in y dimension
  phys->Lz    = prm->Lz_di*phys->di; // size of box in z dimension

  phys->b0 = phys->me*c*phys->wce/phys->ec; // Asymptotic magnetic field strength
  phys->n0 = phys->me*phys->eps0*phys->wpe*phys->wpe/(phys->ec*phys->ec);  // Peak electron (ion) density
  phys->vdri = 2*c*Ti/(phys->ec*phys->b0*phys->L);   // Ion drift velocity
  phys->vdre = -phys->vdri/(prm->Ti_Te);      // electron drift velocity

  double Lx = phys->Lx, Ly = phys->Ly, Lz = phys->Lz, L = phys->L;
  double Npe_sheet = 2*phys->n0*Lx*Ly*L*tanh(0.5*Lz/L); // N physical e's in sheet
  double Npe_back  = prm->nb_n0*phys->n0 * Ly*Lz*Lx;          // N physical e's in backgrnd
  double Npe       = Npe_sheet + Npe_back;
  phys->Ne         = prm->nppc * prm->gdims[0] * prm->gdims[1] * prm->gdims[2];  // total macro electrons in box
  phys->Ne_sheet   = phys->Ne*Npe_sheet/Npe;
  phys->Ne_back    = phys->Ne*Npe_back/Npe;
  phys->Ne_sheet   = simulation->trunc_granular(phys->Ne_sheet,nproc); // Make it divisible by nproc
  phys->Ne_back    = simulation->trunc_granular(phys->Ne_back, nproc);  // Make it divisible by nproc
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
}

// ----------------------------------------------------------------------
// user_init_diagnostics

static void user_init_diagnostics(globals_diag *diag, params *prm, globals_physics *phys)
{
  diag->rtoggle              = 0;
  diag->quota_sec = prm->quota*3600;  // Run quota in seconds

  diag->interval = int(prm->t_intervali/(phys->wci*phys->dt));
  diag->fields_interval = diag->interval;
  diag->ehydro_interval = diag->interval;
  diag->Hhydro_interval = diag->interval;
  diag->eparticle_interval = 8*diag->interval;
  diag->Hparticle_interval = 8*diag->interval;

  diag->energies_interval = 50;
}

// ----------------------------------------------------------------------
// user_init_grid

static void user_init_grid(vpic_simulation *simulation, params *prm, globals_physics *phys)
{
#define rank simulation->rank
  grid_t *grid = simulation->grid;
  
  // Setup basic grid parameters
  grid->dx = phys->Lx / prm->gdims[0];
  grid->dy = phys->Ly / prm->gdims[1];
  grid->dz = phys->Lz / prm->gdims[2];
  grid->dt = phys->dt;
  grid->cvac = phys->c;
  grid->eps0 = phys->eps0;

  // Define the grid
  simulation->define_periodic_grid(0       , -0.5*phys->Ly, -0.5*phys->Lz,      // Low corner
				   phys->Lx,  0.5*phys->Ly,  0.5*phys->Lz,      // High corner
				   prm->gdims[0], prm->gdims[1], prm->gdims[2], // Resolution
				   prm->np[0], prm->np[1], prm->np[2]);         // Topology

  // Determine which domains area along the boundaries - Use macro from
  // grid/partition.c.
  
# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                 \
    int _ix, _iy, _iz;                                                  \
    _ix  = (rank);                /* ix = ix+gpx*( iy+gpy*iz ) */       \
    _iy  = _ix/int(prm->np[0]);   /* iy = iy+gpy*iz */			\
    _ix -= _iy*int(prm->np[0]);   /* ix = ix */				\
    _iz  = _iy/int(prm->np[1]);   /* iz = iz */				\
    _iy -= _iz*int(prm->np[2]);   /* iy = iy */				\
    (ix) = _ix;                                                         \
    (iy) = _iy;                                                         \
    (iz) = _iz;                                                         \
  } END_PRIMITIVE

  int ix, iy, iz, left=0,right=0,top=0,bottom=0;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );
  if ( ix ==0 ) left=1;
  if ( ix ==prm->np[0]-1 ) right=1;
  if ( iz ==0 ) bottom=1;
  if ( iz ==prm->np[1]-1 ) top=1;

  // ***** Set Field Boundary Conditions *****
  if (prm->open_bc_x) {
    sim_log("Absorbing fields on X-boundaries");
    if (left ) simulation->set_domain_field_bc( BOUNDARY(-1,0,0), absorb_fields );
    if (right) simulation->set_domain_field_bc( BOUNDARY( 1,0,0), absorb_fields );
  }
  
  sim_log("Conducting fields on Z-boundaries");
  if (bottom) simulation->set_domain_field_bc( BOUNDARY(0,0,-1), pec_fields );
  if (top   ) simulation->set_domain_field_bc( BOUNDARY( 0,0,1), pec_fields );

  // ***** Set Particle Boundary Conditions *****
  if (prm->driven_bc_z) {
    sim_log("Absorb particles on Z-boundaries");
    if (bottom) simulation->set_domain_particle_bc( BOUNDARY(0,0,-1), absorb_particles );
    if (top   ) simulation->set_domain_particle_bc( BOUNDARY(0,0, 1), absorb_particles );
  } else {
    sim_log("Reflect particles on Z-boundaries");
    if (bottom) simulation->set_domain_particle_bc( BOUNDARY(0,0,-1), reflect_particles );
    if (top   ) simulation->set_domain_particle_bc( BOUNDARY(0,0, 1), reflect_particles );
  }
  if (prm->open_bc_x) {
    sim_log("Absorb particles on X-boundaries");
    if (left)  simulation->set_domain_particle_bc( BOUNDARY(-1,0,0), absorb_particles );
    if (right) simulation->set_domain_particle_bc( BOUNDARY( 1,0,0), absorb_particles );
  }

#undef rank
}

// ----------------------------------------------------------------------
// user_init_log

static void user_init_log(vpic_simulation *simulation, params *prm,
			  globals_physics *phys, globals_diag *diag)
{
#define rank simulation->rank
#define nproc simulation->nproc

  sim_log( "***********************************************" );
  sim_log("* Topology:                       " << prm->np[0]
	  << " " << prm->np[1] << " " << prm->np[2]);
  sim_log ( "tanhf = " << phys->tanhf );
  sim_log ( "L_di   = " << prm->L_di );
  sim_log ( "rhoi/L   = " << phys->rhoi_L );
  sim_log ( "Ti/Te = " << prm->Ti_Te ) ;
  sim_log ( "nb/n0 = " << prm->nb_n0 ) ;
  sim_log ( "wpe/wce = " << prm->wpe_wce );
  sim_log ( "mi/me = " << prm->mi_me );
  sim_log ( "theta = " << prm->theta );
  sim_log ( "Lpert/Lx = " << prm->Lpert_Lx );
  sim_log ( "dbz/b0 = " << prm->dbz_b0 );
  sim_log ( "taui = " << prm->taui );
  sim_log ( "t_intervali = " << prm->t_intervali );
  sim_log ( "interval = " << diag->interval );
  sim_log ( "num_step = " << simulation->num_step );
  sim_log ( "Lx/di = " << phys->Lx/phys->di );
  sim_log ( "Lx/de = " << phys->Lx/phys->de );
  sim_log ( "Ly/di = " << phys->Ly/phys->di );
  sim_log ( "Ly/de = " << phys->Ly/phys->de );
  sim_log ( "Lz/di = " << phys->Lz/phys->di );
  sim_log ( "Lz/de = " << phys->Lz/phys->de );
  sim_log ( "nx = " << prm->gdims[0] );
  sim_log ( "ny = " << prm->gdims[1] );
  sim_log ( "nz = " << prm->gdims[2] );
  sim_log ( "courant = " << phys->c*phys->dt/phys->dg );
  sim_log ( "nproc = " << nproc ()  );
  sim_log ( "nppc = " << prm->nppc );
  sim_log ( "b0 = " << phys->b0 );
  sim_log ( "v_A (based on nb) = " << phys->v_A );
  sim_log ( "di = " << phys->di );
  sim_log ( "Ne = " << phys->Ne );
  sim_log ( "Ne_sheet = " << phys->Ne_sheet );
  sim_log ( "Ne_back = " << phys->Ne_back );
  sim_log ( "total # of particles = " << 2*phys->Ne );
  sim_log ( "dt*wpe = " << phys->wpe*phys->dt );
  sim_log ( "dt*wce = " << phys->wce*phys->dt );
  sim_log ( "dt*wci = " << phys->wci*phys->dt );
  sim_log ( "energies_interval: " << diag->energies_interval );
  sim_log ( "dx/de = " << phys->Lx/(phys->de*prm->gdims[0]) );
  sim_log ( "dy/de = " << phys->Ly/(phys->de*prm->gdims[1]) );
  sim_log ( "dz/de = " << phys->Lz/(phys->de*prm->gdims[2]) );
  sim_log ( "dx/rhoi = " << (phys->Lx/prm->gdims[0])/(phys->vthi/phys->wci)  );
  sim_log ( "dx/rhoe = " << (phys->Lx/prm->gdims[0])/(phys->vthe/phys->wce)  );
  sim_log ( "L/debye = " << phys->L/(phys->vthe/phys->wpe)  );
  sim_log ( "dx/debye = " << (phys->Lx/prm->gdims[0])/(phys->vthe/phys->wpe)  );
  sim_log ( "n0 = " << phys->n0 );
  sim_log ( "vthi/c = " << phys->vthi/phys->c );
  sim_log ( "vthe/c = " << phys->vthe/phys->c );
  sim_log ( "vdri/c = " << phys->vdri/phys->c );
  sim_log ( "vdre/c = " << phys->vdre/phys->c );
  sim_log ( "Open BC in x?   = " << prm->open_bc_x );
  sim_log ( "Driven BC in z? = " << prm->driven_bc_z );

  // Dump simulation information to file "info"
  if (rank() == 0 ) {
    FileIO fp_info;
    if ( ! (fp_info.open("info", io_write)==ok) ) ERROR(("Cannot open file."));
    fp_info.print("           ***** Simulation parameters ***** \n");
    fp_info.print("              L/di   =               %e\n", prm->L_di);
    fp_info.print("              L/de   =               %e\n", phys->L/phys->de);
    fp_info.print("              rhoi/L =               %e\n", phys->rhoi_L);
    fp_info.print("              Ti/Te  =               %e\n", prm->Ti_Te );
    fp_info.print("              Tbi/Ti =               %e\n", prm->Tbi_Ti );
    fp_info.print("              Tbe/Te =               %e\n", prm->Tbe_Te );
    fp_info.print("              nb/n0 =                %e\n", prm->nb_n0 );
    fp_info.print("              wpe/wce =              %e\n", prm->wpe_wce );
    fp_info.print("              mi/me =                %e\n", prm->mi_me );
    fp_info.print("              theta =                %e\n", prm->theta );
    fp_info.print("              taui =                 %e\n", prm->taui );
    fp_info.print("              num_step =             %i\n", simulation->num_step );
    fp_info.print("              Lx/de =                %e\n", phys->Lx/phys->de );
    fp_info.print("              Ly/de =                %e\n", phys->Ly/phys->de );
    fp_info.print("              Lz/de =                %e\n", phys->Lz/phys->de );
    fp_info.print("              Lx/di =                %e\n", phys->Lx/phys->di );
    fp_info.print("              Ly/di =                %e\n", phys->Ly/phys->di );
    fp_info.print("              Lz/di =                %e\n", phys->Lz/phys->di );
    fp_info.print("              nx =                   %e\n", prm->gdims[0] );
    fp_info.print("              ny =                   %e\n", prm->gdims[1] );
    fp_info.print("              nz =                   %e\n", prm->gdims[2] );
    fp_info.print("              courant =              %e\n", phys->c*phys->dt/phys->dg );
    fp_info.print("              nproc =                %i\n", nproc() );
    fp_info.print("              nppc =                 %e\n", prm->nppc );
    fp_info.print("              b0 =                   %e\n", phys->b0 );
    fp_info.print("              v_A (based on nb) =    %e\n", phys->v_A );
    fp_info.print("              di =                   %e\n", phys->di );
    fp_info.print("              Ne =                   %e\n", phys->Ne );
    fp_info.print("              Ne_sheet =             %e\n", phys->Ne_sheet );
    fp_info.print("              Ne_back =              %e\n", phys->Ne_back );
    fp_info.print("              total # of particles = %e\n", 2*phys->Ne );
    fp_info.print("              dt*wpe =               %e\n", phys->wpe*phys->dt );
    fp_info.print("              dt*wce =               %e\n", phys->wce*phys->dt );
    fp_info.print("              dt*wci =               %e\n", phys->wci*phys->dt );
    fp_info.print("              energies_interval:     %i\n", diag->energies_interval);
    fp_info.print("              dx/de =                %e\n", phys->Lx/(phys->de*prm->gdims[0]) );
    fp_info.print("              dy/de =                %e\n", phys->Ly/(phys->de*prm->gdims[1]) );
    fp_info.print("              dz/de =                %e\n", phys->Lz/(phys->de*prm->gdims[2]) );
    fp_info.print("              L/debye =              %e\n", phys->L/(phys->vthe/phys->wpe) );
    fp_info.print("              dx/rhoi =              %e\n", (phys->Lx/prm->gdims[0])/(phys->vthi/phys->wci) );
    fp_info.print("              dx/rhoe =              %e\n", (phys->Lx/prm->gdims[0])/(phys->vthe/phys->wce) );
    fp_info.print("              dx/debye =             %e\n", (phys->Lx/prm->gdims[0])/(phys->vthe/phys->wpe) );
    fp_info.print("              n0 =                   %e\n", phys->n0 );
    fp_info.print("              vthi/c =               %e\n", phys->vthi/phys->c );
    fp_info.print("              vthe/c =               %e\n", phys->vthe/phys->c );
    fp_info.print("              vdri/c =               %e\n", phys->vdri/phys->c );
    fp_info.print("              vdre/c =               %e\n", phys->vdre/phys->c );
    fp_info.print("              ***************************\n");
    fp_info.close();
  }

#undef rank
#undef nproc
}

// ----------------------------------------------------------------------
// user_load_fields

static void user_load_fields(vpic_simulation *simulation, params *prm, globals_physics *phys)
{
  double cs = cos(prm->theta), sn = sin(prm->theta);
  double L = phys->L, Lx = phys->Lx, Lz = phys->Lz, Lpert = phys->Lpert;
  double b0 = phys->b0, bg = prm->bg, dbx = phys->dbx, dbz = phys->dbz;
  grid_t *grid = simulation->grid;
#define rank simulation->rank
#define field simulation->field
  
  sim_log( "Loading fields" );
  set_region_field( everywhere, 0, 0, 0,                    // Electric field
    cs*b0*tanh(z/L)+dbx*cos(2.0*M_PI*(x-0.5*Lx)/Lpert)*sin(M_PI*z/Lz), //Bx
    -sn*b0*tanh(z/L) + b0*bg, //By
    dbz*cos(M_PI*z/Lz)*sin(2.0*M_PI*(x-0.5*Lx)/Lpert) ); // Bz

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

#undef field
#undef rank
}

// ----------------------------------------------------------------------
// user_load_particles

static void user_load_particles(vpic_simulation *simulation, params *prm,
				globals_physics *phys, species_t *electron, species_t *ion)
{
  double cs = cos(prm->theta), sn = sin(prm->theta);
  double Ne_sheet = phys->Ne_sheet, vthe = phys->vthe, vthi = phys->vthi;
  double weight_s = phys->weight_s;
  double tanhf = phys->tanhf, L = phys->L;
  double gdre = phys->gdre, udre = phys->udre, gdri = phys->gdri, udri = phys->udri;
  double Ne_back = phys->Ne_back, vtheb = phys->vtheb, vthib = phys->vthib;
  double weight_b = phys->weight_b;
  grid_t *grid = simulation->grid;
#define rank simulation->rank
#define nproc simulation->nproc
#define rng simulation->rng
#define uniform simulation->uniform
#define normal simulation->normal
#define inject_particle simulation->inject_particle
#define seed_entropy simulation->seed_entropy
  
  // LOAD PARTICLES

  sim_log( "Loading particles" );

  // Do a fast load of the particles

  seed_entropy( rank() );  //Generators desynchronized
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  // Load Harris population

  sim_log( "-> Main Harris Sheet" );

  repeat ( Ne_sheet/nproc() ) {
    double x, y, z, ux, uy, uz, d0 ;

    do {
      z = L*atanh(uniform( rng(0), -1, 1)*tanhf);
    } while( z<= zmin || z>=zmax );
    x = uniform( rng(0), xmin, xmax );
    y = uniform( rng(0), ymin, ymax );

    // inject_particles() will return an error for particles no on this
    // node and will not inject particle locally

    ux = normal( rng(0), 0, vthe);
    uy = normal( rng(0), 0, vthe);
    uz = normal( rng(0), 0, vthe);
    d0 = gdre*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udre;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

    inject_particle(electron, x, y, z, ux, uy, uz, weight_s, 0, 0 );

    ux = normal( rng(0), 0, vthi);
    uy = normal( rng(0), 0, vthi);
    uz = normal( rng(0), 0, vthi);
    d0 = gdri*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udri;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

    inject_particle(ion, x, y, z, ux, uy, uz, weight_s, 0, 0 );
  }

  sim_log( "-> Background Population" );

  repeat ( Ne_back/nproc() ) {

    double x = uniform( rng(0), xmin, xmax );
    double y = uniform( rng(0), ymin, ymax );
    double z = uniform( rng(0), zmin, zmax );

    inject_particle( electron, x, y, z,
                     normal( rng(0), 0, vtheb),
                     normal( rng(0), 0, vtheb),
                     normal( rng(0), 0, vtheb),
                     weight_b, 0, 0);

    inject_particle( ion, x, y, z,
                     normal( rng(0), 0, vthib),
                     normal( rng(0), 0, vthib),
                     normal( rng(0), 0, vthib),
                     weight_b, 0 ,0 );
  }

  sim_log( "Finished loading particles" );

#undef rank
#undef nproc
#undef rng
#undef uniform
#undef normal
#undef inject_particle
#undef seed_entropy
}

// ----------------------------------------------------------------------
// user_setup_diagnostics

static void user_setup_diagnostics(vpic_simulation *simulation, globals_diag *diag,
				   species_t *electron, species_t *ion)
{
#define create_field_list simulation->create_field_list
#define create_hydro_list simulation->create_hydro_list
#define rank simulation->rank

  /*--------------------------------------------------------------------------
   * New dump definition
   *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Set data output format
   *
   * This option allows the user to specify the data format for an output
   * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
   * format is the native storage format for data in VPIC.  For field data,
   * this looks something like:
   *
   *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
   *
   * Banded data format stores all data of a particular state variable as a
   * contiguous array, and is easier for ParaView to process efficiently.
   * Banded data looks like:
   *
   *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
   *
   *------------------------------------------------------------------------*/

  diag->fdParams.format = band;
  sim_log ( "Fields output format = band" );

  diag->hedParams.format = band;
  sim_log ( "Electron species output format = band" );

  diag->hHdParams.format = band;
  sim_log ( "Ion species output format = band" );

  /*--------------------------------------------------------------------------
   * Set stride
   *
   * This option allows data down-sampling at output.  Data are down-sampled
   * in each dimension by the stride specified for that dimension.  For
   * example, to down-sample the x-dimension of the field data by a factor
   * of 2, i.e., half as many data will be output, select:
   *
   *   global->fdParams.stride_x = 2;
   *
   * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
   * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
   * Setting the strides in x and y to equal 2 results in an output grid of
   * nx = 4, ny = 4, with actual extents 6x6.
   *
   * G G G G G G G G G
   * G X X X X X X X G
   * G X X X X X X X G         G G G G G G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G   ==>   G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G G G G G G
   * G G G G G G G G G
   *
   * Note that grid extents in each dimension must be evenly divisible by
   * the stride for that dimension:
   *
   *   nx = 150;
   *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
   *
   *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
   *------------------------------------------------------------------------*/

  // relative path to fields data from global header
  sprintf(diag->fdParams.baseDir, "fields");

  // base file name for fields output
  sprintf(diag->fdParams.baseFileName, "fields");

  diag->fdParams.stride_x = 1;
  diag->fdParams.stride_y = 1;
  diag->fdParams.stride_z = 1;

  // add field parameters to list
  diag->outputParams.push_back(&diag->fdParams);

  sim_log ( "Fields x-stride " << diag->fdParams.stride_x );
  sim_log ( "Fields y-stride " << diag->fdParams.stride_y );
  sim_log ( "Fields z-stride " << diag->fdParams.stride_z );

  // relative path to electron species data from global header
  sprintf(diag->hedParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(diag->hedParams.baseFileName, "ehydro");

  diag->hedParams.stride_x = 1;
  diag->hedParams.stride_y = 1;
  diag->hedParams.stride_z = 1;

  // add electron species parameters to list
  diag->outputParams.push_back(&diag->hedParams);

  sim_log ( "Electron species x-stride " << diag->hedParams.stride_x );
  sim_log ( "Electron species y-stride " << diag->hedParams.stride_y );
  sim_log ( "Electron species z-stride " << diag->hedParams.stride_z );

  // relative path to electron species data from global header
  sprintf(diag->hHdParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(diag->hHdParams.baseFileName, "Hhydro");

  diag->hHdParams.stride_x = 1;
  diag->hHdParams.stride_y = 1;
  diag->hHdParams.stride_z = 1;

  sim_log ( "Ion species x-stride " << diag->hHdParams.stride_x );
  sim_log ( "Ion species y-stride " << diag->hHdParams.stride_y );
  sim_log ( "Ion species z-stride " << diag->hHdParams.stride_z );

  // add electron species parameters to list
  diag->outputParams.push_back(&diag->hHdParams);

  /*--------------------------------------------------------------------------
   * Set output fields
   *
   * It is now possible to select which state-variables are output on a
   * per-dump basis.  Variables are selected by passing an or-list of
   * state-variables by name.  For example, to only output the x-component
   * of the electric field and the y-component of the magnetic field, the
   * user would call output_variables like:
   *
   *   global->fdParams.output_variables( ex | cby );
   *
   * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
   * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
   * 'output_variables' WILL HAVE NO EFFECT.
   *
   * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
   * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
   *
   * For convenience, the output variable 'all' is defined:
   *
   *   global->fdParams.output_variables( all );
   *------------------------------------------------------------------------*/
  /* CUT AND PASTE AS A STARTING POINT
   * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE

   output_variables( all );

   output_variables( electric | div_e_err | magnetic | div_b_err |
                     tca      | rhob      | current  | rhof |
                     emat     | nmat      | fmat     | cmat );

   output_variables( current_density  | charge_density |
                     momentum_density | ke_density     | stress_tensor );
   */

  //global->fdParams.output_variables( electric | magnetic );
  //global->hedParams.output_variables( current_density | charge_density
  //                                    | stress_tensor );
  //global->hHdParams.output_variables( current_density | charge_density );
  //                                    | stress_tensor );

  diag->fdParams.output_variables( all );
  diag->hedParams.output_variables( all );
  diag->hHdParams.output_variables( all );

  /*--------------------------------------------------------------------------
   * Convenience functions for simlog output
   *------------------------------------------------------------------------*/

  char varlist[512];
  create_field_list(varlist, diag->fdParams);

  sim_log ( "Fields variable list: " << varlist );

  create_hydro_list(varlist, diag->hedParams);

  sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, diag->hHdParams);

  sim_log ( "Ion species variable list: " << varlist );

#undef create_field_list
#undef create_hydro_list
#undef rank
}

#define should_dump(x)                                                  \
  (diag->x##_interval>0 && remainder(step(), diag->x##_interval) == 0)

begin_diagnostics {

  /*--------------------------------------------------------------------------
   * NOTE: YOU CANNOT DIRECTLY USE C FILE DESCRIPTORS OR SYSTEM CALLS ANYMORE
   *
   * To create a new directory, use:
   *
   *   dump_mkdir("full-path-to-directory/directoryname")
   *
   * To open a file, use: FileIO class
   *
   * Example for file creation and use:
   *
   *   // declare file and open for writing
   *   // possible modes are: io_write, io_read, io_append,
   *   // io_read_write, io_write_read, io_append_read
   *   FileIO fileIO;
   *   FileIOStatus status;
   *   status= fileIO.open("full-path-to-file/filename", io_write);
   *
   *   // formatted ASCII  output
   *   fileIO.print("format string", varg1, varg2, ...);
   *
   *   // binary output
   *   // Write n elements from array data to file.
   *   // T is the type, e.g., if T=double
   *   // fileIO.write(double * data, size_t n);
   *   // All basic types are supported.
   *   fileIO.write(T * data, size_t n);
   *
   *   // close file
   *   fileIO.close();
   *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Data output directories
   * WARNING: The directory list passed to "global_header" must be
   * consistent with the actual directories where fields and species are
   * output using "field_dump" and "hydro_dump".
   *
   * DIRECTORY PATHES SHOULD BE RELATIVE TO
   * THE LOCATION OF THE GLOBAL HEADER!!!
   *------------------------------------------------------------------------*/

  struct globals_diag *diag = &global->diag;
  struct params *prm = &global->prm;

  /*--------------------------------------------------------------------------
   * Normal rundata dump
   *------------------------------------------------------------------------*/
  if(step()==0) {
    dump_mkdir("fields");
    dump_mkdir("hydro");
    dump_mkdir("rundata");
    dump_mkdir("injectors");
    dump_mkdir("restart1");  // 1st backup
    dump_mkdir("restart2");  // 2nd backup
    dump_mkdir("particle");

    dump_grid("rundata/grid");
    dump_materials("rundata/materials");
    dump_species("rundata/species");
    global_header("global", diag->outputParams);
  } // if

  /*--------------------------------------------------------------------------
   * Normal rundata energies dump
   *------------------------------------------------------------------------*/
  if(should_dump(energies)) {
    dump_energies("rundata/energies", step() == 0 ? 0 : 1);
  } // if

  /*--------------------------------------------------------------------------
   * Field data output
   *------------------------------------------------------------------------*/

  if(step() == -1 || should_dump(fields)) field_dump(diag->fdParams);

  /*--------------------------------------------------------------------------
   * Electron species output
   *------------------------------------------------------------------------*/

  if(should_dump(ehydro)) hydro_dump("electron", diag->hedParams);

  /*--------------------------------------------------------------------------
   * Ion species output
   *------------------------------------------------------------------------*/

  if(should_dump(Hhydro)) hydro_dump("ion", diag->hHdParams);

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  if(step() && !(step()%prm->restart_interval)) {
    if(!diag->rtoggle) {
      diag->rtoggle = 1;
      checkpt("restart1/restart", 0);
    }
    else {
      diag->rtoggle = 0;
      checkpt("restart2/restart", 0);
    } // if
  } // if

  // Dump particle data

  char subdir[36];
  if ( should_dump(eparticle) && step() !=0
       && step() > 56*(diag->fields_interval)  ) {
    // if ( should_dump(eparticle) && step() !=0 ) {
    sprintf(subdir,"particle/T.%lld",step());
    dump_mkdir(subdir);
    sprintf(subdir,"particle/T.%lld/eparticle",step());
    dump_particles("electron", subdir);
  }

  if ( should_dump(Hparticle) && step() !=0
       && step() > 56*(diag->fields_interval)  ) {
    sprintf(subdir,"particle/T.%lld/Hparticle",step());
    dump_particles("ion", subdir);
  }

  // Shut down simulation when wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will
  // be synchronized across processors. Note that this is only checked every
  // few timesteps to eliminate the expensive mp_elapsed call from every
  // timestep. mp_elapsed has an ALL_REDUCE in it!

  if (( step()>0 && prm->quota_check_interval>0
        && (step()&prm->quota_check_interval)==0)
      || (diag->write_end_restart) ) {
    if ( (diag->write_end_restart) ) {

      diag->write_end_restart = 0; // reset restart flag

      sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
      double dumpstart = uptime();

      checkpt("restart0/restart",0);

      mp_barrier(  ); // Just to be safe
      sim_log( "Restart dump restart completed." );
      double dumpelapsed = uptime() - dumpstart;
      sim_log("Restart duration "<< dumpelapsed);
      exit(0); // Exit or abort?
    }
    if( uptime() > diag->quota_sec ) diag->write_end_restart = 1;
  }

} // end diagnostics

begin_particle_injection {
}

begin_current_injection {
}

begin_field_injection {
}

begin_particle_collisions {
}
