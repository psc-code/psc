//////////////////////////////////////////////////////
//
//   Harris Sheet Reconnection - Open Boundary Model
//
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////

// structure to hold the data for energy diagnostics
struct edata {
  species_id sp_id;       /* species id */
  double       vth;       /* thermal energy */
  char  fname[256];       /* file to save data */
};

begin_globals {

  int restart_interval;
  int energies_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  int quota_check_interval;  // How frequently to check if quota exceeded

  int rtoggle;               // enables save of last 2 restart dumps for safety
  double quota_sec;          // Run quota in seconds
  double b0;                 // B0
  double bg;                 // Guide field
  double v_A;
  double topology_x;       // domain topology
  double topology_y;
  double topology_z;

  // parameters for the collision model
  int    ee_collisions; // Flag to signal we want to do e-e collisions.
  int    ei_collisions; // Flag to signal we want to do e-i collisions.
  int    ii_collisions; // Flag to signal we want to do i-i collisions.

  double cvar; // Base variance (dimensionless) used in particle collision
  double nppc_max; // Max number of particles/cell (used to def key array size).
  int tstep_coll;  // Collision interval (=multiple of sort interval).
  double Z;
  double mi_me;
  double wpewce;

  // Variables for Open BC Model
  int open_bc_x; // Flag to signal we want to do open boundary condition in x
  int driven_bc_z; // Flag to signal we want to do driven boundary condition in z
  double nb;      // Background density
  int nsp;        // Number of Species
  double vth[2];  // Thermal velocity of Harris components
  double vthb[2]; // Thermal velocity of background components
  double q[2];    // Species charge
  double L_de;    // Initial Harris sheet thickness
  double uf[2];   // Initial Fluid Drift in Harris
  double nfac; // Normalization factor to convert particles per cell to density
  double rin[3];  // Relaxation parameter for inflow boundary moments
  double rout[3]; // Relaxation parameter for outlfow boundary moments
  double sort[2]; // Intervals where we know particles are sorted
  double edrive;  // Drive field for inflow boundary
  double tdrive;
  int left,right,top,bottom;  // Keep track of boundary domains
  // Moments for bottom injectors
  double *nbot, *ubot, *pbot, *bbot, *fbot;
  // Moments for top injectors
  double *ntop, *utop, *ptop, *btop, *ftop;
  // Moments for left injectors
  double *nleft, *uleft, *pleft, *bleft, *fleft;
  // Moments for right injectors
  double *nright, *uright, *pright, *bright, *fright;

  // Output variables
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;

  // Variables for the energy diagnostics
  edata ede;             // parameters for electron species
  edata edi;             // parameters for ion species
  double emax;           // maximum energy (in units of me*c**2)
  int nex;               // number of energy bins

  //Vadim: modified restart machinary
  int write_restart;     // global flag for all to write restart files
  int write_end_restart; // global flag for all to write restart files
};

begin_initialization {
  // use natural PIC units
  double ec   = 1;         // Charge normalization
  double me   = 1;         // Mass normalization
  double c    = 1;         // Speed of light
  double de   = 1;         // Length normalization (electron inertial length)
  double eps0 = 1;         // Permittivity of space

  double cfl_req   = 0.99;  // How close to Courant should we try to run
  double wpedt_max = 0.36;  // Max dt allowed if Courant not too restrictive
  double damp      = 0.0;   // Level of radiation damping
  int rng_seed     = 1;     // Random number seed increment

  // Physics parameters
  double mi_me   = 25.0;  // Ion mass / electron mass
  double L_di    = 0.5;    // Sheet thickness / ion inertial length
  double Ti_Te   = 5.0;    // Ion temperature / electron temperature
  double Z   = 1.0;      // Ion charge
  double nb_n0   = 0.228;   // background plasma density
  double Tbe_Te  = 0.7598;  // Ratio of background T_e to Harris T_e
  double Tbi_Ti  = 0.3039;  // Ratio of background T_i to Harris T_i
  double wpe_wce = 2.0;    // electron plasma freq / electron cyclotron freq
  double bg = 0.0;
  double theta   = 0;      // B0 = Bx
  double taui    = 1.1;//100;    // simulation wci's to run

  double Lpert_Lx = 1.; // wavelength of perturbation in terms of Lx
  double dbz_b0 = 0.03; // perturbation in Bz relative to B0

  int restart_interval = 8000;
  double t_intervali = 1; // output interval in terms of 1/wci

  double quota   = 11.0;   // run quota in hours
  double quota_sec = quota*3600;  // Run quota in seconds

  double cs   = cos(theta);
  double sn   = sin(theta);

  //derived qunatities
  double mi = me*mi_me;       // Ion mass
  double Te = me*c*c/(2*eps0*wpe_wce*wpe_wce*(1+Ti_Te)); // Electron temperature
  double Ti = Te*Ti_Te;       // Ion temperature
  double vthe = sqrt(Te/me);                        // Electron thermal velocity
  double vthi = sqrt(Ti/mi);  // Ion thermal velocity
  double vtheb = sqrt(Tbe_Te*Te/me);  // normalized background e thermal vel.
  double vthib = sqrt(Tbi_Ti*Ti/mi);  // normalized background ion thermal vel.
  double wci  = 1.0/(mi_me*wpe_wce);  // Ion cyclotron frequency
  double wce  = wci*mi_me;            // Electron cyclotron freqeuncy
  double wpe  = wce*wpe_wce;          // electron plasma frequency
  double wpi  = wpe/sqrt(mi_me);      // ion plasma frequency
  double di   = c/wpi;                // ion inertial length
  double L    = L_di*di;              // Harris sheet thickness
  double rhoi_L = sqrt(Ti_Te/(1.0+Ti_Te))/L_di;
  double v_A= (wci/wpi)/sqrt(nb_n0); // based on nb

  double ion_sort_interval = 1000; // Injector moments also updated
  double electron_sort_interval = 1000; // Injector moments also updated

  // Parameters for Open BC model
  // Relaxation - density, velocity + particle flux, pressure tensor
  int open_bc_x = 0;
  int driven_bc_z = 0;
  double rin[3] =  {0.000, 0.06, 0.000};
  double rout[3] = {0.002, 0.002, 0.002};
  double edrive = 0.08*v_A/(wpe_wce);    // edrive = 0 gives undriven limit
  //double edrive = 0.0099;    // Setting edrive = 0 will give undriven limit

  double tdrive = 32000.0;
  double sort_interval = 10;  // Injector moments also updated at this interval

  // Numerical parameters
  double nppc  = 100; // Average number of macro particle per cell per species

  double Lx    = 25.6*di;      // size of box in x dimension
  double Ly    = 1*di; // size of box in y dimension
  double Lz    = 12.8*di;      // size of box in z dimension

  double topology_x = 4;  // Number of domains in x, y, and z
  double topology_y = 1;
  double topology_z = 1;  // For load balance, keep "1" or "2" for Harris sheet

  double nx = 64;
  double ny = 1;
  double nz = 32;

  double hx = Lx/nx;
  double hy = Ly/ny;
  double hz = Lz/nz;

  double b0 = me*c*wce/ec; // Asymptotic magnetic field strength
  double n0 = me*eps0*wpe*wpe/(ec*ec);  // Peak electron (ion) density
  double vdri = 2*c*Ti/(ec*b0*L);   // Ion drift velocity
  double vdre = -vdri/(Ti_Te);      // electron drift velocity

  double Npe_sheet = 2*n0*Lx*Ly*L*tanh(0.5*Lz/L); // N physical e's in sheet
  double Npe_back  = nb_n0*n0*Ly*Lz*Lx;           // N physical e's in backgrnd
  double Npe       = Npe_sheet + Npe_back;
  double Ne        = nppc*nx*ny*nz;  // total macro electrons in box
  double Ne_sheet  = Ne*Npe_sheet/Npe;
  double Ne_back   = Ne*Npe_back/Npe;
  Ne_sheet = trunc_granular(Ne_sheet,nproc()); // Make it divisible by nproc
  Ne_back  = trunc_granular(Ne_back,nproc());  // Make it divisible by nproc
  Ne = Ne_sheet + Ne_back;
  double qe_s = -ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  double qi_s =  ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  double weight_s = ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  double qe_b = -ec*Npe_back/Ne_back;  // Charge per macro electron
  double qi_b =  ec*Npe_back/Ne_back;  // Charge per macro electron
  double weight_b =  ec*Npe_back/Ne_back;  // Charge per macro electron
  double nfac = qi_s/(hx*hy*hz);       // Convert density to particles per cell

  double gdri = 1/sqrt(1-vdri*vdri/(c*c));  // gamma of ion drift frame
  double gdre = 1/sqrt(1-vdre*vdre/(c*c)); // gamma of electron drift frame
  double udri = vdri*gdri;                 // 4-velocity of ion drift frame
  double udre = vdre*gdre;                 // 4-velocity of electron drift frame
  double tanhf = tanh(0.5*Lz/L);
  double Lpert = Lpert_Lx*Lx; // wavelength of perturbation
  double dbz =  dbz_b0*b0; // Perturbation in Bz relative to Bo (Only change here)
  double dbx = -dbz*Lpert/(2.0*Lz); // Set Bx perturbation so that div(B) = 0

  // Determine the time step
  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz);  // courant length
  double dt = cfl_req*dg/c;                       // courant limited time step
  if( wpe*dt>wpedt_max)
    dt=wpedt_max/wpe;  // override timestep if plasma frequency limited

  // Intervals for output
  int energies_interval = 50;
  int interval = int(t_intervali/(wci*dt));
  int fields_interval = interval;
  int ehydro_interval = interval;
  int Hhydro_interval = interval;
  int eparticle_interval = 8*interval;
  int Hparticle_interval = 8*interval;
  int quota_check_interval     = 100;

  // Collision parameters
  // In CGS, variance of tan theta =
  // 2 pi e^4 n_e dt_coll loglambda / (m_ab^2 u^3)

  int ii_collisions = 0; // do collisions between the corresponding species
  int ee_collisions = 0;
  int ei_collisions = 0;

  int tstep_coll = (int) sort_interval;  // How frequently to do collisions
  double dt_coll = dt*(tstep_coll);      // in (1/wpe)
  double nuei_wce = 0.05;
  double cvar = dt_coll*3.0*sqrt(2.0*M_PI)/4.0*pow(vthe,3)*nuei_wce/wpe_wce;
  // Max possible number of particles/cell of ea. species
  // (to size the key array in collision handler)
  double nppc_max   = 20*nppc;

  // Determine which domains area along the boundaries - Use macro from
  // grid/partition.c.

# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                 \
    int _ix, _iy, _iz;                                                  \
    _ix  = (rank);                /* ix = ix+gpx*( iy+gpy*iz ) */       \
    _iy  = _ix/int(topology_x);   /* iy = iy+gpy*iz */                  \
    _ix -= _iy*int(topology_x);   /* ix = ix */                         \
    _iz  = _iy/int(topology_y);   /* iz = iz */                         \
    _iy -= _iz*int(topology_y);   /* iy = iy */                         \
    (ix) = _ix;                                                         \
    (iy) = _iy;                                                         \
    (iz) = _iz;                                                         \
  } END_PRIMITIVE

  int ix, iy, iz, left=0,right=0,top=0,bottom=0;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );
  if ( ix ==0 ) left=1;
  if ( ix ==topology_x-1 ) right=1;
  if ( iz ==0 ) bottom=1;
  if ( iz ==topology_z-1 ) top=1;

  ///////////////////////////////////////////////
  // Setup high level simulation parameters
  num_step             = int(taui/(wci*dt));
  status_interval      = 100;
  sync_shared_interval = status_interval/2;
  clean_div_e_interval = status_interval/2;
  clean_div_b_interval = status_interval/2;

  global->restart_interval     = restart_interval;
  global->energies_interval    = energies_interval;
  global->fields_interval      = fields_interval;
  global->ehydro_interval      = ehydro_interval;
  global->Hhydro_interval      = Hhydro_interval;
  global->eparticle_interval   = eparticle_interval;
  global->Hparticle_interval   = Hparticle_interval;
  global->quota_check_interval = quota_check_interval;
  global->quota_sec            = quota_sec;
  global->rtoggle              = 0;

  global->b0  = b0;
  global->bg  = bg;
  global->v_A  = v_A;
  global->edrive  = edrive;
  global->tdrive  = tdrive;

  global->topology_x  = topology_x;
  global->topology_y  = topology_y;
  global->topology_z  = topology_z;

  global->left = left;
  global->right = right;
  global->top = top;
  global->bottom = bottom;

  // Parameters for the open boundary model
  global->open_bc_x = open_bc_x;
  global->driven_bc_z = driven_bc_z;
  global->nsp = 2;
  global->nb  = nb_n0;
  global->rin[0]  = rin[0];
  global->rin[1]  = rin[1];
  global->rin[2]  = rin[2];
  global->rout[0]  = rout[0];
  global->rout[1]  = rout[1];
  global->rout[2]  = rout[2];
  global->vth[0]  = sqrt(2)*vthe;
  global->vth[1]  = sqrt(2)*vthi;
  global->vthb[0]  = sqrt(2)*vtheb;
  global->vthb[1]  = sqrt(2)*vthib;
  global->q[0]  = weight_s;
  global->q[1]  = weight_s;
  global->uf[0]  = udre;
  global->uf[1]  = udri;
  global->sort[0]  = sort_interval;
  global->sort[1]  = sort_interval;
  global->nfac  = nfac;
  global->L_de  = L;

  // Collision model parameters
  global->ee_collisions            = ee_collisions;
  global->ii_collisions            = ii_collisions;
  global->ei_collisions            = ei_collisions;

  global->cvar                     = cvar;
  global->nppc_max                 = nppc_max;
  global->tstep_coll               = tstep_coll;
  global->mi_me                    = mi_me;
  global->Z                        = Z;
  global->wpewce                   = wpe_wce;
  global->nfac                     = nfac;

  //////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;
  //grid->damp = damp;

  // Define the grid
  define_periodic_grid(  0, -0.5*Ly, -0.5*Lz,    // Low corner
                         Lx, 0.5*Ly, 0.5*Lz,     // High corner
                         nx, ny, nz,             // Resolution
                         topology_x, topology_y, topology_z); // Topology

  // ***** Set Field Boundary Conditions *****
  //sim_log("Absorbing fields on X & Z-boundaries");
  //if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), absorb_fields );
  //if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY( 0,0,1), absorb_fields );
 if (global->open_bc_x)
 {
  sim_log("Absorbing fields on X-boundaries");
  if ( ix==0 )
    set_domain_field_bc( BOUNDARY(-1,0,0), absorb_fields );
  if ( ix==topology_x-1 )
    set_domain_field_bc( BOUNDARY( 1,0,0), absorb_fields );
 }

  sim_log("Conducting fields on Z-boundaries");
  if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), pec_fields );
  if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY( 0,0,1), pec_fields );
  //if ( ix==0 )            set_domain_field_bc( BOUNDARY(-1,0,0), pec_fields );
  //if ( ix==topology_x-1 ) set_domain_field_bc( BOUNDARY( 1,0,0), pec_fields );

  // ***** Set Particle Boundary Conditions *****

 if (global->driven_bc_z)
 {
  sim_log("Absorb particles on Z-boundaries");
  if ( iz==0 )
    set_domain_particle_bc( BOUNDARY(0,0,-1), absorb_particles );
  if ( iz==topology_z-1 )
    set_domain_particle_bc( BOUNDARY(0,0,1), absorb_particles );
 } else {
  sim_log("Reflect particles on Z-boundaries");
  if ( iz==0 )
    set_domain_particle_bc( BOUNDARY(0,0,-1), reflect_particles );
  if ( iz==topology_z-1 )
    set_domain_particle_bc( BOUNDARY(0,0,1), reflect_particles );
 }
 if (global->open_bc_x)
 {
  sim_log("Absorb particles on X-boundaries");
  if ( ix==0 )
    set_domain_particle_bc( BOUNDARY(-1,0,0), absorb_particles );
  if ( ix==topology_x-1 )
    set_domain_particle_bc( BOUNDARY(1,0,0), absorb_particles );
 }

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

  // Define resistive layer surrounding boundary --> set thickness=0
  // to eliminate this feature
  double thickness = 0;
#define resistive_layer ((global->open_bc_x && x < hx*thickness) || \
                          (global->open_bc_x && x > Lx-hx*thickness) \
                         || z <-Lz/2+hz*thickness  || z > Lz/2-hz*thickness )

  if (thickness > 0)
  {
    sim_log("Setting resistive layer of thickness "<< thickness);
    set_region_material(resistive_layer, resistive, resistive);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Setup the species

  sim_log("Setting up species. ");
  double nmax = 2.0*Ne/nproc();
  double nmovers = 0.1*nmax;
  double sort_method = 1;   // 0=in place and 1=out of place
  species_t *electron = define_species("electron", -ec, me, nmax, nmovers,
    electron_sort_interval, sort_method);
  species_t *ion = define_species("ion", ec, mi, nmax, nmovers,
    ion_sort_interval, sort_method);

  ///////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log("* Topology:                       " << topology_x
    << " " << topology_y << " " << topology_z);
  sim_log ( "tanhf = " << tanhf );
  sim_log ( "L_di   = " << L_di );
  sim_log ( "rhoi/L   = " << rhoi_L );
  sim_log ( "Ti/Te = " << Ti_Te ) ;
  sim_log ( "nb/n0 = " << nb_n0 ) ;
  sim_log ( "wpe/wce = " << wpe_wce );
  sim_log ( "mi/me = " << mi_me );
  sim_log ( "theta = " << theta );
  sim_log ( "Lpert/Lx = " << Lpert_Lx );
  sim_log ( "dbz/b0 = " << dbz_b0 );
  sim_log ( "taui = " << taui );
  sim_log ( "t_intervali = " << t_intervali );
  sim_log ( "intervali = " << interval );
  sim_log ( "num_step = " << num_step );
  sim_log ( "Lx/di = " << Lx/di );
  sim_log ( "Lx/de = " << Lx/de );
  sim_log ( "Ly/di = " << Ly/di );
  sim_log ( "Ly/de = " << Ly/de );
  sim_log ( "Lz/di = " << Lz/di );
  sim_log ( "Lz/de = " << Lz/de );
  sim_log ( "nx = " << nx );
  sim_log ( "ny = " << ny );
  sim_log ( "nz = " << nz );
  sim_log ( "damp = " << damp );
  sim_log ( "courant = " << c*dt/dg );
  sim_log ( "nproc = " << nproc ()  );
  sim_log ( "nppc = " << nppc );
  sim_log ( "b0 = " << b0 );
  sim_log ( "v_A (based on nb) = " << v_A );
  sim_log ( "di = " << di );
  sim_log ( "Ne = " << Ne );
  sim_log ( "Ne_sheet = " << Ne_sheet );
  sim_log ( "Ne_back = " << Ne_back );
  sim_log ( "total # of particles = " << 2*Ne );
  sim_log ( "dt*wpe = " << wpe*dt );
  sim_log ( "dt*wce = " << wce*dt );
  sim_log ( "dt*wci = " << wci*dt );
  sim_log ( "energies_interval: " << energies_interval );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/rhoi = " << (Lx/nx)/(vthi/wci)  );
  sim_log ( "dx/rhoe = " << (Lx/nx)/(vthe/wce)  );
  sim_log ( "L/debye = " << L/(vthe/wpe)  );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "n0 = " << n0 );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );
  sim_log ( "vdri/c = " << vdri/c );
  sim_log ( "vdre/c = " << vdre/c );
  sim_log ("* nu/wce:                        " << nuei_wce);
  sim_log ("* nu*dt_coll:                    " << nuei_wce/wpe_wce*dt_coll);
  sim_log ( "ee_collisions   = " << ee_collisions );
  sim_log ( "ei_collisions   = " << ei_collisions );
  sim_log ( "ii_collisions   = " << ii_collisions );
  sim_log ( "Open BC in x?   = " << open_bc_x );
  sim_log ( "Driven BC in z? = " << driven_bc_z );

  // Dump simulation information to file "info"
  if (rank() == 0 ) {
    FileIO fp_info;
    if ( ! (fp_info.open("info", io_write)==ok) ) ERROR(("Cannot open file."));
    fp_info.print("           ***** Simulation parameters ***** \n");
    fp_info.print("              L/di   =               %e\n", L_di);
    fp_info.print("              L/de   =               %e\n", L/de);
    fp_info.print("              rhoi/L =               %e\n", rhoi_L);
    fp_info.print("              Ti/Te  =               %e\n", Ti_Te );
    fp_info.print("              Tbi/Ti =               %e\n", Tbi_Ti );
    fp_info.print("              Tbe/Te =               %e\n", Tbe_Te );
    fp_info.print("              nb/n0 =                %e\n", nb_n0 );
    fp_info.print("              wpe/wce =              %e\n", wpe_wce );
    fp_info.print("              mi/me =                %e\n", mi_me );
    fp_info.print("              theta =                %e\n", theta );
    fp_info.print("              taui =                 %e\n", taui );
    fp_info.print("              num_step =             %i\n", num_step );
    fp_info.print("              Lx/de =                %e\n", Lx/de );
    fp_info.print("              Ly/de =                %e\n", Ly/de );
    fp_info.print("              Lz/de =                %e\n", Lz/de );
    fp_info.print("              Lx/di =                %e\n", Lx/di );
    fp_info.print("              Ly/di =                %e\n", Ly/di );
    fp_info.print("              Lz/di =                %e\n", Lz/di );
    fp_info.print("              nx =                   %e\n", nx );
    fp_info.print("              ny =                   %e\n", ny );
    fp_info.print("              nz =                   %e\n", nz );
    fp_info.print("              damp =                 %e\n", damp );
    fp_info.print("              courant =              %e\n", c*dt/dg );
    fp_info.print("              nproc =                %i\n", nproc() );
    fp_info.print("              nppc =                 %e\n", nppc );
    fp_info.print("              b0 =                   %e\n", b0 );
    fp_info.print("              v_A (based on nb) =    %e\n", v_A );
    fp_info.print("              di =                   %e\n", di );
    fp_info.print("              Ne =                   %e\n", Ne );
    fp_info.print("              Ne_sheet =             %e\n", Ne_sheet );
    fp_info.print("              Ne_back =              %e\n", Ne_back );
    fp_info.print("              total # of particles = %e\n", 2*Ne );
    fp_info.print("              dt*wpe =               %e\n", wpe*dt );
    fp_info.print("              dt*wce =               %e\n", wce*dt );
    fp_info.print("              dt*wci =               %e\n", wci*dt );
    fp_info.print("              energies_interval:     %i\n",
      energies_interval);
    fp_info.print("              dx/de =                %e\n", Lx/(de*nx) );
    fp_info.print("              dy/de =                %e\n", Ly/(de*ny) );
    fp_info.print("              dz/de =                %e\n", Lz/(de*nz) );
    fp_info.print("              L/debye =              %e\n", L/(vthe/wpe) );
    fp_info.print("              dx/rhoi =              %e\n",
      (Lx/nx)/(vthi/wci) );
    fp_info.print("              dx/rhoe =              %e\n",
      (Lx/nx)/(vthe/wce) );
    fp_info.print("              dx/debye =             %e\n",
      (Lx/nx)/(vthe/wpe) );
    fp_info.print("              n0 =                   %e\n", n0 );
    fp_info.print("              vthi/c =               %e\n", vthi/c );
    fp_info.print("              vthe/c =               %e\n", vthe/c );
    fp_info.print("              vdri/c =               %e\n", vdri/c );
    fp_info.print("              vdre/c =               %e\n", vdre/c );
    fp_info.print(" tstep_coll:                    %i\n", tstep_coll);
    fp_info.print(" nu/wce:                        %g\n", nuei_wce);
    fp_info.print(" nu*dt_coll:                    %g\n",
      nuei_wce/wpe_wce*dt_coll);
    fp_info.print("              ***************************\n");
    fp_info.close();
  }

  ////////////////////////////
  // Load fields

  sim_log( "Loading fields" );
  set_region_field( everywhere, 0, 0, 0,                    // Electric field
    cs*b0*tanh(z/L)+dbx*cos(2.0*M_PI*(x-0.5*Lx)/Lpert)*sin(M_PI*z/Lz), //Bx
    -sn*b0*tanh(z/L) + b0*bg, //By
    dbz*cos(M_PI*z/Lz)*sin(2.0*M_PI*(x-0.5*Lx)/Lpert) ); // Bz

  // Localized Perturbation to lauch a light wave

//    # define R2 ((x-0.5*Lx)*(x-0.5*Lx) + z*z)/(L*L)
//    # define PERT  0.2*tanh(R2)/cosh(R2)
//      set_region_field( everywhere, 0, 0, 0,                    // Electric field
//        cs*b0*tanh(z/L) - z*PERT, //Bx
//        b0*bg, //By
//        (x-0.5*Lx)*PERT ); // Bz

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

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

  global->fdParams.format = band;
  sim_log ( "Fields output format = band" );

  global->hedParams.format = band;
  sim_log ( "Electron species output format = band" );

  global->hHdParams.format = band;
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
  sprintf(global->fdParams.baseDir, "fields");

  // base file name for fields output
  sprintf(global->fdParams.baseFileName, "fields");

  global->fdParams.stride_x = 1;
  global->fdParams.stride_y = 1;
  global->fdParams.stride_z = 1;

  // add field parameters to list
  global->outputParams.push_back(&global->fdParams);

  sim_log ( "Fields x-stride " << global->fdParams.stride_x );
  sim_log ( "Fields y-stride " << global->fdParams.stride_y );
  sim_log ( "Fields z-stride " << global->fdParams.stride_z );

  // relative path to electron species data from global header
  sprintf(global->hedParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(global->hedParams.baseFileName, "ehydro");

  global->hedParams.stride_x = 1;
  global->hedParams.stride_y = 1;
  global->hedParams.stride_z = 1;

  // add electron species parameters to list
  global->outputParams.push_back(&global->hedParams);

  sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
  sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
  sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

  // relative path to electron species data from global header
  sprintf(global->hHdParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(global->hHdParams.baseFileName, "Hhydro");

  global->hHdParams.stride_x = 1;
  global->hHdParams.stride_y = 1;
  global->hHdParams.stride_z = 1;

  sim_log ( "Ion species x-stride " << global->hHdParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hHdParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hHdParams.stride_z );

  // add electron species parameters to list
  global->outputParams.push_back(&global->hHdParams);

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

  global->fdParams.output_variables( all );
  global->hedParams.output_variables( all );
  global->hHdParams.output_variables( all );

  /*--------------------------------------------------------------------------
   * Convenience functions for simlog output
   *------------------------------------------------------------------------*/

  char varlist[512];
  create_field_list(varlist, global->fdParams);

  sim_log ( "Fields variable list: " << varlist );

  create_hydro_list(varlist, global->hedParams);

  sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, global->hHdParams);

  sim_log ( "Ion species variable list: " << varlist );

  /* ---------------------------------------------
     Add parameters for the energy diagnostics
     --------------------------------------------- */

  global->ede.sp_id = electron->id;
  global->ede.vth = sqrt(2.0)*vthe;
  sprintf(global->ede.fname,global->hedParams.baseFileName);

  global->edi.sp_id = ion->id;
  global->edi.vth = sqrt(2.0)*vthi;
  sprintf(global->edi.fname, global->hHdParams.baseFileName);

  global->nex  = 6;
  global->emax = 120;

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

#define should_dump(x)                                                  \
  (global->x##_interval>0 && remainder(step(), global->x##_interval) == 0)

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

  const int nsp=global->nsp;
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;

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
    global_header("global", global->outputParams);
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

  if(step() == -1 || should_dump(fields)) field_dump(global->fdParams);

  /*--------------------------------------------------------------------------
   * Electron species output
   *------------------------------------------------------------------------*/

  if(should_dump(ehydro)) hydro_dump("electron", global->hedParams);

  /*--------------------------------------------------------------------------
   * Ion species output
   *------------------------------------------------------------------------*/

  if(should_dump(Hhydro)) hydro_dump("ion", global->hHdParams);

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  if(step() && !(step()%global->restart_interval)) {
    if(!global->rtoggle) {
      global->rtoggle = 1;
      checkpt("restart1/restart", 0);
    }
    else {
      global->rtoggle = 0;
      checkpt("restart2/restart", 0);
    } // if
  } // if

  // Dump particle data

  char subdir[36];
  if ( should_dump(eparticle) && step() !=0
       && step() > 56*(global->fields_interval)  ) {
    // if ( should_dump(eparticle) && step() !=0 ) {
    sprintf(subdir,"particle/T.%d",step());
    dump_mkdir(subdir);
    sprintf(subdir,"particle/T.%d/eparticle",step());
    dump_particles("electron", subdir);
  }

  if ( should_dump(Hparticle) && step() !=0
       && step() > 56*(global->fields_interval)  ) {
    sprintf(subdir,"particle/T.%d/Hparticle",step());
    dump_particles("ion", subdir);
  }

  // Shut down simulation when wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will
  // be synchronized across processors. Note that this is only checked every
  // few timesteps to eliminate the expensive mp_elapsed call from every
  // timestep. mp_elapsed has an ALL_REDUCE in it!

  if (( step()>0 && global->quota_check_interval>0
        && (step()&global->quota_check_interval)==0)
      || (global->write_end_restart) ) {
    if ( (global->write_end_restart) ) {

      global->write_end_restart = 0; // reset restart flag

      sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
      double dumpstart = uptime();

      checkpt("restart0/restart",0);

      mp_barrier(  ); // Just to be safe
      sim_log( "Restart dump restart completed." );
      double dumpelapsed = uptime() - dumpstart;
      sim_log("Restart duration "<< dumpelapsed);
      exit(0); // Exit or abort?
    }
    if( uptime() > global->quota_sec ) global->write_end_restart = 1;
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
