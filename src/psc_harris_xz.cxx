
#include <psc.h>
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
#include "balance.hxx"
#include "fields3d.hxx"
#include "setup_particles.hxx"

#include <psc_particles_as_single.h>
#include <psc_particles_vpic.h>

#include "rngpool_iface.h"

#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static RngPool *rngpool; // FIXME, should be member (of struct psc, really)

// ----------------------------------------------------------------------

#define VAR(x) (void *)offsetof(struct psc_harris, x)
static struct param psc_harris_descr[] = {
  { "wpedt_max"             , VAR(prm.wpedt_max)             , PARAM_DOUBLE(.36)  },
  { "wpe_wce"               , VAR(prm.wpe_wce)               , PARAM_DOUBLE(2.),
    .help = "electron plasma freq / electron cyclotron freq" },
  { "mi_me"                 , VAR(prm.mi_me)                 , PARAM_DOUBLE(25.),
    .help = "ion mass over electron mass" },
  { "Lx_di"                 , VAR(prm.Lx_di)                 , PARAM_DOUBLE(25.6),
    .help = "x-size of simulatin domain in terms of d_i" },
  { "Ly_di"                 , VAR(prm.Ly_di)                 , PARAM_DOUBLE(1.),
    .help = "y-size of simulatin domain in terms of d_i" },
  { "Lz_di"                 , VAR(prm.Lz_di)                 , PARAM_DOUBLE(12.8),
    .help = "z-size of simulatin domain in terms of d_i" },
  { "nppc"                   , VAR(prm.nppc)                 , PARAM_DOUBLE(100.),
    .help = "average number of macro particle per cell per species" },

  { "ion_sort_interval"     , VAR(prm.ion_sort_interval)     , PARAM_INT(1000)    },
  { "electron_sort_interval", VAR(prm.electron_sort_interval), PARAM_INT(1000)    },

  { "taui"                  , VAR(prm.taui)                  , PARAM_DOUBLE(100.),
    .help = "simulation wci's to run" },
  { "t_intervali"           , VAR(prm.t_intervali)           , PARAM_DOUBLE(1.),
    .help = "output interval in terms of 1/wci" },
  { "output_field_interval" , VAR(prm.output_field_interval) , PARAM_DOUBLE(1.),
    .help = "field output interval in terms of 1/wci" },
  { "output_particle_interval", VAR(prm.output_particle_interval), PARAM_DOUBLE(0.),
    .help = "particle output interval in terms of 1/wci" },

  { "L_di"                  , VAR(prm.L_di)                  , PARAM_DOUBLE(.5),
    .help = "Sheet thickness / ion inertial length" },
  { "Ti_Te"                 , VAR(prm.Ti_Te)                 , PARAM_DOUBLE(5.),
    .help = "Ion temperature / electron temperature" },
  { "nb_n0"                 , VAR(prm.nb_n0)                 , PARAM_DOUBLE(.228),
    .help = "background plasma density" },
  { "Tbe_Te"                , VAR(prm.Tbe_Te)                , PARAM_DOUBLE(.7598),
    .help = "Ratio of background T_e to Harris T_e" },
  { "Tbi_Ti"                , VAR(prm.Tbi_Ti)                , PARAM_DOUBLE(.3039),
    .help = "Ratio of background T_i to Harris T_i" },
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
// psc_harris_create

static void
psc_harris_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nicell = 1;
  psc->prm.cfl = 0.99;

  psc->prm.stats_every = 100;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_CONDUCTING_WALL;
  psc->domain.bnd_fld_hi[2] = BND_FLD_CONDUCTING_WALL;
 
  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_REFLECTING;
  psc->domain.bnd_part_hi[2] = BND_PART_REFLECTING;

  psc_method_set_type(psc->method, "vpic");

  psc_sort_set_type(psc->sort, "vpic");
  // FIXME: the "vpic" sort actually keeps track of per-species sorting intervals
  // internally
  psc_sort_set_param_int(psc->sort, "every", 1);

  psc_collision_set_type(psc->collision, "vpic");

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc->push_particles, "vpic");
  psc_bnd_particles_set_type(psc->bnd_particles, "vpic");

  psc_push_fields_set_type(psc->push_fields, "vpic");
  psc_bnd_set_type(psc->bnd, "vpic");
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "vpic");

  psc_marder_set_type(psc->marder, "vpic");
  // FIXME, marder "vpic" manages its own cleaning intervals
  psc_marder_set_param_int(psc->marder, "every_step", 1);
  psc_marder_set_param_int(psc->marder, "clean_div_e_interval", 50);
  psc_marder_set_param_int(psc->marder, "clean_div_b_interval", 50);
  psc_marder_set_param_int(psc->marder, "sync_shared_interval", 50);
  psc_marder_set_param_int(psc->marder, "num_div_e_round", 2);
  psc_marder_set_param_int(psc->marder, "num_div_b_round", 2);
}

// FIXME, helper should go somewhere...

static inline double trunc_granular( double a, double b )
{
  return b * (int)(a/b);
}

// ----------------------------------------------------------------------
// psc_harris_setup_ic

static void
psc_harris_setup_ic(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  struct vpic_harris_params *prm = &sub->prm;

  sub->n_global_patches = psc->domain.np[0] * psc->domain.np[1] * psc->domain.np[2];
  
  assert(psc->domain.np[2] <= 2); // For load balance, keep "1" or "2" for Harris sheet

  // FIXME, the general normalization stuff should be shared somehow

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
  int *gdims       = psc->domain.gdims;
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

  psc->domain.length[0] = phys->Lx;
  psc->domain.length[1] = phys->Ly;
  psc->domain.length[2] = phys->Lz;

  psc->domain.corner[0] = 0.;
  psc->domain.corner[1] = -.5 * phys->Ly;
  psc->domain.corner[2] = -.5 * phys->Lz;
}

// ----------------------------------------------------------------------
// psc_harris_setup_domain

static void
psc_harris_setup_domain(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  MPI_Comm comm = psc_comm(psc);

  psc_setup_coeff(psc); // FIXME -- in the middle of things here, will be done again later
  psc_setup_domain(psc);
  
  // Setup basic grid parameters
  double dx[3], xl[3], xh[3];
  for (int d = 0; d < 3; d++) {
    dx[d] = psc->domain.length[d] / psc->domain.gdims[d];
    xl[d] = psc->domain.corner[d];
    xh[d] = xl[d] + psc->domain.length[d];
  }
  Simulation_setup_grid(sub->sim, dx, phys->dt, phys->c, phys->eps0);

  // Define the grid
  Simulation_define_periodic_grid(sub->sim, xl, xh, psc->domain.gdims, psc->domain.np);

  int p = 0;
  bool left = psc_at_boundary_lo(psc, p, 0);
  bool right = psc_at_boundary_hi(psc, p, 0);

  bool bottom = psc_at_boundary_lo(psc, p, 2);
  bool top = psc_at_boundary_hi(psc, p, 2);

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
    if (bottom) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0,-1), BND_PART_ABSORBING);
    if (top   ) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0, 1), BND_PART_ABSORBING);
  } else {
    mpi_printf(comm, "Reflect particles on Z-boundaries\n");
    if (bottom) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0,-1), BND_PART_REFLECTING);
    if (top   ) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(0,0, 1), BND_PART_REFLECTING);
  }
  if (sub->prm.open_bc_x) {
    mpi_printf(comm, "Absorb particles on X-boundaries\n");
    if (left)  Simulation_set_domain_particle_bc(sub->sim, BOUNDARY(-1,0,0), BND_PART_ABSORBING);
    if (right) Simulation_set_domain_particle_bc(sub->sim, BOUNDARY( 1,0,0), BND_PART_ABSORBING);
  }
}

// ----------------------------------------------------------------------
// psc_harris_setup_fields

static void
psc_harris_setup_fields(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  MPI_Comm comm = psc_comm(psc);

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
  assert(psc->nr_patches > 0);
  struct globals_physics *phys = &sub->phys;
  Simulation_set_region_resistive_harris(sub->sim, &sub->prm, phys, psc->patch[0].dx,
					 0., resistive);
#endif
}

// ----------------------------------------------------------------------
// psc_harris_setup_species

static void
psc_harris_setup_species(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  MPI_Comm comm = psc_comm(psc);

  mpi_printf(comm, "Setting up species.\n");
  double nmax = sub->prm.overalloc * phys->Ne / sub->n_global_patches;
  double nmovers = .1 * nmax;
  double sort_method = 1;   // 0=in place and 1=out of place

  psc_set_kinds(psc, {{-phys->ec, phys->me, "e"}, {phys->ec, phys->mi, "i"}});

  Simulation_define_species(sub->sim, "electron", -phys->ec, phys->me, nmax, nmovers,
			    sub->prm.electron_sort_interval, sort_method);
  Simulation_define_species(sub->sim, "ion", phys->ec, phys->mi, nmax, nmovers,
			    sub->prm.ion_sort_interval, sort_method);
}

// ----------------------------------------------------------------------
// psc_harris_setup_log

static void
psc_harris_setup_log(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  MPI_Comm comm = psc_comm(psc);

  mpi_printf(comm, "***********************************************\n");
  mpi_printf(comm, "* Topology: %d x %d x %d\n",
	     psc->domain.np[0], psc->domain.np[1], psc->domain.np[2]);
  mpi_printf(comm, "tanhf    = %g\n", phys->tanhf);
  mpi_printf(comm, "L_di     = %g\n", sub->prm.L_di);
  mpi_printf(comm, "rhoi/L   = %g\n", phys->rhoi_L);
  mpi_printf(comm, "Ti/Te    = %g\n", sub->prm.Ti_Te) ;
  mpi_printf(comm, "nb/n0    = %g\n", sub->prm.nb_n0) ;
  mpi_printf(comm, "wpe/wce  = %g\n", sub->prm.wpe_wce);
  mpi_printf(comm, "mi/me    = %g\n", sub->prm.mi_me);
  mpi_printf(comm, "theta    = %g\n", sub->prm.theta);
  mpi_printf(comm, "Lpert/Lx = %g\n", sub->prm.Lpert_Lx);
  mpi_printf(comm, "dbz/b0   = %g\n", sub->prm.dbz_b0);
  mpi_printf(comm, "taui     = %g\n", sub->prm.taui);
  mpi_printf(comm, "t_intervali = %g\n", sub->prm.t_intervali);
  mpi_printf(comm, "num_step = %d\n", psc->prm.nmax);
  mpi_printf(comm, "Lx/di = %g\n", phys->Lx/phys->di);
  mpi_printf(comm, "Lx/de = %g\n", phys->Lx/phys->de);
  mpi_printf(comm, "Ly/di = %g\n", phys->Ly/phys->di);
  mpi_printf(comm, "Ly/de = %g\n", phys->Ly/phys->de);
  mpi_printf(comm, "Lz/di = %g\n", phys->Lz/phys->di);
  mpi_printf(comm, "Lz/de = %g\n", phys->Lz/phys->de);
  mpi_printf(comm, "nx = %d\n", psc->domain.gdims[0]);
  mpi_printf(comm, "ny = %d\n", psc->domain.gdims[1]);
  mpi_printf(comm, "nz = %d\n", psc->domain.gdims[2]);
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
  mpi_printf(comm, "dx/de = %g\n", phys->Lx/(phys->de*psc->domain.gdims[0]));
  mpi_printf(comm, "dy/de = %g\n", phys->Ly/(phys->de*psc->domain.gdims[1]));
  mpi_printf(comm, "dz/de = %g\n", phys->Lz/(phys->de*psc->domain.gdims[2]));
  mpi_printf(comm, "dx/rhoi = %g\n", (phys->Lx/psc->domain.gdims[0])/(phys->vthi/phys->wci));
  mpi_printf(comm, "dx/rhoe = %g\n", (phys->Lx/psc->domain.gdims[0])/(phys->vthe/phys->wce));
  mpi_printf(comm, "L/debye = %g\n", phys->L/(phys->vthe/phys->wpe));
  mpi_printf(comm, "dx/debye = %g\n", (phys->Lx/psc->domain.gdims[0])/(phys->vthe/phys->wpe));
  mpi_printf(comm, "n0 = %g\n", phys->n0);
  mpi_printf(comm, "vthi/c = %g\n", phys->vthi/phys->c);
  mpi_printf(comm, "vthe/c = %g\n", phys->vthe/phys->c);
  mpi_printf(comm, "vdri/c = %g\n", phys->vdri/phys->c);
  mpi_printf(comm, "vdre/c = %g\n", phys->vdre/phys->c);
  mpi_printf(comm, "Open BC in x?   = %d\n", sub->prm.open_bc_x);
  mpi_printf(comm, "Driven BC in z? = %d\n", sub->prm.driven_bc_z);
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
 
// ----------------------------------------------------------------------
// psc_harris_setup

static void
psc_harris_setup(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  MPI_Comm comm = psc_comm(psc);

  psc_harris_setup_ic(psc);

  // Determine the time step
  phys->dg = courant_length(psc->domain.length, psc->domain.gdims);
  phys->dt = psc->prm.cfl * phys->dg / phys->c; // courant limited time step
  if (phys->wpe * phys->dt > sub->prm.wpedt_max) {
    phys->dt = sub->prm.wpedt_max / phys->wpe;  // override timestep if plasma frequency limited
  }
  psc->dt = phys->dt;

  psc->prm.nmax = (int) (sub->prm.taui / (phys->wci*phys->dt)); // number of steps from taui

  if (strcmp(psc_method_type(psc->method), "vpic") != 0) {
    psc_setup_super(psc);
    psc_harris_setup_log(psc);
    return;
  }
  bool split;
  psc_method_get_param_bool(psc->method, "split", &split);
  if (!split) {
    psc_setup_super(psc);
    psc_harris_setup_log(psc);
    return;
  }

  sub->sim = Simulation_create();
  psc_method_set_param_ptr(psc->method, "sim", sub->sim);
  // set high level VPIC simulation parameters
  // FIXME, will be unneeded eventually
  Simulation_set_params(sub->sim, psc->prm.nmax, psc->prm.stats_every,
			psc->prm.stats_every / 2, psc->prm.stats_every / 2,
			psc->prm.stats_every / 2);
  psc_harris_setup_domain(psc);
  psc_harris_setup_fields(psc);
  psc_harris_setup_species(psc);
  psc_harris_setup_log(psc);

  int interval = (int) (sub->prm.t_intervali / (phys->wci*phys->dt));
  Simulation_diagnostics_init(sub->sim, interval);

  psc->n_state_fields = VPIC_MFIELDS_N_COMP;
  psc->ibn[0] = psc->ibn[1] = psc->ibn[2] = 1;

  // partition and initial balancing
  auto n_prts_by_patch_old = psc_method_setup_partition(psc->method, psc);
  auto balance = PscBalanceBase{psc->balance};
  auto n_prts_by_patch_new = balance.initial(psc, n_prts_by_patch_old);

  psc->particles = PscMparticlesCreate(mrc_domain_comm(psc->mrc_domain), psc->grid(),
				       psc->prm.particles_base).mprts();

  psc->flds = PscMfieldsCreate(mrc_domain_comm(psc->mrc_domain), psc->grid(),
			       psc->n_state_fields, psc->ibn, psc->prm.fields_base).mflds();

  SetupParticles<mparticles_t::sub_t>::setup_particles(psc, n_prts_by_patch_new);

  psc_set_ic_fields(psc);
  
  Simulation_diagnostics_setup(sub->sim);

  mpi_printf(comm, "*** Finished with user-specified initialization ***\n");
  
  psc_setup_member_objs(psc);

  if (sub->prm.output_field_interval > 0) {
    struct psc_output_fields *out;
    mrc_obj_for_each_child(out, psc->output_fields_collection, struct psc_output_fields) {
      psc_output_fields_set_param_int(out, "pfield_step",
				      (int) (sub->prm.output_field_interval / (phys->wci*phys->dt)));
    }
  }

  if (sub->prm.output_particle_interval > 0) {
    psc_output_particles_set_param_int(psc->output_particles, "every_step",
				      (int) (sub->prm.output_particle_interval / (phys->wci*phys->dt)));
  }
}

// ----------------------------------------------------------------------
// psc_harris_init_field

static double
psc_harris_init_field(struct psc *psc, double crd[3], int m)
{
  struct psc_harris *sub = psc_harris(psc);
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
    //inject_particle(electron, &prt, 0., 0);

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
    create           = psc_harris_create;
    setup            = psc_harris_setup;
    destroy          = psc_harris_destroy;
    read             = psc_harris_read;
    init_field       = psc_harris_init_field;
    setup_particles  = psc_harris_setup_particles;
  }
} psc_harris_ops;

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_harris_ops);
}
