
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

#include <psc_particles_as_single.h>
#include <psc_particles_vpic.h>

#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>

#include "libpsc/vpic/vpic_iface.h"

static void psc_harris_set_ic_particles(struct psc *psc);

// ----------------------------------------------------------------------

void *
rng(int i)
{
  return NULL;
}

double
uniform(void *dummy, double lo, double hi)
{
  return lo + (hi - lo) * random() / ((float) RAND_MAX + 1);
}

double
normal(void *dummy, double mu, double sigma)
{
  float ran1, ran2;
  do {
    ran1 = random() / ((float) RAND_MAX + 1);
    ran2 = random() / ((float) RAND_MAX + 1);
  } while (ran1 >= 1.f || ran2 >= 1.f);
	      
  return mu + sigma * sqrtf(-2.f * logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
}

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
  vpic_simulation_setup_grid(dx, phys->dt, phys->c, phys->eps0);

  // Define the grid
  vpic_simulation_define_periodic_grid(xl, xh, psc->domain.gdims, psc->domain.np);

  int p = 0;
  bool left = psc_at_boundary_lo(psc, p, 0);
  bool right = psc_at_boundary_hi(psc, p, 0);

  bool bottom = psc_at_boundary_lo(psc, p, 2);
  bool top = psc_at_boundary_hi(psc, p, 2);

  // ***** Set Field Boundary Conditions *****
  if (sub->prm.open_bc_x) {
    mpi_printf(comm, "Absorbing fields on X-boundaries\n");
    if (left ) vpic_simulation_set_domain_field_bc(BOUNDARY(-1,0,0), BND_FLD_ABSORBING);
    if (right) vpic_simulation_set_domain_field_bc(BOUNDARY( 1,0,0), BND_FLD_ABSORBING);
  }
  
  mpi_printf(comm, "Conducting fields on Z-boundaries\n");
  if (bottom) vpic_simulation_set_domain_field_bc(BOUNDARY(0,0,-1), BND_FLD_CONDUCTING_WALL);
  if (top   ) vpic_simulation_set_domain_field_bc(BOUNDARY(0,0, 1), BND_FLD_CONDUCTING_WALL);

  // ***** Set Particle Boundary Conditions *****
  if (sub->prm.driven_bc_z) {
    mpi_printf(comm, "Absorb particles on Z-boundaries\n");
    if (bottom) vpic_simulation_set_domain_particle_bc(BOUNDARY(0,0,-1), BND_PART_ABSORBING);
    if (top   ) vpic_simulation_set_domain_particle_bc(BOUNDARY(0,0, 1), BND_PART_ABSORBING);
  } else {
    mpi_printf(comm, "Reflect particles on Z-boundaries\n");
    if (bottom) vpic_simulation_set_domain_particle_bc(BOUNDARY(0,0,-1), BND_PART_REFLECTING);
    if (top   ) vpic_simulation_set_domain_particle_bc(BOUNDARY(0,0, 1), BND_PART_REFLECTING);
  }
  if (sub->prm.open_bc_x) {
    mpi_printf(comm, "Absorb particles on X-boundaries\n");
    if (left)  vpic_simulation_set_domain_particle_bc(BOUNDARY(-1,0,0), BND_PART_ABSORBING);
    if (right) vpic_simulation_set_domain_particle_bc(BOUNDARY( 1,0,0), BND_PART_ABSORBING);
  }
}

// ----------------------------------------------------------------------
// psc_harris_setup_fields

static void
psc_harris_setup_fields(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  MPI_Comm comm = psc_comm(psc);

  mpi_printf(comm, "Setting up materials.\n");

  vpic_simulation_define_material("vacuum", 1., 1., 0., 0.);
  struct material *resistive =
    vpic_simulation_define_material("resistive", 1., 1., 1., 0.);

  vpic_simulation_define_field_array(NULL, 0.);

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

  //////////////////////////////////////////////////////////////////////////////
  // Finalize Field Advance

  mpi_printf(comm, "Finalizing Field Advance\n");

  assert(psc->nr_patches > 0);
  vpic_simulation_set_region_resistive_harris(&sub->prm, phys, psc->patch[0].dx,
					      0., resistive);
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
  double nmax = 2.*phys->Ne / sub->n_global_patches;
  double nmovers = .1*nmax;
  double sort_method = 1;   // 0=in place and 1=out of place

  struct psc_kind kinds[2] = {
    [KIND_ELECTRON] = { .name = "e", .q = -phys->ec, .m = phys->me, },
    [KIND_ION     ] = { .name = "i", .q =  phys->ec, .m = phys->mi, },
  };
  psc_set_kinds(psc, 2, kinds);

  vpic_simulation_define_species("electron", -phys->ec, phys->me, nmax, nmovers,
				 sub->prm.electron_sort_interval, sort_method);
  vpic_simulation_define_species("ion", phys->ec, phys->mi, nmax, nmovers,
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
// psc_setup_base_mprts
//
// FIXME, duplicated

static void
psc_setup_base_mprts(struct psc *psc)
{
  psc->particles = psc_mparticles_create(mrc_domain_comm(psc->mrc_domain));
  psc_mparticles_set_type(psc->particles, psc->prm.particles_base);
  psc_mparticles_set_name(psc->particles, "mparticles");
  int nr_patches;
  mrc_domain_get_patches(psc->mrc_domain, &nr_patches);
  psc_mparticles_set_param_int(psc->particles, "nr_patches", nr_patches);
  if (psc->prm.particles_base_flags == 0) {
    psc->prm.particles_base_flags = psc_push_particles_get_mp_flags(ppsc->push_particles);
  }
  psc_mparticles_set_param_int(psc->particles, "flags", psc->prm.particles_base_flags);
  psc_mparticles_setup(psc->particles);
}

// ----------------------------------------------------------------------
// psc_setup_base_mflds
//
// FIXME, duplicated

static void
psc_setup_base_mflds(struct psc *psc)
{
  psc->flds = psc_mfields_create(mrc_domain_comm(psc->mrc_domain));
  psc_mfields_list_add(&psc_mfields_base_list, &psc->flds);
  psc_mfields_set_type(psc->flds, psc->prm.fields_base);
  psc_mfields_set_name(psc->flds, "mfields");
  psc_mfields_set_param_obj(psc->flds, "domain", psc->mrc_domain);
  psc_mfields_set_param_int(psc->flds, "nr_fields", psc->n_state_fields);
  psc_mfields_set_param_int3(psc->flds, "ibn", psc->ibn);
  psc_mfields_setup(psc->flds);
}


// ----------------------------------------------------------------------
// courant length
//
// FIXME, the dt calculating should be consolidated with what regular PSC does

static inline double
courant_length(double length[3], int gdims[3])
{
  double inv_sum = 0.;
  for (int d = 0; d > 3; d++) {
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

  bool split;
  psc_method_get_param_bool(psc->method, "split", &split);
  if (split) {
  psc_harris_setup_ic(psc);

  // Determine the time step
  phys->dg = courant_length(psc->domain.length, psc->domain.gdims);
  phys->dt = psc->prm.cfl * phys->dg / phys->c; // courant limited time step
  if (phys->wpe * phys->dt > sub->prm.wpedt_max) {
    phys->dt = sub->prm.wpedt_max / phys->wpe;  // override timestep if plasma frequency limited
  }
  psc->dt = phys->dt;

  psc->prm.nmax = (int) (sub->prm.taui / (phys->wci*phys->dt)); // number of steps from taui

  psc_marder_set_param_int(psc->marder, "clean_div_e_interval", psc->prm.stats_every / 2);
  psc_marder_set_param_int(psc->marder, "clean_div_b_interval", psc->prm.stats_every / 2);
  psc_marder_set_param_int(psc->marder, "sync_shared_interval", psc->prm.stats_every / 2);
  psc_marder_set_param_int(psc->marder, "num_div_e_round", 2);
  psc_marder_set_param_int(psc->marder, "num_div_b_round", 2);
    
  // set high level VPIC simulation parameters
  // FIXME, will be unneeded eventually
  vpic_simulation_new();
  vpic_simulation_set_params(psc->prm.nmax,
			     psc->prm.stats_every,
			     psc->prm.stats_every / 2,
			     psc->prm.stats_every / 2,
			     psc->prm.stats_every / 2);
  psc_harris_setup_domain(psc);
  psc_harris_setup_fields(psc);
  psc_harris_setup_species(psc);
  psc_harris_setup_log(psc);
  }

  int interval = (int) (sub->prm.t_intervali / (phys->wci*phys->dt));
  vpic_diagnostics_init(interval);
  
  psc->n_state_fields = VPIC_MFIELDS_N_COMP;
  psc->ibn[0] = psc->ibn[1] = psc->ibn[2] = 1;

  // partition and initial balancing
  int *n_prts_by_patch = calloc(psc->nr_patches, sizeof(*n_prts_by_patch));
  psc_method_setup_partition(psc->method, psc, n_prts_by_patch);
  psc_balance_initial(psc->balance, psc, &n_prts_by_patch);
  free(n_prts_by_patch);

  psc_setup_base_mprts(psc);
  psc_setup_base_mflds(psc);

  psc_harris_set_ic_particles(psc);

  psc_set_ic_fields(psc);
  
  vpic_diagnostics_setup();

  mpi_printf(comm, "*** Finished with user-specified initialization ***\n");
  
  psc_setup_member_objs(psc);
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
// psc_harris_set_ic_particles
//
// set particles x^{n+1/2}, p^{n+1/2}

typedef struct psc_particle_inject {
  double x[3];
  double u[3];
  double w;
  int kind;
} particle_inject_t;

void
inject_particle(struct species *sp, const struct psc_particle_inject *prt,
		double age, int update_rhob)
{
  vpic_simulation_inject_particle(sp, prt->x[0], prt->x[1], prt->x[2],
				  prt->u[0], prt->u[1], prt->u[2], prt->w, age, update_rhob);
}

static void
psc_harris_set_ic_particles(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;
  MPI_Comm comm = psc_comm(psc);

  struct psc_mparticles *mprts_vpic = psc_mparticles_get_as(psc->particles, "vpic", MP_DONT_COPY | MP_DONT_RESIZE);
  
  double cs = cos(sub->prm.theta), sn = sin(sub->prm.theta);
  double Ne_sheet = phys->Ne_sheet, vthe = phys->vthe, vthi = phys->vthi;
  double weight_s = phys->weight_s;
  double tanhf = phys->tanhf, L = phys->L;
  double gdre = phys->gdre, udre = phys->udre, gdri = phys->gdri, udri = phys->udri;
  double Ne_back = phys->Ne_back, vtheb = phys->vtheb, vthib = phys->vthib;
  double weight_b = phys->weight_b;
  //grid_t *grid = simulation->grid;
  int n_global_patches = sub->n_global_patches;
  
  // LOAD PARTICLES

  mpi_printf(comm, "Loading particles\n");

  // Do a fast load of the particles

  int rank;
  MPI_Comm_rank(comm, &rank);
  vpic_simulation_seed_entropy(rank);  //Generators desynchronized
  struct rng *rng = vpic_simulation_rng(0);

  assert(psc->nr_patches > 0);
  struct psc_patch *patch = &psc->patch[0];
  double xmin = patch->xb[0], xmax = patch->xb[0] + patch->dx[0] * patch->ldims[0];
  double ymin = patch->xb[1], ymax = patch->xb[1] + patch->dx[1] * patch->ldims[1];
  double zmin = patch->xb[2], zmax = patch->xb[2] + patch->dx[2] * patch->ldims[2];

  struct species *electron = vpic_simulation_find_species("electron");
  struct species *ion = vpic_simulation_find_species("ion");
  
  // Load Harris population

  mpi_printf(comm, "-> Main Harris Sheet\n");

  for (int64_t n = 0; n < Ne_sheet / n_global_patches; n++) {
    double x, y, z, ux, uy, uz, d0;

    particle_inject_t prt;

    do {
      z = L*atanh(vpic_simulation_uniform(rng, -1., 1.)*tanhf);
    } while (z <= zmin || z >= zmax);
    x = vpic_simulation_uniform(rng, xmin, xmax);
    y = vpic_simulation_uniform(rng, ymin, ymax);

    // inject_particles() will return an error for particles no on this
    // node and will not inject particle locally

    ux = vpic_simulation_normal(rng, 0, vthe);
    uy = vpic_simulation_normal(rng, 0, vthe);
    uz = vpic_simulation_normal(rng, 0, vthe);
    d0 = gdre*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udre;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

    prt.x[0] = x; prt.x[1] = y; prt.x[2] = z;
    prt.u[0] = ux; prt.u[1] = uy; prt.u[2] = uz;
    prt.w = weight_s;
    prt.kind = KIND_ELECTRON;
    inject_particle(electron, &prt, 0., 0);

    ux = vpic_simulation_normal(rng, 0, vthi);
    uy = vpic_simulation_normal(rng, 0, vthi);
    uz = vpic_simulation_normal(rng, 0, vthi);
    d0 = gdri*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udri;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

    prt.u[0] = ux; prt.u[1] = uy; prt.u[2] = uz;
    prt.kind = KIND_ION;
    inject_particle(ion, &prt, 0., 0);
  }

  mpi_printf(comm, "-> Background Population\n");

  for (int64_t n = 0; n < Ne_back / n_global_patches; n++) {
    particle_inject_t prt;
    double x = vpic_simulation_uniform(rng, xmin, xmax);
    double y = vpic_simulation_uniform(rng, ymin, ymax);
    double z = vpic_simulation_uniform(rng, zmin, zmax);

    prt.x[0] = x; prt.x[1] = y; prt.x[2] = z;
    prt.u[0] = vpic_simulation_normal(rng, 0, vtheb);
    prt.u[1] = vpic_simulation_normal(rng, 0, vtheb);
    prt.u[2] = vpic_simulation_normal(rng, 0, vtheb);
    prt.w = weight_b;
    prt.kind = KIND_ELECTRON;
    inject_particle(electron, &prt, 0., 0);
    
    prt.u[0] = vpic_simulation_normal(rng, 0, vthib);
    prt.u[1] = vpic_simulation_normal(rng, 0, vthib);
    prt.u[2] = vpic_simulation_normal(rng, 0, vthib);
    prt.kind = KIND_ION;
    inject_particle(ion, &prt, 0., 0);
  }

  mpi_printf(comm, "Finished loading particles\n");

  psc_mparticles_put_as(mprts_vpic, psc->particles, 0);
}


#define inject_particle(knd, x, y, z, ux, uy, uz, weight, a, b) do {	\
    double Vi = 1./(psc->patch[0].dx[0] * psc->patch[0].dx[1] * psc->patch[0].dx[2]); \
    particle_t prt;							\
    prt.xi = x - xmin;							\
    prt.yi = y - ymin;							\
    prt.zi = z - zmin;							\
    prt.pxi = ux;							\
    prt.pyi = uy;							\
    prt.pzi = uz;							\
    prt.qni_wni = psc->kinds[knd].q * weight * Vi;			\
    prt.kind = knd;							\
    mparticles_patch_push_back(mprts, p, prt);				\
  } while (0)


// ----------------------------------------------------------------------
// psc_harris_setup_particles

static void
psc_harris_setup_particles(struct psc *psc, int *nr_particles_by_patch, bool count_only)
{
  struct psc_harris *sub = psc_harris(psc);
  struct globals_physics *phys = &sub->phys;

  for (int p = 0; p < psc->nr_patches; p++) {
    nr_particles_by_patch[p] = phys->Ne / sub->n_global_patches;
  }

  if (count_only) {
    return;
  }

  double theta = sub->prm.theta;
  double L = phys->L;
  double tanhf = phys->tanhf;
  double vthe = phys->vthe, vthi = phys->vthi;
  double vtheb = phys->vtheb, vthib = phys->vthib;
  double gdre = phys->gdre, gdri = phys->gdri;
  double udre = phys->udre, udri = phys->udri;
  double weight_s = phys->weight_s;
  double weight_b = phys->weight_b;
  int electron = KIND_ELECTRON, ion = KIND_ION;
  double cs = cos(theta), sn = sin(theta);

  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, PARTICLE_TYPE, MP_DONT_COPY);

  for (int p = 0; p < psc->nr_patches; p++) {
    struct psc_patch *patch = &psc->patch[p];
    double xmin = patch->xb[0], xmax = patch->xb[0] + patch->ldims[0] * patch->dx[0];
    double ymin = patch->xb[1], ymax = patch->xb[1] + patch->ldims[1] * patch->dx[1];
    double zmin = patch->xb[2], zmax = patch->xb[2] + patch->ldims[2] * patch->dx[2];
  
    for (int n = 0; n < phys->Ne_sheet / sub->n_global_patches; n++) {
      double x, y, z, ux, uy, uz, d0 ;

      do {
	z = L*atanh(uniform(rng(0), -1, 1)*tanhf);
      } while( z<= zmin || z>=zmax );
      x = uniform( rng(0), xmin, xmax );
      y = uniform( rng(0), ymin, ymax );
      
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

    for (int n = 0; n < phys->Ne_back / sub->n_global_patches; n++) {

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

  }

  psc_mparticles_put_as(mprts, psc->particles, 0);
}

// ----------------------------------------------------------------------
// psc_harris_read

static void
psc_harris_read(struct psc *psc, struct mrc_io *io)
{
  psc_read_super(psc, io);
}

// ======================================================================
// psc_harris_ops

struct psc_ops psc_harris_ops = {
  .name             = "harris",
  .size             = sizeof(struct psc_harris),
  .param_descr      = psc_harris_descr,
  .create           = psc_harris_create,
  .setup            = psc_harris_setup,
  .read             = psc_harris_read,
  .init_field       = psc_harris_init_field,
  .setup_particles  = psc_harris_setup_particles,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_harris_ops);
}
