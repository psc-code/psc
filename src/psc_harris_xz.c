
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_collision.h>
#include <psc_balance.h>
#include <psc_marder.h>
#include <psc_sort.h>
#include <psc_collision.h>
#include <psc_bnd_particles.h>
#include <psc_checks.h>
#include <psc_event_generator.h>

#include <psc_particles_as_single.h>

#include <mrc_params.h>

#include <math.h>

#include "libpsc/vpic/vpic_iface.h"

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

struct psc_harris {
  // params
  
  // state
  int n_global_patches;

  double Npe_sheet, Npe_back, Npe;
  double Ne_sheet, Ne_back, Ne;
  double weight_s, weight_b;
  double vthe, vthi;
  double vtheb, vthib;
  double L;
  double gdre;
  double gdri;
  double udre;
  double udri;
  double tanhf;
  double sn, cs;
  double b0, bg;
  double dbx, dbz;
  double Lx, Ly, Lz;
  double Lpert;
};

#define psc_harris(psc) mrc_to_subobj(psc, struct psc_harris)

#define VAR(x) (void *)offsetof(struct psc_harris, x)
static struct param psc_harris_descr[] = {
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_harris_create

static void
psc_harris_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
  psc->prm.nicell = 1;
  psc->prm.cfl = 0.99;

  // will be set to actual values in psc_harris_setup()
  psc->domain.length[0] = 1.; // no x dependence 
  psc->domain.length[1] = 1.;
  psc->domain.length[2] = 1.;

  psc->domain.gdims[0] = 64;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 32;

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

  psc_sort_set_type(psc->sort, "vpic");
  // FIXME: the "vpic" sort actually keeps track of per-species sorting intervals
  // internally
  psc_sort_set_param_int(psc->sort, "every", 1);

  psc_collision_set_type(psc->collision, "vpic");

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc->push_particles, "vpic");
  psc_bnd_particles_set_type(psc->bnd_particles, "vpic");

  psc_push_fields_set_type(psc->push_fields, "vpic");
  psc_marder_set_type(psc->marder, "vpic");
  // FIXME, marder "vpic" manages its own cleaning intervals
  psc_marder_set_param_int(psc->marder, "every_step", 1);
}

static inline double trunc_granular( double a, double b )
{
  return b * (int)(a/b);
}

// ----------------------------------------------------------------------
// psc_harris_setup

static void
psc_harris_setup(struct psc *psc)
{
  struct psc_harris *sub = psc_harris(psc);

  // use natural PIC units
  double ec   = 1;         // Charge normalization
  double me   = 1;         // Mass normalization
  double c    = 1;         // Speed of light
  /* double de   = 1;         // Length normalization (electron inertial length) */
  double eps0 = 1;         // Permittivity of space

#if 0
  double cfl_req   = 0.99;  // How close to Courant should we try to run
#endif
  double wpedt_max = 0.36;  // Max dt allowed if Courant not too restrictive
#if 0
  double damp      = 0.0;   // Level of radiation damping
  int rng_seed     = 1;     // Random number seed increment
#endif

  // Physics parameters
  double mi_me   = 25.0;  // Ion mass / electron mass
  double L_di    = 0.5;    // Sheet thickness / ion inertial length
  double Ti_Te   = 5.0;    // Ion temperature / electron temperature
  /* double Z   = 1.0;      // Ion charge */
  double nb_n0   = 0.228;   // background plasma density
  double Tbe_Te  = 0.7598;  // Ratio of background T_e to Harris T_e
  double Tbi_Ti  = 0.3039;  // Ratio of background T_i to Harris T_i
  double wpe_wce = 2.0;    // electron plasma freq / electron cyclotron freq
  double bg = 0.0;
  double theta   = 0;      // B0 = Bx

  double cs   = cos(theta);
  double sn   = sin(theta);

  //derived quantities
  double mi = me*mi_me;       // Ion mass
  double Te = me*c*c/(2*eps0*wpe_wce*wpe_wce*(1+Ti_Te)); // Electron temperature
  double Ti = Te*Ti_Te;       // Ion temperature
  double vthe = sqrt(Te/me);                        // Electron thermal velocity
  double vthi = sqrt(Ti/mi);  // Ion thermal velocity
  double vtheb = sqrt(Tbe_Te*Te/me);  // normalized background e thermal vel. */
  double vthib = sqrt(Tbi_Ti*Ti/mi);  // normalized background ion thermal vel. */
  double wci  = 1.0/(mi_me*wpe_wce);  // Ion cyclotron frequency
  double wce  = wci*mi_me;            // Electron cyclotron freqeuncy
  double wpe  = wce*wpe_wce;          // electron plasma frequency
  double wpi  = wpe/sqrt(mi_me);      // ion plasma frequency
  double di   = c/wpi;                // ion inertial length
  double L    = L_di*di;              // Harris sheet thickness
  /* double rhoi_L = sqrt(Ti_Te/(1.0+Ti_Te))/L_di; */
  /* double v_A= (wci/wpi)/sqrt(nb_n0); // based on nb */

  // Numerical parameters
  double nppc  = 100; // Average number of macro particle per cell per species

  double Lx    = 25.6*di;      // size of box in x dimension
  double Ly    = 1*di; // size of box in y dimension
  double Lz    = 12.8*di;      // size of box in z dimension

  sub->L = L;
  sub->Lx = Lx;
  sub->Ly = Ly;
  sub->Lz = Lz;

  int *gdims = psc->domain.gdims, *np = psc->domain.np;
  sub->n_global_patches = np[0] * np[1] * np[2];
  
  double b0 = me*c*wce/ec; // Asymptotic magnetic field strength
  double n0 = me*eps0*wpe*wpe/(ec*ec);  // Peak electron (ion) density
  double vdri = 2*c*Ti/(ec*b0*L);   // Ion drift velocity
  double vdre = -vdri/(Ti_Te);      // electron drift velocity
  
  double Npe_sheet = 2*n0*Lx*Ly*L*tanh(0.5*Lz/L); // N physical e's in sheet
  double Npe_back  = nb_n0*n0*Ly*Lz*Lx;           // N physical e's in backgrnd
  double Npe       = Npe_sheet + Npe_back;
  double Ne        = nppc*gdims[0]*gdims[1]*gdims[2];  // total macro electrons in box
  double Ne_sheet  = Ne*Npe_sheet/Npe;
  double Ne_back   = Ne*Npe_back/Npe;
  Ne_sheet = trunc_granular(Ne_sheet, sub->n_global_patches); // Make it divisible by nproc
  Ne_back  = trunc_granular(Ne_back, sub->n_global_patches);  // Make it divisible by nproc
  Ne = Ne_sheet + Ne_back;
  //double qe_s = -ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  //double qi_s =  ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  double weight_s = ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  //double qe_b = -ec*Npe_back/Ne_back;  // Charge per macro electron
  //double qi_b =  ec*Npe_back/Ne_back;  // Charge per macro electron
  double weight_b =  ec*Npe_back/Ne_back;  // Charge per macro electron

  sub->b0 = b0;
  sub->bg = bg;
  sub->Npe_sheet = Npe_sheet;
  sub->Npe_back = Npe_back;
  sub->Npe = Npe;
  sub->Ne_sheet = Ne_sheet;
  sub->Ne_back = Ne_back;
  sub->Ne = Ne;
  sub->weight_s = weight_s;
  sub->weight_b = weight_b;
  sub->vthe = vthe;
  sub->vthi = vthi;
  sub->vtheb = vtheb;
  sub->vthib = vthib;
  sub->cs = cs;
  sub->sn = sn;

  psc->domain.length[0] = Lx;
  psc->domain.length[1] = Ly;
  psc->domain.length[2] = Lz;

  psc->domain.corner[0] = 0.;
  psc->domain.corner[1] = -.5 * Ly;
  psc->domain.corner[2] = -.5 * Lz;
  
  double gdri = 1/sqrt(1-vdri*vdri/(c*c));  // gamma of ion drift frame
  double gdre = 1/sqrt(1-vdre*vdre/(c*c)); // gamma of electron drift frame
  double udri = vdri*gdri;                 // 4-velocity of ion drift frame
  double udre = vdre*gdre;                 // 4-velocity of electron drift frame

  sub->gdre = gdre;
  sub->gdri = gdri;
  sub->udre = udre;
  sub->udri = udri;

  sub->L = L;
  sub->tanhf = tanh(0.5*Lz/L);
  sub->Lpert = Lx; // wavelength of perturbation
  sub->dbz =  0.03*b0; // Perturbation in Bz relative to Bo (Only change here)
  sub->dbx = -sub->dbz * sub->Lpert/(2.0*Lz); // Set Bx perturbation so that div(B) = 0

  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "b0 = %g\n", b0);
  mpi_printf(comm, "n0 = %g\n", n0);
  mpi_printf(comm, "Ti = %g Te = %g\n", Ti, Te);
  mpi_printf(comm, "vthi = %g vthe = %g\n", vthi, vthe);
  mpi_printf(comm, "vdri = %g vdre = %g\n", vdri, vdre);
  mpi_printf(comm, "Npe = %g, Npe_sheet = %g Npe_back = %g\n",
	     sub->Npe, sub->Npe_sheet, sub->Npe_back);
  mpi_printf(comm, "Ne = %g, Ne_sheet = %g Ne_back = %g\n",
	     sub->Ne, sub->Ne_sheet, sub->Ne_back);
  mpi_printf(comm, "weight_s = %g, weight_b = %g\n",
	     sub->weight_s, sub->weight_b);
  
  // initializes fields, particles, etc.
  psc_setup_super(psc);

  if (wpe*psc->dt > wpedt_max) {
    psc->dt = wpedt_max/wpe;  // override timestep if plasma frequency limited
    mprintf("::: dt reduced to %g\n", psc->dt);
  }
  
  //  psc->dt = 0.0714471;

#if 0
  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "dt = %g, dy = %g dz = %g\n", psc->dt,
	     psc->patch[0].dx[1], psc->patch[0].dx[2]);
  mpi_printf(comm, "d_e = %g, d_i = %g\n", 1., d_i);
  mpi_printf(comm, "v_A = %g\n", harris->B0 / sqrt(harris->mi_over_me));
  double om_ci = harris->B0 / harris->mi_over_me;
  double om_ce = harris->B0;
  mpi_printf(comm, "om_ci = %g om_ce = %g\n", om_ci, om_ce);
  mpi_printf(comm, "om_pi = %g om_pe = %g\n", 1. / sqrt(harris->mi_over_me), 1.);
  mpi_printf(comm, "Ti = %g Te = %g\n", harris->Ti, harris->Te);
  double vthi = sqrt(2 * harris->Ti / harris->mi_over_me);
  double vthe = sqrt(2 * harris->Te);
  mpi_printf(comm, "vthi = %g vthe = %g\n", vthi, vthe);
  mpi_printf(comm, "rhoi = %g\n", vthi / om_ci);
  mpi_printf(comm, "L / rhoi = %g\n", harris->LLL / (vthi / om_ci));
  mpi_printf(comm, "pert A0 / (B0 d_i) = %g\n", harris->AA / (harris->B0 * d_i));
  mpi_printf(comm, "lambda_De = %g\n", sqrt(harris->Te));
#endif
}

// ----------------------------------------------------------------------
// psc_harris_init_field

static double
psc_harris_init_field(struct psc *psc, double crd[3], int m)
{
  struct psc_harris *sub = psc_harris(psc);
  double cs = sub->cs, sn = sub->sn;
  double b0 = sub->b0, bg = sub->bg, dbx = sub->dbx, dbz = sub->dbz;
  double L = sub->L, Lx = sub->Lx, Lz = sub->Lz, Lpert = sub->Lpert;
  double x = crd[0], z = crd[2];

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

static void
psc_harris_setup_particles(struct psc *psc, int *nr_particles_by_patch, bool count_only)
{
  struct psc_harris *sub = psc_harris(psc);
  for (int p = 0; p < psc->nr_patches; p++) {
    nr_particles_by_patch[p] = sub->Ne / sub->n_global_patches;
  }

  if (count_only) {
    return;
  }

  double L = sub->L;
  double tanhf = sub->tanhf;
  double cs = sub->cs, sn = sub->sn;
  double vthe = sub->vthe, vthi = sub->vthi;
  double vtheb = sub->vtheb, vthib = sub->vthib;
  double gdre = sub->gdre, gdri = sub->gdri;
  double udre = sub->udre, udri = sub->udri;
  double weight_s = sub->weight_s;
  double weight_b = sub->weight_b;
  int electron = KIND_ELECTRON, ion = KIND_ION;

  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, PARTICLE_TYPE, MP_DONT_COPY);

  for (int p = 0; p < psc->nr_patches; p++) {
    struct psc_patch *patch = &psc->patch[p];
    double xmin = patch->xb[0], xmax = patch->xb[0] + patch->ldims[0] * patch->dx[0];
    double ymin = patch->xb[1], ymax = patch->xb[1] + patch->ldims[1] * patch->dx[1];
    double zmin = patch->xb[2], zmax = patch->xb[2] + patch->ldims[2] * patch->dx[2];
  
    for (int n = 0; n < sub->Ne_sheet / sub->n_global_patches; n++) {
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

    for (int n = 0; n < sub->Ne_back / sub->n_global_patches; n++) {

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
// PSC/VPIC coupling that should eventually move out of this file

// ----------------------------------------------------------------------
// psc_harris_output

static void
psc_harris_output(struct psc *psc)
{
  // currently, vpic does its first output as part of its own initialize,
  // so let's not do it again here
  if (psc->timestep > 0) {
    vpic_inc_step(psc->timestep);
    vpic_print_status();
    vpic_diagnostics();
  }
}

// ----------------------------------------------------------------------
// psc_harris_step

static void
psc_harris_step(struct psc *psc)
{
  psc_sort_run(psc->sort, psc->particles);
  psc_collision_run(psc->collision, psc->particles);

  psc_checks_continuity_before_particle_push(psc->checks, psc);

  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);

  psc_push_fields_push_H(psc->push_fields, psc->flds, .5);

  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);

  psc_event_generator_run(psc->event_generator, psc->particles, psc->flds);
  
  psc_push_fields_push_E(psc->push_fields, psc->flds, 1.); // FIXME, should be _b2
  psc_push_fields_push_H(psc->push_fields, psc->flds, .5); // FIXME, should be _a

  psc_checks_continuity_after_particle_push(psc->checks, psc);

  psc_marder_run(psc->marder, psc->flds, psc->particles);

  psc_checks_gauss(psc->checks, psc);

  // Fields are updated ... load the interpolator for next time step and
  // particle diagnostics in user_diagnostics if there are any particle
  // species to worry about
  psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
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
  .step             = psc_harris_step,
  .output           = psc_harris_output,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  //  return psc_main(&argc, &argv, &psc_harris_ops);

  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);

  mrc_class_register_subclass(&mrc_class_psc, &psc_harris_ops);

  int from_checkpoint = -1;
  mrc_params_get_option_int("from_checkpoint", &from_checkpoint);

  struct psc *psc;

  if (from_checkpoint < 0) {
    // regular start-up (not from checkpoint)

    // psc_create() will create the psc object, create the sub-objects
    // (particle, field pusher and many others) and set the parameter defaults.
    // It will also set the psc subtype defaults and call psc_subtype_create(),
    // which will change some of the general defaults to match this case.
    psc = psc_create(MPI_COMM_WORLD);
    
    // psc_set_from_options() will override general and bubble psc parameters
    // if given on the command line. It will also call
    // psc_bubble_set_from_options()
    psc_set_from_options(psc);
    
    // psc_setup() will set up the various sub-objects (particle pusher, ...)
    // and set up the initial domain partition, the particles and the fields.
    // The standard implementation, used here, will set particles using
    // psc_bubble_init_npt and the fields using setup_field()
    psc_setup(psc);

    struct vpic_simulation_info info;
    vpic_base_init(&info);

    psc->prm.nmax = info.num_step;
    psc_marder_set_param_int(psc->marder, "clean_div_e_interval",
			     info.clean_div_e_interval);
    psc_marder_set_param_int(psc->marder, "clean_div_b_interval",
			     info.clean_div_b_interval);
    psc_marder_set_param_int(psc->marder, "sync_shared_interval",
			     info.sync_shared_interval);
    
  } else {
    // get psc object from checkpoint file

    psc = psc_read_checkpoint(MPI_COMM_WORLD, from_checkpoint);

    int checkpoint_nmax = -1;
    mrc_params_get_option_int("checkpoint_nmax", &checkpoint_nmax);
    if (checkpoint_nmax >= 0) {
      psc->prm.nmax = checkpoint_nmax;
    }

    int checkpoint_wallclock_limit = -1;
    mrc_params_get_option_int("checkpoint_wallclock_limit", &checkpoint_wallclock_limit);
    if (checkpoint_wallclock_limit >= 0) {
      psc->prm.wallclock_limit = checkpoint_wallclock_limit;
    }
  }

  // psc_view() will just print a whole lot of info about the psc object and
  // sub-objects, in particular all the parameters.
  psc_view(psc);
  
  // psc_integrate() uses the standard implementation, which does the regular
  // classic PIC time integration loop
  psc_integrate(psc);
  
  // psc_destroy() just cleans everything up when we're done.
  psc_destroy(psc);
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
  
}
