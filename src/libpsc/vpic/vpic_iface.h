
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

#include <stdbool.h>

#include "../bk_mparticles_iface.h" // FIXME, path

#ifdef __cplusplus

#include "simulation.h"

#include "VpicFieldArrayBase.h"
#include "PscFieldArrayBase.h"
#include "VpicFieldArrayLocalOps.h"
#include "PscFieldArrayLocalOps.h"
#include "VpicFieldArray.h"
#include "PscFieldArray.h"

#include "VpicParticlesBase.h"
#include "VpicParticlesOps.h"
#include "PscParticlesOps.h"

#include "VpicInterpolator.h"
#include "VpicInterpolatorOps.h"
#include "PscInterpolatorOps.h"

#include "VpicAccumulator.h"
#include "VpicAccumulatorOps.h"
#include "PscAccumulatorOps.h"

#include "VpicHydroArrayBase.h"
#include "VpicHydroArray.h"
#include "PscHydroArray.h"

#include "VpicDiag.h"
#include "NoneDiag.h"

#include "VpicSimulationBase.h"
#include "PscSimulationBase.h"

#endif

BEGIN_C_DECLS

// ----------------------------------------------------------------------
// vpic_mfields

enum {
  VPIC_MFIELDS_EX = 0,
  VPIC_MFIELDS_EY = 1,
  VPIC_MFIELDS_EZ = 2,
  VPIC_MFIELDS_DIV_E_ERR = 3,
  VPIC_MFIELDS_BX = 4,
  VPIC_MFIELDS_BY = 5,
  VPIC_MFIELDS_BZ = 6,
  VPIC_MFIELDS_DIV_B_ERR = 7,
  VPIC_MFIELDS_TCAX = 8,
  VPIC_MFIELDS_TCAY = 9,
  VPIC_MFIELDS_TCAZ = 10,
  VPIC_MFIELDS_RHOB = 11,
  VPIC_MFIELDS_JFX = 12,
  VPIC_MFIELDS_JFY = 13,
  VPIC_MFIELDS_JFZ = 14,
  VPIC_MFIELDS_RHOF = 15,
  VPIC_MFIELDS_N_COMP = 20,
};

enum {
  VPIC_HYDRO_N_COMP = 16,
};

#ifdef __cplusplus

#if 1
typedef PscFieldArrayBase FieldArrayBase;
typedef PscFieldArrayLocalOps<FieldArrayBase> FieldArrayLocalOps;
typedef PscFieldArray<FieldArrayBase, FieldArrayLocalOps> FieldArray;
#else
typedef VpicFieldArrayBase FieldArrayBase;
typedef VpicFieldArray<FieldArrayBase> FieldArray;
#endif

typedef VpicInterpolator Interpolator;
typedef PscInterpolatorOps<Interpolator, FieldArrayBase> InterpolatorOps;

typedef VpicAccumulator Accumulator;
typedef PscAccumulatorOps<Accumulator, FieldArrayBase> AccumulatorOps;

typedef VpicHydroArrayBase HydroArrayBase;
typedef PscHydroArray<HydroArrayBase> HydroArray;

typedef VpicParticlesBase ParticlesBase;  
typedef PscParticlesOps<ParticlesBase, FieldArrayBase, Interpolator, Accumulator> ParticlesOps;
typedef ParticlesBase Particles;

typedef VpicDiagOps<FieldArray, Particles, Interpolator, HydroArray> DiagOps;
//typedef VpicSimulationBase SimulationBase;
typedef PscSimulationBase<FieldArray, Particles, Interpolator, Accumulator, HydroArray> SimulationBase;


typedef VpicSimulation<FieldArray, ParticlesOps, InterpolatorOps, AccumulatorOps, HydroArray,
		       SimulationBase, DiagOps> Simulation;

#else

typedef struct FieldArray_ FieldArray;
typedef struct Particles_ Particles;
typedef struct Simulation_ Simulation;
typedef struct HydroArray_ HydroArray;

#endif

void vpic_mfields_accumulate_rho_p(FieldArray *vmflds, Particles *mprts);

// ----------------------------------------------------------------------
// vpic_mparticles

struct vpic_mparticles_prt {
  float dx[3];
  int i;
  float ux[3];
  float w;
  int kind;
};

struct psc_particle_inject;

Particles *vpic_mparticles_new_from_simulation(Simulation *sim);
int vpic_mparticles_get_nr_particles(Particles *vmprts);
void vpic_mparticles_reserve_all(Particles *vmprts, int n_patches,
				 int *n_prts_by_patch);
void vpic_mparticles_resize_all(Particles *vmprts, int n_patches,
				 int *n_prts_by_patch);
void vpic_mparticles_get_size_all(Particles *vmprts, int n_patches,
				  int *n_prts_by_patch);
void vpic_mparticles_get_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*put_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx);
void vpic_mparticles_set_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*get_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx);
void vpic_mparticles_push_back(Particles *vmprts, const struct vpic_mparticles_prt *prt);
void vpic_mparticles_get_grid_nx_dx(Particles *vmprts, int *nx, float *dx);
void vpic_mparticles_sort(Particles *vmprts, int step);

void vpic_mparticles_copy_to_single_by_kind(Particles *vmprts, bk_mparticles *bkmprts);
void vpic_mparticles_copy_from_single_by_kind(Particles *vmprts, bk_mparticles *bkmprts);

// ----------------------------------------------------------------------
// vpic_push_particles

struct vpic_push_particles;

struct vpic_push_particles *vpic_push_particles_new_from_Simulation(Simulation *sim);
void vpic_push_particles_push_mprts(struct vpic_push_particles *vpushp,
				    Particles *vmprts,
				    FieldArray *vmflds);
void vpic_push_particles_stagger_mprts(struct vpic_push_particles *vpushp,
				       Particles *vmprts,
				       FieldArray *vmflds);
void vpic_push_particles_prep(struct vpic_push_particles *vpushp,
			      Particles *mprts,
			      FieldArray *vmflds);

// ----------------------------------------------------------------------
// Simulation

struct vpic_simulation_info;
struct field_array;

Simulation *Simulation_create();
void Simulation_delete(Simulation *sim);

void Simulation_set_params(Simulation *vpic,
			   int num_step, int status_interval, int sync_shared_interval,
			   int clean_div_e_interval, int clean_div_b_interval);

void Simulation_setup_grid(Simulation *sim, double dx[3], double dt,
			   double cvac, double eps0);
void Simulation_define_periodic_grid(Simulation *sim, double xl[3],
				     double xh[3], int gdims[3], int np[3]);
void Simulation_set_domain_field_bc(Simulation *sim, int boundary, int fbc);
void Simulation_set_domain_particle_bc(Simulation *sim, int boundary, int pbc);

struct material *Simulation_define_material(Simulation *sim, const char *name,
					    double eps, double mu,
					    double sigma, double zeta);
void Simulation_define_field_array(Simulation *sim, double damp);
struct species * Simulation_define_species(Simulation *sim, const char *name, double q, double m,
					   double max_local_np, double max_local_nm,
					   double sort_interval, double sort_out_of_place);
void Simulation_inject_particle(Simulation *sim, Particles *vmprts, int p,
				const struct psc_particle_inject *prt);
void Simulation_initialize(Simulation *sim, Particles *vmprts, FieldArray *vmflds);
void Simulation_collision_run(Simulation *sim);
void Simulation_field_injection(Simulation *sim);
void Simulation_moments_run(Simulation *sim, HydroArray *mflds, Particles *vmprts, int kind);
void Simulation_advance_b(Simulation *sim, FieldArray *vmflds, double frac);
void Simulation_advance_e(Simulation *sim, FieldArray *vmflds, double frac);


void Simulation_diagnostics_init(Simulation *sim, int interval);
void Simulation_diagnostics_setup(Simulation *sim);

void Simulation_get_info(Simulation *sim, struct vpic_simulation_info *info);
void Simulation_print_status(Simulation *sim);
void Simulation_diagnostics(Simulation *sim);
void Simulation_inc_step(Simulation *sim, int step);
void Simulation_user_initialization(Simulation *sim);

// Harris specific
struct psc_harris;
struct vpic_harris_params;
struct globals_physics;

void Simulation_diagnostics_run(Simulation *sim);
void Simulation_set_region_resistive_harris(Simulation *sim,
					    struct vpic_harris_params *prm,
					    struct globals_physics *phys,
					    double dx[3],
					    double thickness,
					    struct material *resistive);



void Simulation_rngPool_seed(Simulation *sim, int base);
struct Rng *Simulation_rngPool_get(Simulation *sim, int n);

double Rng_uniform(struct Rng *rng, double lo, double hi);
double Rng_normal(struct Rng *rng, double mu, double sigma);

// ----------------------------------------------------------------------
// vpic_harris_params

struct vpic_harris_params {
  // general
  double wpedt_max;

  double wpe_wce;                 // electron plasma freq / electron cyclotron freq
  double mi_me;                   // Ion mass / electron mass
  
  double Lx_di, Ly_di, Lz_di;     // Size of box in d_i

  int ion_sort_interval;
  int electron_sort_interval;
  double taui;                    // simulation wci's to run
  double t_intervali;             // output interval in terms of 1/wci
  double output_field_interval;   // field output interval in terms of 1/wci
  double output_particle_interval;// particle output interval in terms of 1/wci

  // Harris
  double L_di;                    // Sheet thickness / ion inertial length
  double Ti_Te;                   // Ion temperature / electron temperature
  double nb_n0;                   // background plasma density
  double Tbe_Te;                  // Ratio of background T_e to Harris T_e
  double Tbi_Ti;                  // Ratio of background T_i to Harris T_i
  double bg;                      // Guide field
  double theta;

  double Lpert_Lx;                // wavelength of perturbation in terms of Lx
  double dbz_b0;                  // perturbation in Bz relative to B0
  double nppc;                    // Average number of macro particle per cell per species
  bool open_bc_x;                 // Flag to signal we want to do open boundary condition in x
  bool driven_bc_z;               // Flag to signal we want to do driven boundary condition in z

  double overalloc;               // Overallocation factor (> 1) for particle arrays
};

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

// ----------------------------------------------------------------------
// psc_harris

struct psc_harris {
  // must be first because of our hacky way of passing arguments to
  // vpic initialization deck
  // params
  struct vpic_harris_params prm;
  
  // state
  struct globals_physics phys;
  int n_global_patches;

  Simulation *sim;
};

#define psc_harris(psc) mrc_to_subobj(psc, struct psc_harris)

// ----------------------------------------------------------------------
// vpic_kind_info

struct vpic_kind_info {
  double q;
  double m;
  char *name;
};

// ----------------------------------------------------------------------
// vpic_simulation_info
//
// returned from vpic_simulation_init

struct vpic_simulation_info {
  int num_step;
  double dt;
  int nx[3];
  double dx[3];
  double x0[3];
  double x1[3];

  int n_kinds;
  struct vpic_kind_info *kinds;

  int clean_div_e_interval;
  int clean_div_b_interval;
  int sync_shared_interval;
  int num_div_e_round;
  int num_div_b_round;

  int status_interval;
};

void vpic_base_init(int *pargc, char ***pargv);

// FIXME, replicated
#define BOUNDARY(i,j,k) (13+(i)+3*(j)+9*(k)) /* FORTRAN -1:1,-1:1,-1:1 */

END_C_DECLS

#ifndef mprintf

#include <mpi.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

#endif

#endif
