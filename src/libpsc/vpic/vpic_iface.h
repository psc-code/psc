
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

#include "../bk_mparticles_iface.h" // FIXME, path

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

struct vpic_mparticles;

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

struct vpic_mfields;

struct vpic_mfields *vpic_mfields_create();
void vpic_mfields_ctor_from_simulation_fields(struct vpic_mfields *vmflds);
void vpic_mfields_ctor_from_simulation_hydro(struct vpic_mfields *vmflds);
float *vpic_mfields_get_data(struct vpic_mfields *mflds, int *ib, int *im);
double vpic_mfields_synchronize_tang_e_norm_b(struct vpic_mfields *mflds);
void vpic_mfields_compute_div_b_err(struct vpic_mfields *vmflds);
double vpic_mfields_compute_rms_div_b_err(struct vpic_mfields *vmflds);
void vpic_mfields_clean_div_b(struct vpic_mfields *vmflds);
void vpic_mfields_compute_div_e_err(struct vpic_mfields *vmflds);
double vpic_mfields_compute_rms_div_e_err(struct vpic_mfields *vmflds);
void vpic_mfields_clean_div_e(struct vpic_mfields *vmflds);
void vpic_mfields_clear_rhof(struct vpic_mfields *vmflds);
void vpic_mfields_accumulate_rho_p(struct vpic_mfields *vmflds, struct vpic_mparticles *mprts);
void vpic_mfields_synchronize_rho(struct vpic_mfields *vmflds);
void vpic_mfields_compute_rhob(struct vpic_mfields *vmflds);
void vpic_mfields_compute_curl_b(struct vpic_mfields *vmflds);

// ----------------------------------------------------------------------
// vpic_mparticles

struct vpic_mparticles_prt {
  float dx[3];
  int i;
  float ux[3];
  float w;
  int kind;
};

struct vpic_mparticles *vpic_mparticles_create();
void vpic_mparticles_ctor_from_simulation(struct vpic_mparticles *vmprts);
int vpic_mparticles_get_nr_particles(struct vpic_mparticles *vmprts);
void vpic_mparticles_reserve_all(struct vpic_mparticles *vmprts, int n_patches,
				 int *n_prts_by_patch);
void vpic_mparticles_resize_all(struct vpic_mparticles *vmprts, int n_patches,
				 int *n_prts_by_patch);
void vpic_mparticles_get_size_all(struct vpic_mparticles *vmprts, int n_patches,
				  int *n_prts_by_patch);
void vpic_mparticles_get_particles(struct vpic_mparticles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*put_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx);
void vpic_mparticles_set_particles(struct vpic_mparticles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*get_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx);
void vpic_mparticles_push_back(struct vpic_mparticles *vmprts, const struct vpic_mparticles_prt *prt);
void vpic_mparticles_get_grid_nx_dx(struct vpic_mparticles *vmprts, int *nx, float *dx);

void vpic_mparticles_copy_to_single_by_kind(struct vpic_mparticles *vmprts, bk_mparticles *bkmprts);
void vpic_mparticles_copy_from_single_by_kind(struct vpic_mparticles *vmprts, bk_mparticles *bkmprts);

// ----------------------------------------------------------------------
// vpic_push_particles

struct vpic_push_particles;

struct vpic_push_particles *vpic_push_particles_create();
void vpic_push_particles_ctor_from_simulation(struct vpic_push_particles *vpushp);
void vpic_push_particles_push_mprts(struct vpic_push_particles *vpushp,
				    struct vpic_mparticles *vmprts,
				    struct vpic_mfields *vmflds);
void vpic_push_particles_stagger_mprts(struct vpic_push_particles *vpushp,
				       struct vpic_mparticles *vmprts,
				       struct vpic_mfields *vmflds);
void vpic_push_particles_prep(struct vpic_push_particles *vpushp,
			      struct vpic_mparticles *mprts,
			      struct vpic_mfields *vmflds);

// ----------------------------------------------------------------------
// vpic_push_fields

void vpic_push_fields_advance_b(struct vpic_mfields *vmflds, double frac);
void vpic_push_fields_advance_e(struct vpic_mfields *vmflds, double frac);

// ----------------------------------------------------------------------
// other (may want an object eventually)

void vpic_sort_run(struct vpic_mparticles *vmprts, int step);
void vpic_collision_run();
void vpic_emitter();
void vpic_current_injection();
void vpic_field_injection();
void vpic_moments_run(struct vpic_mfields *mflds, struct vpic_mparticles *vmprts, int kind);

// ----------------------------------------------------------------------
// vpic_params

struct vpic_params {
  double cfl_req;

  int np[3];    // domain topology
  int gdims[3]; // global dims

  int status_interval;

  int quota_check_interval;  // How frequently to check if quota exceeded
  double quota;              // Run quota in hours

  int restart_interval;
};

struct vpic_kind_info {
  double q;
  double m;
  char *name;
};

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
void vpic_simulation_init_split(struct vpic_params *vpic_prm,
				struct psc_harris *harris,
				struct vpic_simulation_info *info);
void vpic_simulation_init(struct vpic_simulation_info *info);


void vpic_print_status();
void vpic_diagnostics();
void vpic_diagnostics_split(struct vpic_params *vpic_prm, struct psc_harris *harris);
void vpic_inc_step(int step);


#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#ifndef mprintf

#include <mpi.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

#endif

#endif
