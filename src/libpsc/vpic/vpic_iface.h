
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

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
void vpic_mparticles_get_size_all(struct vpic_mparticles *vmprts, int n_patches,
				  int *n_prts_by_patch);
void vpic_mparticles_get_particles(struct vpic_mparticles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*put_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx);
void vpic_mparticles_get_grid_nx_dx(struct vpic_mparticles *vmprts, int *nx, float *dx);

// ----------------------------------------------------------------------
// vpic_push_particles

struct vpic_push_particles;

struct vpic_push_particles *vpic_push_particles_create();
void vpic_push_particles_ctor_from_simulation(struct vpic_push_particles *vpushp);
void vpic_push_particles_push_mprts(struct vpic_push_particles *vpushp,
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

struct vpic_kind_info {
  double q;
  double m;
  char *name;
};

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
void vpic_simulation_init(struct vpic_simulation_info *info);
void vpic_simulation_init2(struct vpic_push_particles *vpushp, struct vpic_mfields *vmflds,
			   struct vpic_mparticles *vmprts);


void vpic_print_status();
void vpic_diagnostics();
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
