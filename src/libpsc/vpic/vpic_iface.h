
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

#include "vpic_config.h"

#include <stdbool.h>

#include <mpi.h>
#include <mrc_common.h>

#include "../bk_mparticles_iface.h" // FIXME, path

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

template <typename F>
void vpic_mparticles_get_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   F setter);
template <typename F>
void vpic_mparticles_set_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   F getter);

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
				     double xh[3], const int gdims[3], const int np[3]);
void Simulation_set_domain_field_bc(Simulation *sim, int boundary, int fbc);
void Simulation_set_domain_particle_bc(Simulation *sim, int boundary, int pbc);

struct material *Simulation_define_material(Simulation *sim, const char *name,
					    double eps, double mu,
					    double sigma, double zeta);
void Simulation_define_field_array(Simulation *sim, double damp);
struct species * Simulation_define_species(Simulation *sim, const char *name, double q, double m,
					   double max_local_np, double max_local_nm,
					   double sort_interval, double sort_out_of_place);

Particles* Simulation_get_particles(Simulation *sim);
int Simulation_mprts_get_nr_particles(Simulation* sim, Particles* vmprts);
void Simulation_mprts_reserve_all(Simulation* sim, Particles* vmprts, int n_patches,
				  const uint* n_prts_by_patch);
void Simulation_mprts_resize_all(Simulation* sim, Particles* vmprts, int n_patches,
				 const uint* n_prts_by_patch);
void Simulation_mprts_push_back(Simulation* sim, Particles* vmprts, const struct vpic_mparticles_prt *prt);

void Simulation_inject_particle(Simulation *sim, Particles *vmprts, int p,
				const struct psc_particle_inject *prt);

void Simulation_initialize(Simulation *sim, Particles *vmprts, FieldArray *vmflds);
void Simulation_moments_run(Simulation *sim, HydroArray *mflds, Particles *vmprts, int kind);
void Simulation_accumulate_rho_p(Simulation *sim, Particles *mprts, FieldArray *vmflds);

void Simulation_diagnostics_init(Simulation *sim, int interval);
void Simulation_diagnostics_setup(Simulation *sim);
void Simulation_diagnostics_run(Simulation *sim);

// Harris specific
void Simulation_set_region_resistive_harris(Simulation *sim,
					    struct vpic_harris_params *prm,
					    struct globals_physics *phys,
					    double dx[3],
					    double thickness,
					    struct material *resistive);



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

#ifndef mprintf

#include <mpi.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

#endif

// ======================================================================
// vpic_mparticles implementation

// ----------------------------------------------------------------------
// vpic_mparticles_get_particles

template<typename F>
void vpic_mparticles_get_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   F getter)
{
  unsigned int v_off = 0;
  for (auto sp = vmprts->cbegin(); sp != vmprts->cend(); ++sp) {
    unsigned int v_n_prts = sp->np;

    unsigned int nb = std::max(v_off, off), ne = std::min(v_off + v_n_prts, off + n_prts);
    for (unsigned int n = nb; n < ne; n++) {
      Particles::Particle* p = &sp->p[n - v_off];
#if 0
      int i = p->i;
      int im[3] = { sp->g->nx + 2, sp->g->ny + 2, sp->g->nz + 2 };
      int i3[3];
      i3[2] = i / (im[0] * im[1]); i -= i3[2] * (im[0] * im[1]);
      i3[1] = i / im[0]; i-= i3[1] * im[0];
      i3[0] = i;
      if (!(i3[2] >= 1 && i3[2] <= sp->g->nz)) {
	mprintf("i3 %d %d %d\n", i3[0], i3[1], i3[2]);
	assert(0);
      }
#endif
      struct vpic_mparticles_prt prt;
      prt.dx[0] = p->dx;
      prt.dx[1] = p->dy;
      prt.dx[2] = p->dz;
      prt.i     = p->i;
      prt.ux[0] = p->ux;
      prt.ux[1] = p->uy;
      prt.ux[2] = p->uz;
      prt.w     = p->w;
      prt.kind  = sp->id;
      getter(&prt, n - off);
    }

    v_off += v_n_prts;
  }
}

template<typename F>
void vpic_mparticles_set_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   F setter)
{
  unsigned int v_off = 0;
  for (auto sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    unsigned int v_n_prts = sp->np;

    unsigned int nb = std::max(v_off, off), ne = std::min(v_off + v_n_prts, off + n_prts);
    for (unsigned int n = nb; n < ne; n++) {
      struct vpic_mparticles_prt prt;
      setter(&prt, n - off);
      Particles::Particle *p = &sp->p[n - v_off];
      p->dx = prt.dx[0];
      p->dy = prt.dx[1];
      p->dz = prt.dx[2];
      p->i  = prt.i;
      p->ux = prt.ux[0];
      p->uy = prt.ux[1];
      p->uz = prt.ux[2];
      p->w  = prt.w;
      assert(prt.kind == sp->id);
    }

    v_off += v_n_prts;
  }
}



#endif
