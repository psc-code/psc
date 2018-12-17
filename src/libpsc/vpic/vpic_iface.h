
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

#include "vpic_config.h"

#include <stdbool.h>

#include <mpi.h>
#include <mrc_common.h>

// ----------------------------------------------------------------------
// vpic_mparticles

struct psc_particle_inject;

// ----------------------------------------------------------------------
// Simulation

struct vpic_simulation_info;
struct field_array;

// Harris specific
void Simulation_set_region_resistive_harris(struct vpic_harris_params *prm,
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

struct vpic_simulation;

// FIXME, replicated
#define BOUNDARY(i,j,k) (13+(i)+3*(j)+9*(k)) /* FORTRAN -1:1,-1:1,-1:1 */

#ifndef mprintf

#include <mpi.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

#endif


#endif
