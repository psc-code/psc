
#ifndef MRC_DDC_PRIVATE_H
#define MRC_DDC_PRIVATE_H

#include "mrc_ddc.h"

#define N_DIR (27)

struct mrc_ddc_sendrecv {
  int ilo[3], ihi[3];
  int rank_nei;
  int len;
  void *buf;
};

struct mrc_ddc_pattern {
  struct mrc_ddc_sendrecv send[N_DIR];
  struct mrc_ddc_sendrecv recv[N_DIR];
};

struct mrc_ddc {
  struct mrc_obj obj;
  // parameters
  int size_of_type;
  int max_n_fields;
  int ibn[3]; // # ghost points

  int rank, size;
  struct mrc_ddc_funcs *funcs;
  MPI_Datatype mpi_type;
};

struct mrc_ddc_ops {
  MRC_OBJ_OPS;
  void (*set_domain)(struct mrc_ddc *ddc, struct mrc_domain *domain);
  void (*fill_ghosts)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  void (*add_ghosts)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  int  (*get_rank_nei)(struct mrc_ddc *ddc, int dir[3]);
};

void libmrc_ddc_register(struct mrc_ddc_ops *ops);

#define to_mrc_ddc(obj) container_of(obj, struct mrc_ddc, obj)

// ======================================================================
// mrc_ddc_simple

struct mrc_ddc_simple {
  // parameters
  int n_proc[3]; // # procs in 3D grid
  int ilo[3], ihi[3]; // local domain (no ghosts)
  int bc[3]; // boundary condition

  int proc[3]; // this proc's position in the 3D proc grid
  struct mrc_ddc_pattern add_ghosts;
  struct mrc_ddc_pattern fill_ghosts;
  MPI_Request send_reqs[N_DIR];
  MPI_Request recv_reqs[N_DIR];
};

void libmrc_ddc_register_simple(void);

// ======================================================================
// mrc_ddc_multi

struct mrc_ddc_multi {
  struct mrc_domain *domain;
  int np[3]; // # patches per direction
  int bc[3]; // boundary condition
  int nr_patches;
  struct mrc_patch *patches;
  int proc[3]; // this proc's position in the 3D proc grid
  struct mrc_ddc_pattern *add_ghosts;
  struct mrc_ddc_pattern *fill_ghosts;
  MPI_Request send_reqs[N_DIR];
  MPI_Request recv_reqs[N_DIR];
};

void libmrc_ddc_register_multi(void);


#endif
