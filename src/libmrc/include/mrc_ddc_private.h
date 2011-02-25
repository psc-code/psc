
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
  MPI_Comm comm;
  struct mrc_ddc_params prm;
  struct mrc_ddc_ops *ops;
  int rank, size;
  int proc[3]; // this proc's position in the 3D proc grid
  MPI_Datatype mpi_type;
  struct mrc_ddc_pattern add_ghosts;
  struct mrc_ddc_pattern fill_ghosts;
  MPI_Request send_reqs[N_DIR];
  MPI_Request recv_reqs[N_DIR];
};

#endif
