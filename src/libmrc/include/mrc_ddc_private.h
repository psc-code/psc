
#ifndef MRC_DDC_PRIVATE_H
#define MRC_DDC_PRIVATE_H

#include "mrc_ddc.h"

#define N_DIR (27)

struct mrc_ddc_dir {
  int dir[3];
  int recv_ib[3], recv_ie[3];
  int send_ib[3], send_ie[3];
  int in[3], n_bnd; // len of bnd box (x,y,z), total
  int recv_rank;
  int send_rank;
  float *recv_buf;
  float *send_buf;
};

struct mrc_ddc {
  struct mrc_ddc_ops *ops;
  MPI_Comm comm;
  struct mrc_ddc_params params;
  int rank;  // rank (in comm)
  int size;  // size (in comm), must == np[0]*np[1]*np[2]
  int ip[3]; // rank split in terms of crd direction
  struct mrc_ddc_dir dir[N_DIR];
  MPI_Request recv_reqs[N_DIR], send_reqs[N_DIR];
  int n_recv_reqs, n_send_reqs;
};

#endif
