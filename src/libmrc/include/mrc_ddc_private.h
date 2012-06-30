
#ifndef MRC_DDC_PRIVATE_H
#define MRC_DDC_PRIVATE_H

#include "mrc_ddc.h"

#define N_DIR (27)

struct mrc_ddc_sendrecv_entry {
  int patch; // patch on this rank
  int nei_patch; // partner patch (partner rank is index in send/recv_entry)
  int dir1;  // direction
  int len;
  int ilo[3];
  int ihi[3];
};

struct mrc_ddc_rank_info {
  // what to send, by rank
  struct mrc_ddc_sendrecv_entry *send_entry;
  int n_send_entries;
  int n_send;

  // what to receive, by rank
  struct mrc_ddc_sendrecv_entry *recv_entry;
  int n_recv_entries;
  int n_recv;

  // for setting up the recv_entry's in the wrong order first
  struct mrc_ddc_sendrecv_entry *recv_entry_;
};

struct mrc_ddc_pattern2 {
  // communication info for each rank (NULL for those we don't communicate with)
  struct mrc_ddc_rank_info *ri;
  // number of ranks we're communicating with (excluding self)
  int n_recv_ranks, n_send_ranks;
  // one request each per rank we're communicating with
  MPI_Request *send_req, *recv_req;
  // total number of (ddc->mpi_type) we're sending / receivng to all ranks
  int n_send, n_recv;
  // buffers with the above sizes
  void *send_buf, *recv_buf;
  void *local_buf;
};

struct mrc_ddc_sendrecv {
  int ilo[3], ihi[3];
  int nei_rank;
  int nei_patch;
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
  MRC_SUBCLASS_OPS(struct mrc_ddc);
  void (*set_domain)(struct mrc_ddc *ddc, struct mrc_domain *domain);
  struct mrc_domain *(*get_domain)(struct mrc_ddc *ddc);
  void (*fill_ghosts)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  void (*add_ghosts)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  void (*get_nei_rank_patch)(struct mrc_ddc *ddc, int p, int dir[3],
			     int *nei_rank, int *nei_patch);
};

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

extern struct mrc_ddc_ops mrc_ddc_simple_ops;

// ======================================================================
// mrc_ddc_multi

struct mrc_ddc_multi {
  struct mrc_domain *domain;
  int np[3]; // # patches per direction
  int bc[3]; // boundary condition
  int nr_patches;
  struct mrc_patch *patches;
  int mpi_rank, mpi_size;
  struct mrc_ddc_pattern2 add_ghosts2;
  struct mrc_ddc_pattern2 fill_ghosts2;
};

extern struct mrc_ddc_ops mrc_ddc_multi_ops;

// ======================================================================
// mrc_ddc_amr

extern struct mrc_ddc_ops mrc_ddc_amr_ops;

#endif
