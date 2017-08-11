
#ifndef DDC_PARTICLES_H
#define DDC_PARTICLES_H

#include "psc.h"

#define N_DIR (27)

struct ddcp_info_by_rank {
  struct ddcp_send_entry {
    int patch; // source patch (source rank is this rank)
    int nei_patch; // target patch (target rank is index in send_entry)
    int dir1;  // direction
    int dir1neg;
  } *send_entry;
  int *send_cnts;
  int n_send_entries;
  int n_send;

  struct ddcp_recv_entry { // needs to be same as send_entry with different order!
    int nei_patch;
    int patch;
    int dir1neg;
    int dir1;
  } *recv_entry;
  int *recv_cnts;
  int n_recv_entries;
  int n_recv;

  int rank;
};

struct ddcp_nei {
  void *send_buf;
  int n_send;
  int n_recv;
  int rank;
  int patch;
  int send_buf_size;
};

struct ddcp_patch {
  int head;
  struct ddcp_nei nei[N_DIR];
  int n_recv;
};

struct ddc_particles {
  int nr_patches;
  struct ddcp_patch *patches;
  int size_of_particle;
  MPI_Datatype mpi_type_real;
  void  (*realloc)(void *mparticles, int p, int new_nr_particles);
  void *(*get_addr)(void *mparticles, int p, int n);
  struct ddcp_info_by_rank *by_rank;
  struct ddcp_info_by_rank *cinfo; // compressed info
  int n_ranks;
  MPI_Request *send_reqs;
  MPI_Request *recv_reqs;

  struct mrc_domain *domain;
};

#endif
