
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
    int n_send;
  } *send_entry;
  int n_send_entries;
  int n_send;

  struct ddcp_recv_entry { // needs to be same as send_entry with different order!
    int nei_patch;
    int patch;
    int dir1neg;
    int dir1;
    int n_recv;
  } *recv_entry;
  int n_recv_entries;
  int n_recv;

  struct ddcp_recv_entry *recv_entry_;
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
  int size_of_real;
  MPI_Datatype mpi_type_real;
  void  (*realloc)(void *mparticles, int p, int new_nr_particles);
  void *(*get_addr)(void *mparticles, int p, int n);
  struct ddcp_info_by_rank *by_rank;
  int n_send_ranks;
  int n_recv_ranks;
  MPI_Request *send_reqs;
  MPI_Request *recv_reqs;
};

struct ddc_particles *ddc_particles_create(struct mrc_ddc *ddc, int size_of_particle,
					   int size_of_real, MPI_Datatype mpi_type_real,
					   void (*realloc)(void *, int, int),
					   void *(*get_addr)(void *, int, int));
void ddc_particles_destroy(struct ddc_particles *ddcp);
void ddc_particles_queue(struct ddc_particles *ddcp, struct ddcp_patch *patch,
			 int dir[3], void *p);
void ddc_particles_comm(struct ddc_particles *ddcp, void *particles);


#endif
