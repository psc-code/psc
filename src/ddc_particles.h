
#ifndef DDC_PARTICLES_H
#define DDC_PARTICLES_H

#include "psc.h"

#define N_DIR (27)

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
};

struct ddc_particles {
  struct ddcp_patch *patches;
  MPI_Request *send_reqs;
  MPI_Request *sendp_reqs;
  MPI_Request *recv_reqs;
  int size_of_particle;
};

struct ddc_particles *ddc_particles_create(struct mrc_ddc *ddc, int size_of_particle);
void ddc_particles_queue(struct ddc_particles *ddcp, struct ddcp_patch *patch,
			 int dir[3], void *p);
void ddc_particles_comm(struct ddc_particles *ddcp, mparticles_base_t *particles);


#endif
