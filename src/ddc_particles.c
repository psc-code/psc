
#include "ddc_particles.h"

#include <mrc_ddc.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================
// ddc_particles

struct ddc_particles *
ddc_particles_create(struct mrc_ddc *ddc, int size_of_particle,
		     int size_of_real, MPI_Datatype mpi_type_real,
		     void  (*realloc)(void *, int, int),
		     void *(*get_addr)(void *, int, int))
{
  struct ddc_particles *ddcp = malloc(sizeof(*ddcp));
  memset(ddcp, 0, sizeof(*ddcp));

  struct mrc_domain *domain = mrc_ddc_get_domain(ddc);
  mrc_domain_get_patches(domain, &ddcp->nr_patches);
  ddcp->size_of_particle = size_of_particle;
  ddcp->size_of_real = size_of_real;
  ddcp->mpi_type_real = mpi_type_real;
  ddcp->realloc = realloc;
  ddcp->get_addr = get_addr;
  ddcp->patches = calloc(ddcp->nr_patches, sizeof(*ddcp->patches));
  ddcp->send_reqs  = calloc(ddcp->nr_patches * N_DIR, sizeof(*ddcp->send_reqs));
  ddcp->sendp_reqs = calloc(ddcp->nr_patches * N_DIR, sizeof(*ddcp->sendp_reqs));
  ddcp->recv_reqs  = calloc(ddcp->nr_patches * N_DIR, sizeof(*ddcp->recv_reqs));
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  ddcp->send_reqs[dir1 + p * N_DIR] = MPI_REQUEST_NULL;
	  ddcp->sendp_reqs[dir1 + p * N_DIR] = MPI_REQUEST_NULL;
	  ddcp->recv_reqs[dir1 + p * N_DIR] = MPI_REQUEST_NULL;
	  
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    nei->rank = -1;
	    continue;
	  }
	  
	  nei->send_buf_size = 8;
	  nei->send_buf = malloc(nei->send_buf_size * size_of_particle);
	  nei->n_send = 0;
	  mrc_ddc_get_nei_rank_patch(ddc, p, dir, &nei->rank, &nei->patch);
	}
      }
    }
  }

  return ddcp;
}

void
ddc_particles_destroy(struct ddc_particles *ddcp)
{
  if (!ddcp)
    return;

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  
	  free(nei->send_buf);
	}
      }
    }
  }
  free(ddcp->patches);
  free(ddcp->send_reqs);
  free(ddcp->sendp_reqs);
  free(ddcp->recv_reqs);
  free(ddcp);
}

void
ddc_particles_queue(struct ddc_particles *ddcp, struct ddcp_patch *patch,
		    int dir[3], void *p)
{
  struct ddcp_nei *nei = &patch->nei[mrc_ddc_dir2idx(dir)];

  if (nei->n_send == nei->send_buf_size) {
    // reallocate a larger buffer, doubling buffer size each time
    assert(nei->send_buf_size > 0);
    nei->send_buf_size *= 2;
    nei->send_buf = realloc(nei->send_buf, 
			    nei->send_buf_size * ddcp->size_of_particle);
  }
  memcpy(nei->send_buf + nei->n_send * ddcp->size_of_particle, p,
	 ddcp->size_of_particle);
  nei->n_send++;
}

void
ddc_particles_comm(struct ddc_particles *ddcp, void *particles)
{
  // FIXME, this is assuming our struct is equiv to an array of real_type
  assert(ddcp->size_of_particle % ddcp->size_of_real == 0);
  int sz = ddcp->size_of_particle / ddcp->size_of_real;
  int dir[3];

  // post receives for # particles we'll receive in the next step
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
#if 0
	  mprintf(":%d irecv # from %d:%d tag %d\n", p, nei->rank, nei->patch,
		  dir1neg + p * N_DIR);
#endif
	  MPI_Irecv(&nei->n_recv, 1, MPI_INT, nei->rank, dir1neg + p * N_DIR,
		    MPI_COMM_WORLD, &ddcp->recv_reqs[dir1 + p * N_DIR]);
	}
      }
    }
  }

  // post sends for # particles and then the actual particles
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
#if 0
	  mprintf(":%d isend # to   %d:%d tag %d: len %d\n", p, nei->rank, nei->patch,
		  dir1 + nei->patch * N_DIR, nei->n_send);
	  mprintf(":%d isend P to   %d:%d tag %d\n", p, nei->rank, nei->patch,
		  dir1 + nei->patch * N_DIR);
#endif
	  MPI_Isend(&nei->n_send, 1, MPI_INT,
		    nei->rank, dir1 + nei->patch * N_DIR, MPI_COMM_WORLD,
		    &ddcp->send_reqs[dir1 + p * N_DIR]);
	  MPI_Isend(nei->send_buf, sz * nei->n_send, ddcp->mpi_type_real,
		    nei->rank, dir1 + nei->patch * N_DIR, MPI_COMM_WORLD,
		    &ddcp->sendp_reqs[dir1 + p * N_DIR]);
	}
      }
    }
  }

  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->recv_reqs, MPI_STATUSES_IGNORE);

  // calc total # of particles
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    int new_n_particles = patch->head; // particles which stayed on this proc

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
	  // add the ones we're about to receive
	  new_n_particles += nei->n_recv;
	}
      }
    }
    ddcp->realloc(particles, p, new_n_particles);
  }

  // post particle receives
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
#if 0
	  mprintf(":%d irecv P from %d:%d tag %d len %d\n", p, nei->rank, nei->patch,
		  dir1neg + p * N_DIR, nei->n_recv);
#endif
	  void *addr = ddcp->get_addr(particles, p, patch->head);
	  MPI_Irecv(addr, sz * nei->n_recv, ddcp->mpi_type_real,
		    nei->rank, dir1neg + p * N_DIR, MPI_COMM_WORLD,
		    &ddcp->recv_reqs[dir1 + p * N_DIR]);
	  patch->head += nei->n_recv;
	}
      }
    }
  }

  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->sendp_reqs, MPI_STATUSES_IGNORE);
}

