
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
	  
	  nei->send_buf_size = 8;
	  nei->send_buf = malloc(nei->send_buf_size * size_of_particle);
	  nei->n_send = 0;

	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    // use this one as buffer for particles that stay in the same patch
	    nei->rank = -1;
	  } else {
	    mrc_ddc_get_nei_rank_patch(ddc, p, dir, &nei->rank, &nei->patch);
	  }
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

// ----------------------------------------------------------------------
// ddc_particles_comm
//
// OPT: could use MPI_Waitany?
// OPT: overall more async
// OPT: 1d instead of 3d loops
// OPT: make the status buffers only as large as needed?

struct info_by_rank {
  struct send_entry {
    int patch; // source patch (source rank is this rank)
    int nei_patch; // target patch (target rank is index in send_entry)
    int dir1;  // direction
    int dir1neg;
    int n_send;
  } *send_entry;
  int n_send_entries;
  int n_send;

  struct recv_entry { // needs to be same as send_entry with different order!
    int nei_patch;
    int patch;
    int dir1neg;
    int dir1;
    int n_recv;
  } *recv_entry;
  int n_recv_entries;
  int n_recv;

  struct recv_entry *recv_entry_;
};

void
ddc_particles_comm(struct ddc_particles *ddcp, void *particles)
{
  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // FIXME, this is assuming our struct is equiv to an array of real_type
  assert(ddcp->size_of_particle % ddcp->size_of_real == 0);
  int sz = ddcp->size_of_particle / ddcp->size_of_real;
  int dir[3];

  // OPT, could be prepared once
  struct info_by_rank *info = calloc(size, sizeof(*info));

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  info[nei->rank].n_recv_entries++;
	}
      }
    }
  }
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      info[r].recv_entry =
	malloc(info[r].n_recv_entries * sizeof(*info[r].recv_entry));
      info[r].n_recv_entries = 0;
    }
  }

#if 1
  // post receives for # particles we'll receive in the next step
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
#if 0
	  mprintf(":%d irecv # from %d:%d tag %d\n", p, nei->rank, nei->patch,
		  dir1neg + p * N_DIR);
#endif
	  MPI_Irecv(&nei->n_recv, 1, MPI_INT, nei->rank, dir1neg + p * N_DIR,
		    comm, &ddcp->recv_reqs[dir1 + p * N_DIR]);
	}
      }
    }
  }
#endif

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  info[nei->rank].n_send_entries++;
	}
      }
    }
  }
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      info[r].send_entry =
	malloc(info[r].n_send_entries * sizeof(*info[r].send_entry));
      info[r].n_send_entries = 0;
    }
  }
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  struct send_entry *se =
	    &info[nei->rank].send_entry[info[nei->rank].n_send_entries++];
	  se->patch = p;
	  se->nei_patch = nei->patch;
	  se->dir1 = dir1;
	  se->dir1neg = dir1neg;
	  se->n_send = nei->n_send;
	  info[nei->rank].n_send += nei->n_send;
	}
      }
    }
  }

#if 1
  // post sends for # particles
  // to remote procs only
  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_send_entries; i++) {
      struct send_entry *se = &info[r].send_entry[i];
#if 0
      mprintf("S %d %d:%d -> %d:%d dir1 %02d n_send %d\n", i, rank,
	      se->patch, r, se->nei_patch, se->dir1, se->n_send);
#endif
      MPI_Isend(&se->n_send, 1, MPI_INT,
		r, se->dir1 + se->nei_patch * N_DIR, comm,
		&ddcp->send_reqs[se->dir1 + se->patch * N_DIR]);
    }
  }
#endif
  
#if 1
  // finish receiving remote # of particles
  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
#endif

  // fill recv_entry
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  struct recv_entry *re =
	    &info[nei->rank].recv_entry[info[nei->rank].n_recv_entries++];
	  re->patch = p;
	  re->nei_patch = nei->patch;
	  re->dir1 = dir1;
	  re->dir1neg = dir1neg;
	  re->n_recv = nei->n_recv;
	}
      }
    }
  }

  int n_recv_ranks = 0, n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      n_recv_ranks++;
    }
    if (info[r].n_send_entries) {
      n_send_ranks++;
    }
  }

  MPI_Request *send_req = malloc(n_send_ranks * sizeof(*send_req));
  MPI_Request *recv_req = malloc(n_recv_ranks * sizeof(*recv_req));

  n_recv_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      info[r].recv_entry_ =
	malloc(info[r].n_recv_entries * sizeof(struct recv_entry));
      MPI_Irecv(info[r].recv_entry_,
		sizeof(struct recv_entry) / sizeof(int) * info[r].n_recv_entries,
		MPI_INT, r, 0, comm, &recv_req[n_recv_ranks++]);
    }
  }  

  n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      MPI_Isend(info[r].send_entry,
		sizeof(struct send_entry) / sizeof(int) * info[r].n_send_entries,
		MPI_INT, r, 0, comm, &send_req[n_send_ranks++]);
    }
  }  

  MPI_Waitall(n_recv_ranks, recv_req, MPI_STATUSES_IGNORE);
  MPI_Waitall(n_send_ranks, send_req, MPI_STATUSES_IGNORE);

  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct recv_entry *re = &info[r].recv_entry[i];
      /* mprintf("R %d %d:%d -> %d:%d dir1 %02d n_recv %d\n", i, r, */
      /* 	      re->nei_patch, rank, re->patch, re->dir1, re->n_recv); */
    }
  }

  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct recv_entry *re = &info[r].recv_entry_[i];
      /* mprintf("r %d %d:%d -> %d:%d dir1 %02d n_recv %d\n", i, r, */
      /* 	      re->nei_patch, rank, re->patch, re->dir1, re->n_recv); */
    }
  }

  // set n_recv
  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      info[r].recv_entry[i] = info[r].recv_entry_[i];
      struct recv_entry *re = &info[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      struct ddcp_nei *nei = &patch->nei[re->dir1];
      nei->n_recv = re->n_recv;
      info[r].n_recv += re->n_recv;
    }
  }

  // post sends for the actual particles
  // to remote procs only
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
#if 0
	  mprintf(":%d isend P to   %d:%d tag %d\n", p, nei->rank, nei->patch,
		  dir1 + nei->patch * N_DIR);
#endif
	  MPI_Isend(nei->send_buf, sz * nei->n_send, ddcp->mpi_type_real,
		    nei->rank, dir1 + nei->patch * N_DIR, comm,
		    &ddcp->sendp_reqs[dir1 + p * N_DIR]);
	}
      }
    }
  }

  // count local # particles
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    patch->n_recv = 0;
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank != rank) {
	    continue;
	  }
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei_send = &ddcp->patches[nei->patch].nei[dir1neg];
	  patch->n_recv += nei_send->n_send;
	}
      }
    }
    mprintf("lp %d %d\n", p, patch->n_recv);
  }

  // add remote # particles
  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct recv_entry *re = &info[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      patch->n_recv += re->n_recv;
    }
  }

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    ddcp->realloc(particles, p, patch->head + patch->n_recv);
  }

  // post particle receives
  // from remote procs only
  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct recv_entry *re = &info[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
#if 0
      mprintf(":%d irecv P from %d:%d tag %d len %d\n", re->patch, r, re->nei_patch,
	      re->dir1neg + re->patch * N_DIR, re->n_recv);
#endif
      void *addr = ddcp->get_addr(particles, re->patch, patch->head);
      if (re->n_recv)
	mprintf("A addr %d:%d %d\n", r, re->patch, patch->head);
      MPI_Irecv(addr, sz * re->n_recv, ddcp->mpi_type_real,
		r, re->dir1neg + re->patch * N_DIR, comm,
		&ddcp->recv_reqs[re->dir1 + re->patch * N_DIR]);
      patch->head += re->n_recv;
    }
  }

  // copy particles from local proc
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank != rank) {
	    continue;
	  }
	  void *addr = ddcp->get_addr(particles, p, patch->head);
	  struct ddcp_nei *nei_send = &ddcp->patches[nei->patch].nei[dir1neg];
	  memcpy(addr, nei_send->send_buf, nei_send->n_send * ddcp->size_of_particle);
	  patch->head += nei_send->n_send;
	}
      }
    }
  }

  // wait for receive of remote particles and clean up
  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->nr_patches * N_DIR, ddcp->sendp_reqs, MPI_STATUSES_IGNORE);

  int n_send = 0, n_recv = 0;
  for (int r = 0; r < size; r++) {
    n_send += info[r].n_send;
    /* if (info[r].n_send) { */
    /*   mprintf("n_send[%d] = %d\n", r, info[r].n_send); */
    /* } */
    n_recv += info[r].n_recv;
    /* if (info[r].n_recv) { */
    /*   mprintf("n_recv[%d] = %d\n", r, info[r].n_recv); */
    /* } */
  }
  
  void *send_buf = malloc(n_send * ddcp->size_of_particle);
  void *p = send_buf;
  n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_send == 0)
      continue;

    void *p0 = p;
    for (int i = 0; i < info[r].n_send_entries; i++) {
      struct send_entry *se = &info[r].send_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[se->patch];
      struct ddcp_nei *nei = &patch->nei[se->dir1];
      memcpy(p, nei->send_buf, se->n_send * ddcp->size_of_particle);
      p += se->n_send * ddcp->size_of_particle;
    }
    MPI_Isend(p0, sz * info[r].n_send, ddcp->mpi_type_real,
	      r, 1, comm, &send_req[n_send_ranks++]);
  }
  assert(p == send_buf + n_send * ddcp->size_of_particle);

  void *recv_buf = malloc(n_recv * ddcp->size_of_particle);
  p = recv_buf;
  n_recv_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv == 0)
      continue;

    MPI_Irecv(p, sz * info[r].n_recv, ddcp->mpi_type_real,
	      r, 1, comm, &recv_req[n_recv_ranks++]);
    p += info[r].n_recv * ddcp->size_of_particle;
  }
  assert(p == recv_buf + n_recv * ddcp->size_of_particle);

  MPI_Waitall(n_recv_ranks, recv_req, MPI_STATUSES_IGNORE);
  MPI_Waitall(n_send_ranks, send_req, MPI_STATUSES_IGNORE);

  int *old_head = malloc(ddcp->nr_patches * sizeof(*old_head));
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    //    mprintf("A head %d: %d\n", p, patch->head);
    old_head[p] = patch->head;
    patch->head -= patch->n_recv;
  }

  p = recv_buf;
  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct recv_entry *re = &info[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      void *addr = ddcp->get_addr(particles, re->patch, patch->head);
      if (re->n_recv)
	mprintf("B addr %d:%d %d\n", r, re->patch, patch->head);
      assert(memcmp(addr, p, re->n_recv * ddcp->size_of_particle) == 0);
      patch->head += re->n_recv;
      p += re->n_recv * ddcp->size_of_particle;
    }
  }
  assert(p == recv_buf + n_recv * ddcp->size_of_particle);

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    patch->head = old_head[p];
  }
  free(old_head);

  free(send_buf);
  free(recv_buf);
  
  free(recv_req);
  free(send_req);

  for (int r = 0; r < size; r++) {
    free(info[r].send_entry);
    free(info[r].recv_entry);
    free(info[r].recv_entry_);
  }
  free(info);
}

