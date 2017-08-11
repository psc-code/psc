
#include "ddc_particles.h"

#include <mrc_ddc.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================
// ddc_particles

struct ddc_particles *
ddc_particles_create(struct mrc_domain *domain, int size_of_particle,
		     int size_of_real, MPI_Datatype mpi_type_real,
		     void  (*realloc)(void *, int, int),
		     void *(*get_addr)(void *, int, int))
{
  struct ddc_particles *ddcp = malloc(sizeof(*ddcp));
  memset(ddcp, 0, sizeof(*ddcp));

  ddcp->domain = domain;
  mrc_domain_get_patches(domain, &ddcp->nr_patches);
  ddcp->size_of_particle = size_of_particle;
  ddcp->size_of_real = size_of_real;
  ddcp->mpi_type_real = mpi_type_real;
  ddcp->realloc = realloc;
  ddcp->get_addr = get_addr;
  ddcp->patches = calloc(ddcp->nr_patches, sizeof(*ddcp->patches));
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  
	  nei->send_buf_size = 8;
	  nei->send_buf = malloc(nei->send_buf_size * size_of_particle);
	  nei->n_send = 0;

	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    // use this one as buffer for particles that stay in the same patch
	    nei->rank = -1;
	  } else {
	    mrc_domain_get_neighbor_rank_patch(domain, p, dir, &nei->rank, &nei->patch);
	  }
	}
      }
    }
  }

  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  ddcp->by_rank = calloc(size, sizeof(*ddcp->by_rank));
  struct ddcp_info_by_rank *info = ddcp->by_rank;

  int dir[3];

  // count how many recv_entries per rank
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

  // alloc recv_entries
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      info[r].recv_entry =
	malloc(info[r].n_recv_entries * sizeof(*info[r].recv_entry));
      info[r].recv_cnts =
	malloc(info[r].n_recv_entries * sizeof(*info[r].recv_cnts));
    }
  }

#if 0
  // set up recv_entries
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
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct recv_entry *re =
	    &info[nei->rank].recv_entry[info[nei->rank].n_recv_entries++];
	  re->patch = p;
	  re->nei_patch = nei->patch;
	  re->dir1 = dir1;
	  re->dir1neg = dir1neg;
	}
      }
    }
  }
#endif

  // count send_entries
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

  // alloc send_entries
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      info[r].send_entry =
	malloc(info[r].n_send_entries * sizeof(*info[r].send_entry));
      info[r].send_cnts = 
	malloc(info[r].n_send_entries * sizeof(*info[r].send_cnts));
      info[r].n_send_entries = 0;
    }
  }

  // set up send_entries
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
	  struct ddcp_send_entry *se =
	    &info[nei->rank].send_entry[info[nei->rank].n_send_entries++];
	  se->patch = p;
	  se->nei_patch = nei->patch;
	  se->dir1 = dir1;
	  se->dir1neg = dir1neg;
	}
      }
    }
  }

  int n_recv_ranks = 0;
  int n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      n_recv_ranks++;
    }
    if (info[r].n_send_entries) {
      n_send_ranks++;
    }
  }

  assert(n_send_ranks == n_recv_ranks);
  ddcp->n_ranks = n_send_ranks;

  ddcp->send_reqs = malloc(ddcp->n_ranks * sizeof(*ddcp->send_reqs));
  ddcp->recv_reqs = malloc(ddcp->n_ranks * sizeof(*ddcp->recv_reqs));

  n_recv_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      MPI_Irecv(info[r].recv_entry,
		sizeof(struct ddcp_recv_entry) / sizeof(int) * info[r].n_recv_entries,
		MPI_INT, r, 111, comm, &ddcp->recv_reqs[n_recv_ranks++]);
    }
  }  

  n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      MPI_Isend(info[r].send_entry,
		sizeof(struct ddcp_send_entry) / sizeof(int) * info[r].n_send_entries,
		MPI_INT, r, 111, comm, &ddcp->send_reqs[n_send_ranks++]);
    }
  }  

  // FIXME / OPT, we're copying alloc'd pointers over,
  // fragile, though correct. info could be free'd here or even
  // earlier
  ddcp->cinfo = malloc(ddcp->n_ranks * sizeof(*ddcp->cinfo));
  struct ddcp_info_by_rank *cinfo = ddcp->cinfo;
  int i = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      assert(info[r].n_send_entries);
      cinfo[i] = info[r];
      cinfo[i].rank = r;
      i++;
    }
  }
  assert(i == ddcp->n_ranks);

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
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];

	  free(nei->send_buf);
	}
      }
    }
  }
  free(ddcp->patches);

  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int size;
  MPI_Comm_size(comm, &size);

  struct ddcp_info_by_rank *info = ddcp->by_rank;
  for (int r = 0; r < size; r++) {
    free(info[r].send_entry);
    free(info[r].recv_entry);
    free(info[r].send_cnts);
    free(info[r].recv_cnts);
  }
  free(ddcp->by_rank);
  free(ddcp->cinfo);
  free(ddcp->recv_reqs);
  free(ddcp->send_reqs);

  free(ddcp);
}

// ----------------------------------------------------------------------
// ddc_particles_comm
//
// OPT: could use MPI_Waitany?
// OPT: overall more async
// OPT: 1d instead of 3d loops
// OPT: make the status buffers only as large as needed?

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

  struct ddcp_info_by_rank *cinfo = ddcp->cinfo;

  for (int r = 0; r < ddcp->n_ranks; r++) {
    MPI_Irecv(cinfo[r].recv_cnts, cinfo[r].n_recv_entries,
	      MPI_INT, cinfo[r].rank, 222, comm, &ddcp->recv_reqs[r]);
  }  

  for (int r = 0; r < ddcp->n_ranks; r++) {
    cinfo[r].n_send = 0;
    for (int i = 0; i < cinfo[r].n_send_entries; i++) {
      struct ddcp_send_entry *se = &cinfo[r].send_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[se->patch];
      struct ddcp_nei *nei = &patch->nei[se->dir1];
      cinfo[r].send_cnts[i] = nei->n_send;
      cinfo[r].n_send += nei->n_send;
    }
    MPI_Isend(cinfo[r].send_cnts, cinfo[r].n_send_entries,
	      MPI_INT, cinfo[r].rank, 222, comm, &ddcp->send_reqs[r]);
  }  

  // overlap: count local # particles
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
  }

#if 0
  // still want to figure out how to avoid sending info i can calc
  // just need to fix order
  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &info[r].recv_entry[i];
      mprintf("R %d %d:%d -> %d:%d dir1 %02d n_recv %d\n", i, r,
      	      re->nei_patch, rank, re->patch, re->dir1, re->n_recv);
    }
  }

  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &info[r].recv_entry_[i];
      mprintf("r %d %d:%d -> %d:%d dir1 %02d n_recv %d\n", i, r,
      	      re->nei_patch, rank, re->patch, re->dir1, re->n_recv);
    }
  }
#endif

  MPI_Waitall(ddcp->n_ranks, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->n_ranks, ddcp->send_reqs, MPI_STATUSES_IGNORE);

  // add remote # particles
  for (int r = 0; r < ddcp->n_ranks; r++) {
    cinfo[r].n_recv = 0;
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      patch->n_recv += cinfo[r].recv_cnts[i];
      cinfo[r].n_recv += cinfo[r].recv_cnts[i];
    }
  }

  // realloc
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    ddcp->realloc(particles, p, patch->head + patch->n_recv);
  }

  // leave room for receives (FIXME? just change order)
  for (int r = 0; r < ddcp->n_ranks; r++) {
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      patch->head += cinfo[r].recv_cnts[i];
    }
  }

  int n_send = 0, n_recv = 0;
  for (int r = 0; r < ddcp->n_ranks; r++) {
    n_send += cinfo[r].n_send;
    n_recv += cinfo[r].n_recv;
  }
  
  // post sends
  void *send_buf = malloc(n_send * ddcp->size_of_particle);
  void *p = send_buf;
  for (int r = 0; r < ddcp->n_ranks; r++) {
    if (cinfo[r].n_send == 0)
      continue;

    void *p0 = p;
    for (int i = 0; i < cinfo[r].n_send_entries; i++) {
      struct ddcp_send_entry *se = &cinfo[r].send_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[se->patch];
      struct ddcp_nei *nei = &patch->nei[se->dir1];
      memcpy(p, nei->send_buf, cinfo[r].send_cnts[i] * ddcp->size_of_particle);
      p += cinfo[r].send_cnts[i] * ddcp->size_of_particle;
    }
    MPI_Isend(p0, sz * cinfo[r].n_send, ddcp->mpi_type_real,
	      cinfo[r].rank, 1, comm, &ddcp->send_reqs[r]);
  }
  assert(p == send_buf + n_send * ddcp->size_of_particle);

  // post receives
  void *recv_buf = malloc(n_recv * ddcp->size_of_particle);
  p = recv_buf;
  for (int r = 0; r < ddcp->n_ranks; r++) {
    if (cinfo[r].n_recv == 0)
      continue;

    MPI_Irecv(p, sz * cinfo[r].n_recv, ddcp->mpi_type_real,
	      cinfo[r].rank, 1, comm, &ddcp->recv_reqs[r]);
    p += cinfo[r].n_recv * ddcp->size_of_particle;
  }
  assert(p == recv_buf + n_recv * ddcp->size_of_particle);

  // overlap: copy particles from local proc
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

  MPI_Waitall(ddcp->n_ranks, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->n_ranks, ddcp->send_reqs, MPI_STATUSES_IGNORE);

  // copy received particles into right place

  int *old_head = malloc(ddcp->nr_patches * sizeof(*old_head));
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    old_head[p] = patch->head;
    patch->head -= patch->n_recv;
  }

  p = recv_buf;
  for (int r = 0; r < ddcp->n_ranks; r++) {
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      void *addr = ddcp->get_addr(particles, re->patch, patch->head);
      memcpy(addr, p, cinfo[r].recv_cnts[i] * ddcp->size_of_particle);
      patch->head += cinfo[r].recv_cnts[i];
      p += cinfo[r].recv_cnts[i] * ddcp->size_of_particle;
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
}

