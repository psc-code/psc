
#include "mrc_ddc_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define to_mrc_ddc_multi(ddc) ((struct mrc_ddc_multi *) (ddc)->obj.subctx)

// ----------------------------------------------------------------------
// mrc_ddc_multi_get_nei_rank_patch

static void
mrc_ddc_multi_get_nei_rank_patch(struct mrc_ddc *ddc, int p, int dir[3],
				 int *nei_rank, int *nei_patch)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  struct mrc_ddc_patch *ddc_patch = &multi->ddc_patches[p];
  int patch_idx_nei[3];
  for (int d = 0; d < 3; d++) {
    patch_idx_nei[d] = ddc_patch->patch_idx[d] + dir[d];
    if (multi->bc[d] == BC_PERIODIC) {
      if (patch_idx_nei[d] < 0) {
	patch_idx_nei[d] += multi->np[d];
      }
      if (patch_idx_nei[d] >= multi->np[d]) {
	patch_idx_nei[d] -= multi->np[d];
      }
    }
    if (patch_idx_nei[d] < 0 || patch_idx_nei[d] >= multi->np[d]) {
      *nei_rank = -1;
      *nei_patch = -1;
      return;
    }
  }
  struct mrc_patch_info info;
  mrc_domain_get_idx3_patch_info(multi->domain, patch_idx_nei, &info);
  *nei_rank = info.rank;
  *nei_patch = info.patch;
}

// ----------------------------------------------------------------------
// ddc_init_outside

static void
ddc_init_outside(struct mrc_ddc *ddc, int p, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  mrc_ddc_multi_get_nei_rank_patch(ddc, p, dir, &sr->nei_rank, &sr->nei_patch);
  if (sr->nei_rank < 0)
    return;

  sr->len = 1;
  int ilo[3], ihi[3];
  for (int d = 0; d < 3; d++) {
    ilo[d] = 0;
    ihi[d] = multi->patches[p].ldims[d];
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = ilo[d] - ddc->ibn[d];
      sr->ihi[d] = ilo[d];
      break;
    case 0:
      sr->ilo[d] = ilo[d];
      sr->ihi[d] = ihi[d];
      break;
    case 1:
      sr->ilo[d] = ihi[d];
      sr->ihi[d] = ihi[d] + ddc->ibn[d];
      break;
    }
    sr->len *= (sr->ihi[d] - sr->ilo[d]);
  }
  sr->buf = malloc(sr->len * ddc->max_n_fields * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// ddc_init_inside

static void
ddc_init_inside(struct mrc_ddc *ddc, int p, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  mrc_ddc_multi_get_nei_rank_patch(ddc, p, dir, &sr->nei_rank, &sr->nei_patch);
  if (sr->nei_rank < 0)
    return;

  sr->len = 1;
  int ilo[3], ihi[3];
  for (int d = 0; d < 3; d++) {
    ilo[d] = 0;
    ihi[d] = multi->patches[p].ldims[d];
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = ilo[d];
      sr->ihi[d] = ilo[d] + ddc->ibn[d];
      break;
    case 0:
      sr->ilo[d] = ilo[d];
      sr->ihi[d] = ihi[d];
      break;
    case 1:
      sr->ilo[d] = ihi[d] - ddc->ibn[d];
      sr->ihi[d] = ihi[d];
      break;
    }
    sr->len *= (sr->ihi[d] - sr->ilo[d]);
  }
  sr->buf = malloc(sr->len * ddc->max_n_fields * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_set_domain

static void
mrc_ddc_multi_set_domain(struct mrc_ddc *ddc, struct mrc_domain *domain)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  multi->domain = domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_get_domain

static struct mrc_domain *
mrc_ddc_multi_get_domain(struct mrc_ddc *ddc)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  return multi->domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_setup

static void
mrc_ddc_multi_setup(struct mrc_ddc *ddc)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  if (ddc->size_of_type == sizeof(float)) {
    ddc->mpi_type = MPI_FLOAT;
  } else if (ddc->size_of_type == sizeof(double)) {
    ddc->mpi_type = MPI_DOUBLE;
  } else {
    assert(0);
  }

  assert(ddc->max_n_fields > 0);
  assert(multi->domain);
  mrc_domain_get_nr_procs(multi->domain, multi->np);
  mrc_domain_get_bc(multi->domain, multi->bc);

  multi->patches = mrc_domain_get_patches(multi->domain,
					  &multi->nr_patches);
  multi->add_ghosts = calloc(multi->nr_patches, sizeof(*multi->add_ghosts));
  multi->fill_ghosts = calloc(multi->nr_patches, sizeof(*multi->fill_ghosts));
  multi->ddc_patches = calloc(multi->nr_patches, sizeof(*multi->ddc_patches));
  for (int p = 0; p < multi->nr_patches; p++) {
    struct mrc_ddc_patch *ddc_patch = &multi->ddc_patches[p];
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(multi->domain, p, &info);
    for (int d = 0; d < 3; d++) {
      ddc_patch->patch_idx[d] = info.idx3[d];
    }

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
	    continue;
	  
	  int dir1 = mrc_ddc_dir2idx(dir);

	  struct mrc_ddc_pattern *add_ghosts = &multi->add_ghosts[p];
	  ddc_init_outside(ddc, p, &add_ghosts->send[dir1], dir);
	  ddc_init_inside(ddc, p, &add_ghosts->recv[dir1], dir);
	  
	  struct mrc_ddc_pattern *fill_ghosts = &multi->fill_ghosts[p];
	  ddc_init_inside(ddc, p, &fill_ghosts->send[dir1], dir);
	  ddc_init_outside(ddc, p, &fill_ghosts->recv[dir1], dir);
	}
      }
    }
  }

  int rank;
  MPI_Comm_rank(mrc_ddc_comm(ddc), &rank);
  MPI_Comm_size(mrc_ddc_comm(ddc), &multi->mpi_size);
  multi->rank_info = calloc(multi->mpi_size, sizeof(*multi->rank_info));

  struct mrc_ddc_rank_info *ri = multi->rank_info;
  struct mrc_ddc_pattern *patt = multi->fill_ghosts;

  multi->n_recv_ranks = 0;
  multi->n_send_ranks = 0;
  multi->n_send = 0;
  multi->n_recv = 0;
  int dir[3];

  // count how many recv_entries per rank
  for (int p = 0; p < multi->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct mrc_ddc_sendrecv *r = &patt[p].recv[dir1];
	  if (r->nei_rank != rank && r->len > 0) {
	    ri[r->nei_rank].n_recv_entries++;
	  }
	}
      }
    }
  }

  // alloc recv_entries
  for (int r = 0; r < multi->mpi_size; r++) {
    if (ri[r].n_recv_entries) {
      ri[r].recv_entry = malloc(ri[r].n_recv_entries * sizeof(*ri[r].recv_entry));
      multi->n_recv_ranks++;
      ri[r].n_recv_entries = 0;
    }
  }

  // set up recv_entries by rank
  for (int p = 0; p < multi->nr_patches; p++) {
    for (dir[2] = 1; dir[2] >= -1; dir[2]--) {
      for (dir[1] = 1; dir[1] >= -1; dir[1]--) {
	for (dir[0] = 1; dir[0] >= -1; dir[0]--) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct mrc_ddc_sendrecv *r = &patt[p].recv[dir1];

	  if (r->nei_rank != rank && r->len > 0) {
	    struct mrc_ddc_recv_entry *re =
	      &ri[r->nei_rank].recv_entry[ri[r->nei_rank].n_recv_entries++];
	    re->patch = p;
	    re->nei_patch = r->nei_patch;
	    re->n_recv = r->len;
	    re->dir1 = dir1;
	    re->dir1neg = dir1neg;
	    ri[r->nei_rank].n_recv += r->len;
	    multi->n_recv += r->len;
	  }
	}
      }
    }
  }

  // count how many send_entries per rank
  for (int p = 0; p < multi->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct mrc_ddc_sendrecv *s = &patt[p].send[dir1];
	  if (s->nei_rank != rank && s->len > 0) {
	    ri[s->nei_rank].n_send_entries++;
	  }
	}
      }
    }
  }

  // alloc send_entries
  for (int r = 0; r < multi->mpi_size; r++) {
    if (ri[r].n_send_entries) {
      ri[r].send_entry = malloc(ri[r].n_send_entries * sizeof(*ri[r].send_entry));
      multi->n_send_ranks++;
      ri[r].n_send_entries = 0;
    }
  }

  // set up send_entries per rank
  for (int p = 0; p < multi->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct mrc_ddc_sendrecv *s = &patt[p].send[dir1];
	  if (s->nei_rank != rank && s->len > 0) {
	    struct mrc_ddc_send_entry *se =
	      &ri[s->nei_rank].send_entry[ri[s->nei_rank].n_send_entries++];
	    se->patch = p;
	    se->nei_patch = s->nei_patch;
	    se->n_send = s->len;
	    se->dir1 = dir1;
	    se->dir1neg = dir1neg;
	    ri[s->nei_rank].n_send += s->len;
	    multi->n_send += s->len;
	  }
	}
      }
    }
  }

  multi->send_req = malloc(multi->n_send_ranks * sizeof(*multi->send_req));
  multi->recv_req = malloc(multi->n_recv_ranks * sizeof(*multi->recv_req));

  multi->n_recv_ranks = 0;
  for (int r = 0; r < multi->mpi_size; r++) {
    if (ri[r].n_recv_entries) {
      ri[r].recv_entry_ =
	malloc(ri[r].n_recv_entries * sizeof(*ri[r].recv_entry));
      MPI_Irecv(ri[r].recv_entry_,
		sizeof(struct mrc_ddc_recv_entry) / sizeof(int) * ri[r].n_recv_entries,
		MPI_INT, r, 0, ddc->obj.comm, &multi->recv_req[multi->n_recv_ranks++]);
    }
  }  

  multi->n_send_ranks = 0;
  for (int r = 0; r < multi->mpi_size; r++) {
    if (ri[r].n_send_entries) {
      MPI_Isend(ri[r].send_entry,
		sizeof(struct mrc_ddc_send_entry) / sizeof(int) * ri[r].n_send_entries,
		MPI_INT, r, 0, ddc->obj.comm, &multi->send_req[multi->n_send_ranks++]);
    }
  }  

  MPI_Waitall(multi->n_recv_ranks, multi->recv_req, MPI_STATUSES_IGNORE);
  MPI_Waitall(multi->n_send_ranks, multi->send_req, MPI_STATUSES_IGNORE);

#if 0
  for (int rr = 0; rr < size; rr++) {
    MPI_Barrier(ddc->obj.comm);
    if (rank == rr) {
      for (int r = 0; r < size; r++) {
	for (int i = 0; i < ri[r].n_recv_entries; i++) {
	  struct recv_entry *re = &ri[r].recv_entry[i];
	  mprintf("R %d:%d -> %d:%d dir %02d len %2d tag %d\n",
		  r, re->nei_patch, rank, re->patch, re->dir1, re->n_recv,
		  re->nei_patch * N_DIR + re->dir1);
	}
	for (int i = 0; i < ri[r].n_recv_entries; i++) {
	  struct recv_entry *re = &ri[r].recv_entry_[i];
	  mprintf("r %d:%d -> %d:%d dir %02d len %2d tag %d\n",
		  r, re->nei_patch, rank, re->patch, re->dir1, re->n_recv,
		  re->nei_patch * N_DIR + re->dir1);
	}
	mprintf("=====\n");
      }
      MPI_Barrier(ddc->obj.comm);
    }
  }
#endif

  // use received recv_entries rather than calculated ones
  for (int r = 0; r < multi->mpi_size; r++) {
    for (int i = 0; i < ri[r].n_recv_entries; i++) {
      ri[r].recv_entry[i] = ri[r].recv_entry_[i];
    }
  }

  multi->recv_buf = malloc(multi->n_recv * ddc->max_n_fields * ddc->size_of_type);
  multi->send_buf = malloc(multi->n_send * ddc->max_n_fields * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_destroy

static void
mrc_ddc_multi_destroy(struct mrc_ddc *ddc)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  for (int p = 0; p < multi->nr_patches; p++) {
    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
	    continue;
	  
	  int dir1 = mrc_ddc_dir2idx(dir);

	  struct mrc_ddc_pattern *add_ghosts = &multi->add_ghosts[p];
	  free(add_ghosts->send[dir1].buf);
	  free(add_ghosts->recv[dir1].buf);
	  
	  struct mrc_ddc_pattern *fill_ghosts = &multi->fill_ghosts[p];
	  free(fill_ghosts->send[dir1].buf);
	  free(fill_ghosts->recv[dir1].buf);
	}
      }
    }
  }
  free(multi->add_ghosts);
  free(multi->fill_ghosts);
  free(multi->ddc_patches);

  free(multi->send_req);
  free(multi->recv_req);

  free(multi->send_buf);
  free(multi->recv_buf);

  for (int r = 0; r < multi->mpi_size; r++) {
    free(multi->rank_info[r].send_entry);
    free(multi->rank_info[r].recv_entry);
    free(multi->rank_info[r].recv_entry_);
  }
  free(multi->rank_info);
}

// ----------------------------------------------------------------------
// ddc_run

static void
ddc_run(struct mrc_ddc *ddc, struct mrc_ddc_pattern *patt, int mb, int me,
	void *ctx,
	void (*to_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx),
	void (*from_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);
  int rank;
  MPI_Comm_rank(mrc_ddc_comm(ddc), &rank);
  struct mrc_ddc_rank_info *ri = multi->rank_info;

  // communicate aggregated buffers
  multi->n_recv_ranks = 0;
  void *p = multi->recv_buf;
  for (int r = 0; r < multi->mpi_size; r++) {
    if (ri[r].n_recv_entries) {
      MPI_Irecv(p, ri[r].n_recv * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &multi->recv_req[multi->n_recv_ranks++]);
      p += ri[r].n_recv * (me - mb) * ddc->size_of_type;
    }
  }  
  assert(p == multi->recv_buf + multi->n_recv * (me - mb) * ddc->size_of_type);
  
  p = multi->send_buf;
  multi->n_send_ranks = 0;
  for (int r = 0; r < multi->mpi_size; r++) {
    if (ri[r].n_send_entries) {
      void *p0 = p;
      for (int i = 0; i < ri[r].n_send_entries; i++) {
	struct mrc_ddc_send_entry *se = &ri[r].send_entry[i];
	struct mrc_ddc_sendrecv *s = &patt[se->patch].send[se->dir1];
	to_buf(mb, me, se->patch, s->ilo, s->ihi, p, ctx);
	p += se->n_send * (me - mb) * ddc->size_of_type;
      }
      MPI_Isend(p0, ri[r].n_send * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &multi->send_req[multi->n_send_ranks++]);
    }
  }  
  assert(p == multi->send_buf + multi->n_send * (me - mb) * ddc->size_of_type);

  // overlap: local exchange
  int dir[3];
  for (int p = 0; p < multi->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
 	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct mrc_ddc_sendrecv *s = &patt[p].send[dir1];
	  if (s->nei_rank == rank && s->len > 0) {
	    struct mrc_ddc_sendrecv *r = &patt[s->nei_patch].recv[dir1neg];
	    // OPT, we could just do the copy without going through
	    // a buffer
	    to_buf(mb, me, p, s->ilo, s->ihi, s->buf, ctx);
	    from_buf(mb, me, s->nei_patch, r->ilo, r->ihi, s->buf, ctx);
	  }
	}
      }
    }
  }

  MPI_Waitall(multi->n_recv_ranks, multi->recv_req, MPI_STATUSES_IGNORE);

  p = multi->recv_buf;
  for (int r = 0; r < multi->mpi_size; r++) {
    for (int i = 0; i < ri[r].n_recv_entries; i++) {
      struct mrc_ddc_recv_entry *re = &ri[r].recv_entry[i];
      struct mrc_ddc_sendrecv *rcv = &patt[re->patch].recv[re->dir1];
      from_buf(mb, me, re->patch, rcv->ilo, rcv->ihi, p, ctx);
      p += re->n_recv * (me - mb) * ddc->size_of_type;
    }
  }

  MPI_Waitall(multi->n_send_ranks, multi->send_req, MPI_STATUSES_IGNORE);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_add_ghosts

static void
mrc_ddc_multi_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  ddc_run(ddc, multi->add_ghosts, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->add_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_fill_ghosts

static void
mrc_ddc_multi_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  ddc_run(ddc, multi->fill_ghosts, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->copy_from_buf);
}

// ======================================================================
// mrc_ddc_multi_ops

struct mrc_ddc_ops mrc_ddc_multi_ops = {
  .name                  = "multi",
  .size                  = sizeof(struct mrc_ddc_multi),
  .setup                 = mrc_ddc_multi_setup,
  .destroy               = mrc_ddc_multi_destroy,
  .set_domain            = mrc_ddc_multi_set_domain,
  .get_domain            = mrc_ddc_multi_get_domain,
  .fill_ghosts           = mrc_ddc_multi_fill_ghosts,
  .add_ghosts            = mrc_ddc_multi_add_ghosts,
  .get_nei_rank_patch    = mrc_ddc_multi_get_nei_rank_patch,
};

// ======================================================================
// mrc_ddc_funcs_m3 for mrc_m3

#include <mrc_fld.h>

static void
mrc_m3_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
		   void *_buf, void *ctx)
{
  //  mprintf("to %d:%d x %d:%d x %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
  struct mrc_m3 *m3 = ctx;
  float *buf = _buf;

  struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf,m - mb, ix,iy,iz) = MRC_M3(m3p,m, ix,iy,iz);
	}
      }
    }
  }
}

static void
mrc_m3_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
		     void *_buf, void *ctx)
{
  //  mprintf("from %d:%d x %d:%d x %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
  struct mrc_m3 *m3 = ctx;
  float *buf = _buf;

  struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_M3(m3p,m, ix,iy,iz) = MRC_DDC_BUF3(buf,m - mb, ix,iy,iz);
	}
      }
    }
  }
}

struct mrc_ddc_funcs mrc_ddc_funcs_m3 = {
  .copy_to_buf   = mrc_m3_copy_to_buf,
  .copy_from_buf = mrc_m3_copy_from_buf,
};

