
#include "mrc_ddc_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_bits.h>

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

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(multi->domain, p, &info);
  int patch_idx_nei[3];
  for (int d = 0; d < 3; d++) {
    patch_idx_nei[d] = info.idx3[d] + dir[d];
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
// mrc_ddc_multi_setup_pattern2

static void
mrc_ddc_multi_setup_pattern2(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
			     struct mrc_ddc_pattern *patt)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  struct mrc_ddc_rank_info *ri = calloc(multi->mpi_size, sizeof(*ri));
  patt2->ri = ri;

  int dir[3];

  // count how many {send,recv}_entries per rank
  for (int p = 0; p < multi->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct mrc_ddc_sendrecv *s = &patt[p].send[dir1];
	  if (s->len > 0) {
	    ri[s->nei_rank].n_send_entries++;
	  }
	  struct mrc_ddc_sendrecv *r = &patt[p].recv[dir1];
	  if (r->len > 0) {
	    ri[r->nei_rank].n_recv_entries++;
	  }
	}
      }
    }
  }

  // alloc {send,recv}_entries
  for (int r = 0; r < multi->mpi_size; r++) {
    if (ri[r].n_send_entries) {
      ri[r].send_entry = malloc(ri[r].n_send_entries * sizeof(*ri[r].send_entry));
      ri[r].n_send_entries = 0;
      if (r != multi->mpi_rank) {
	patt2->n_send_ranks++;
      }
    }
    if (ri[r].n_recv_entries) {
      ri[r].recv_entry = malloc(ri[r].n_recv_entries * sizeof(*ri[r].recv_entry));
      ri[r].recv_entry_ = malloc(ri[r].n_recv_entries * sizeof(*ri[r].recv_entry_));
      ri[r].n_recv_entries = 0;
      if (r != multi->mpi_rank) {
	patt2->n_recv_ranks++;
      }
    }
  }

  // set up {send,recv}_entries per rank
  for (int p = 0; p < multi->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);

	  struct mrc_ddc_sendrecv *s = &patt[p].send[dir1];
	  if (s->len > 0) {
	    struct mrc_ddc_send_entry *se =
	      &ri[s->nei_rank].send_entry[ri[s->nei_rank].n_send_entries++];
	    se->patch = p;
	    se->nei_patch = s->nei_patch;
	    se->n_send = s->len;
	    se->dir1 = dir1;
	    memcpy(se->ilo, s->ilo, 3 * sizeof(*se->ilo));
	    memcpy(se->ihi, s->ihi, 3 * sizeof(*se->ihi));
	    ri[s->nei_rank].n_send += s->len;
	    if (s->nei_rank != multi->mpi_rank) {
	      patt2->n_send += s->len;
	    }
	  }

	  struct mrc_ddc_sendrecv *r = &patt[p].recv[dir1];
	  if (r->len > 0) {
	    struct mrc_ddc_recv_entry *re =
	      &ri[r->nei_rank].recv_entry_[ri[r->nei_rank].n_recv_entries++];
	    re->patch = p;
	    re->nei_patch = r->nei_patch;
	    re->n_recv = r->len;
	    re->dir1 = dir1;
	    memcpy(re->ilo, r->ilo, 3 * sizeof(*re->ilo));
	    memcpy(re->ihi, r->ihi, 3 * sizeof(*re->ihi));
	    ri[r->nei_rank].n_recv += r->len;
	    if (r->nei_rank != multi->mpi_rank) {
	      patt2->n_recv += r->len;
	    }
	  }
	}
      }
    }
  }

  // reorganize recv_entries into right order (same as send on the sending rank)
  for (int r = 0; r < multi->mpi_size; r++) {
    for (int p = 0, cnt = 0; cnt < ri[r].n_recv_entries; p++) {
      for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
	for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	  for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	    int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });

	    for (int i = 0; i < ri[r].n_recv_entries; i++) {
	      struct mrc_ddc_recv_entry *re = &ri[r].recv_entry_[i];
	      if (re->nei_patch == p && re->dir1 == dir1neg) {
		ri[r].recv_entry[cnt++] = *re;
		break;
	      }
	    }
	  }
	}
      }
    }
    free(ri[r].recv_entry_);
  }

  // find buffer size for local exchange
  int local_buf_size = 0;
  for (int i = 0; i < ri[multi->mpi_rank].n_recv_entries; i++) {
    struct mrc_ddc_recv_entry *re = &ri[multi->mpi_rank].recv_entry[i];
    local_buf_size = MAX(local_buf_size,
			 (re->ihi[0] - re->ilo[0]) * 
			 (re->ihi[1] - re->ilo[1]) * 
			 (re->ihi[2] - re->ilo[2]));
  }

  patt2->send_req = malloc(patt2->n_send_ranks * sizeof(*patt2->send_req));
  patt2->recv_req = malloc(patt2->n_recv_ranks * sizeof(*patt2->recv_req));

  patt2->recv_buf = malloc(patt2->n_recv * ddc->max_n_fields * ddc->size_of_type);
  patt2->send_buf = malloc(patt2->n_send * ddc->max_n_fields * ddc->size_of_type);
  patt2->local_buf = malloc(local_buf_size * ddc->max_n_fields * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_destroy_rank_info

static void
mrc_ddc_multi_destroy_pattern2(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  free(patt2->send_req);
  free(patt2->recv_req);

  free(patt2->send_buf);
  free(patt2->recv_buf);
  free(patt2->local_buf);

  for (int r = 0; r < multi->mpi_size; r++) {
    free(patt2->ri[r].send_entry);
    free(patt2->ri[r].recv_entry);
  }
  free(patt2->ri);
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
  struct mrc_ddc_pattern *add_ghosts = calloc(multi->nr_patches, sizeof(*add_ghosts));
  struct mrc_ddc_pattern *fill_ghosts = calloc(multi->nr_patches, sizeof(*fill_ghosts));
  for (int p = 0; p < multi->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(multi->domain, p, &info);

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
	    continue;
	  
	  int dir1 = mrc_ddc_dir2idx(dir);

	  ddc_init_outside(ddc, p, &add_ghosts[p].send[dir1], dir);
	  ddc_init_inside(ddc, p, &add_ghosts[p].recv[dir1], dir);
	  
	  ddc_init_inside(ddc, p, &fill_ghosts[p].send[dir1], dir);
	  ddc_init_outside(ddc, p, &fill_ghosts[p].recv[dir1], dir);
	}
      }
    }
  }

  MPI_Comm_rank(mrc_ddc_comm(ddc), &multi->mpi_rank);
  MPI_Comm_size(mrc_ddc_comm(ddc), &multi->mpi_size);

  mrc_ddc_multi_setup_pattern2(ddc, &multi->fill_ghosts2, fill_ghosts);
  mrc_ddc_multi_setup_pattern2(ddc, &multi->add_ghosts2, add_ghosts);

  free(add_ghosts);
  free(fill_ghosts);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_destroy

static void
mrc_ddc_multi_destroy(struct mrc_ddc *ddc)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  mrc_ddc_multi_destroy_pattern2(ddc, &multi->fill_ghosts2);
  mrc_ddc_multi_destroy_pattern2(ddc, &multi->add_ghosts2);
}

// ----------------------------------------------------------------------
// ddc_run

static void
ddc_run(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
	int mb, int me, void *ctx,
	void (*to_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx),
	void (*from_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);
  struct mrc_ddc_rank_info *ri = patt2->ri;

  // communicate aggregated buffers
  int recv_cnt = 0;
  void *p = patt2->recv_buf;
  for (int r = 0; r < multi->mpi_size; r++) {
    if (r != multi->mpi_rank && ri[r].n_recv_entries) {
      MPI_Irecv(p, ri[r].n_recv * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &patt2->recv_req[recv_cnt++]);
      p += ri[r].n_recv * (me - mb) * ddc->size_of_type;
    }
  }  
  assert(p == patt2->recv_buf + patt2->n_recv * (me - mb) * ddc->size_of_type);
  
  p = patt2->send_buf;
  int send_cnt = 0;
  for (int r = 0; r < multi->mpi_size; r++) {
    if (r != multi->mpi_rank && ri[r].n_send_entries) {
      void *p0 = p;
      for (int i = 0; i < ri[r].n_send_entries; i++) {
	struct mrc_ddc_send_entry *se = &ri[r].send_entry[i];
	to_buf(mb, me, se->patch, se->ilo, se->ihi, p, ctx);
	p += se->n_send * (me - mb) * ddc->size_of_type;
      }
      MPI_Isend(p0, ri[r].n_send * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &patt2->send_req[send_cnt++]);
    }
  }  
  assert(p == patt2->send_buf + patt2->n_send * (me - mb) * ddc->size_of_type);

  // overlap: local exchange
  for (int i = 0; i < ri[multi->mpi_rank].n_send_entries; i++) {
    struct mrc_ddc_send_entry *se = &ri[multi->mpi_rank].send_entry[i];
    struct mrc_ddc_recv_entry *re = &ri[multi->mpi_rank].recv_entry[i];
    to_buf(mb, me, se->patch, se->ilo, se->ihi, patt2->local_buf, ctx);
    from_buf(mb, me, se->nei_patch, re->ilo, re->ihi, patt2->local_buf, ctx);
  }

  MPI_Waitall(recv_cnt, patt2->recv_req, MPI_STATUSES_IGNORE);

  p = patt2->recv_buf;
  for (int r = 0; r < multi->mpi_size; r++) {
    if (r != multi->mpi_rank) {
      for (int i = 0; i < ri[r].n_recv_entries; i++) {
	struct mrc_ddc_recv_entry *re = &ri[r].recv_entry[i];
	from_buf(mb, me, re->patch, re->ilo, re->ihi, p, ctx);
	p += re->n_recv * (me - mb) * ddc->size_of_type;
      }
    }
  }

  MPI_Waitall(send_cnt, patt2->send_req, MPI_STATUSES_IGNORE);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_add_ghosts

static void
mrc_ddc_multi_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  ddc_run(ddc, &multi->add_ghosts2, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->add_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_fill_ghosts

static void
mrc_ddc_multi_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  ddc_run(ddc, &multi->fill_ghosts2, mb, me, ctx,
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

