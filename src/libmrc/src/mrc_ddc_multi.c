
#include "mrc_ddc_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_bits.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

// ======================================================================
// mrc_ddc_multi

// we dynamically create patterns for mrc_fld exchange up to this many layers of
// ghost points


#define mrc_ddc_multi(ddc) mrc_to_subobj(ddc, struct mrc_ddc_multi)

// ----------------------------------------------------------------------
// mrc_ddc_multi_get_nei_rank_patch

static void
mrc_ddc_multi_get_nei_rank_patch(struct mrc_ddc *ddc, int p, int dir[3],
				 int *nei_rank, int *nei_patch)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(sub->domain, p, &info);
  int patch_idx_nei[3];
  for (int d = 0; d < 3; d++) {
    patch_idx_nei[d] = info.idx3[d] + dir[d];
    if (sub->bc[d] == BC_PERIODIC) {
      if (patch_idx_nei[d] < 0) {
	patch_idx_nei[d] += sub->np[d];
      }
      if (patch_idx_nei[d] >= sub->np[d]) {
	patch_idx_nei[d] -= sub->np[d];
      }
    }
    if (patch_idx_nei[d] < 0 || patch_idx_nei[d] >= sub->np[d]) {
      *nei_rank = -1;
      *nei_patch = -1;
      return;
    }
  }
  mrc_domain_get_level_idx3_patch_info(sub->domain, 0, patch_idx_nei, &info);
  *nei_rank = info.rank;
  *nei_patch = info.patch;
}

// ----------------------------------------------------------------------
// ddc_init_outside

static void
ddc_init_outside(struct mrc_ddc *ddc, int p, struct mrc_ddc_sendrecv *sr,
		 int dir[3], int ibn[3])
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
    sr->nei_rank = -1;
    return;
  }

  mrc_ddc_multi_get_nei_rank_patch(ddc, p, dir, &sr->nei_rank, &sr->nei_patch);
  if (sr->nei_rank < 0)
    return;

  sr->len = 1;
  int ilo[3], ihi[3];
  for (int d = 0; d < 3; d++) {
    ilo[d] = 0;
    ihi[d] = sub->patches[p].ldims[d];
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = ilo[d] - ibn[d];
      sr->ihi[d] = ilo[d];
      break;
    case 0:
      sr->ilo[d] = ilo[d];
      sr->ihi[d] = ihi[d];
      break;
    case 1:
      sr->ilo[d] = ihi[d];
      sr->ihi[d] = ihi[d] + ibn[d];
      break;
    }
    sr->len *= (sr->ihi[d] - sr->ilo[d]);
  }
}

// ----------------------------------------------------------------------
// ddc_init_inside

static void
ddc_init_inside(struct mrc_ddc *ddc, int p, struct mrc_ddc_sendrecv *sr,
		int dir[3], int ibn[3])
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
    sr->nei_rank = -1;
    return;
  }

  mrc_ddc_multi_get_nei_rank_patch(ddc, p, dir, &sr->nei_rank, &sr->nei_patch);
  if (sr->nei_rank < 0)
    return;

  sr->len = 1;
  int ilo[3], ihi[3];
  for (int d = 0; d < 3; d++) {
    ilo[d] = 0;
    ihi[d] = sub->patches[p].ldims[d];
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = ilo[d];
      sr->ihi[d] = ilo[d] + ibn[d];
      break;
    case 0:
      sr->ilo[d] = ilo[d];
      sr->ihi[d] = ihi[d];
      break;
    case 1:
      sr->ilo[d] = ihi[d] - ibn[d];
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
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  sub->domain = domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_get_domain

static struct mrc_domain *
mrc_ddc_multi_get_domain(struct mrc_ddc *ddc)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  return sub->domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_free_buffers

static void
mrc_ddc_multi_free_buffers(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2)
{
  free(patt2->send_buf);
  free(patt2->recv_buf);
  free(patt2->local_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_alloc_buffers

static void
mrc_ddc_multi_alloc_buffers(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
			    int n_fields)
{
  if (ddc->size_of_type > patt2->max_size_of_type ||
      n_fields > patt2->max_n_fields) {

    if (ddc->size_of_type > patt2->max_size_of_type) {
      patt2->max_size_of_type = ddc->size_of_type;
    }
    if (n_fields > patt2->max_n_fields) {
      patt2->max_n_fields = n_fields;
    }

    mrc_ddc_multi_free_buffers(ddc, patt2);

    patt2->recv_buf = malloc(patt2->n_recv * patt2->max_n_fields * patt2->max_size_of_type);
    patt2->send_buf = malloc(patt2->n_send * patt2->max_n_fields * patt2->max_size_of_type);
    patt2->local_buf = malloc(patt2->local_buf_size * patt2->max_n_fields * patt2->max_size_of_type);
  }
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_setup_pattern2

static void
mrc_ddc_multi_setup_pattern2(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
			     void (*init_send)(struct mrc_ddc *, int p,
					       struct mrc_ddc_sendrecv *, int [3], int[3]),
			     void (*init_recv)(struct mrc_ddc *, int p,
					       struct mrc_ddc_sendrecv *, int [3], int[3]),
			     int ibn[3])
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  struct mrc_ddc_rank_info *ri = calloc(sub->mpi_size, sizeof(*ri));
  patt2->ri = ri;

  int dir[3];

  // count how many {send,recv}_entries per rank
  for (int p = 0; p < sub->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  struct mrc_ddc_sendrecv sr;
	  init_send(ddc, p, &sr, dir, ibn);
	  if (sr.nei_rank >= 0) {
	    ri[sr.nei_rank].n_send_entries++;
	  }
	  init_recv(ddc, p, &sr, dir, ibn);
	  if (sr.nei_rank >= 0) {
	    ri[sr.nei_rank].n_recv_entries++;
	  }
	}
      }
    }
  }

  // alloc {send,recv}_entries
  for (int r = 0; r < sub->mpi_size; r++) {
    if (ri[r].n_send_entries) {
      ri[r].send_entry = malloc(ri[r].n_send_entries * sizeof(*ri[r].send_entry));
      ri[r].n_send_entries = 0;
      if (r != sub->mpi_rank) {
	patt2->n_send_ranks++;
      }
    }
    if (ri[r].n_recv_entries) {
      ri[r].recv_entry = malloc(ri[r].n_recv_entries * sizeof(*ri[r].recv_entry));
      ri[r].recv_entry_ = malloc(ri[r].n_recv_entries * sizeof(*ri[r].recv_entry_));
      ri[r].n_recv_entries = 0;
      if (r != sub->mpi_rank) {
	patt2->n_recv_ranks++;
      }
    }
  }

  // set up {send,recv}_entries per rank
  for (int p = 0; p < sub->nr_patches; p++) {
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);

	  struct mrc_ddc_sendrecv sr;
	  init_send(ddc, p, &sr, dir, ibn);
	  if (sr.nei_rank >= 0) {
	    struct mrc_ddc_sendrecv_entry *se =
	      &ri[sr.nei_rank].send_entry[ri[sr.nei_rank].n_send_entries++];
	    se->patch = p;
	    se->nei_patch = sr.nei_patch;
	    se->len = sr.len;
	    se->dir1 = dir1;
	    memcpy(se->ilo, sr.ilo, 3 * sizeof(*se->ilo));
	    memcpy(se->ihi, sr.ihi, 3 * sizeof(*se->ihi));
	    ri[sr.nei_rank].n_send += sr.len;
	    if (sr.nei_rank != sub->mpi_rank) {
	      patt2->n_send += sr.len;
	    }
	  }

	  init_recv(ddc, p, &sr, dir, ibn);
	  if (sr.nei_rank >= 0) {
	    struct mrc_ddc_sendrecv_entry *re =
	      &ri[sr.nei_rank].recv_entry_[ri[sr.nei_rank].n_recv_entries++];
	    re->patch = p;
	    re->nei_patch = sr.nei_patch;
	    re->len = sr.len;
	    re->dir1 = dir1;
	    memcpy(re->ilo, sr.ilo, 3 * sizeof(*re->ilo));
	    memcpy(re->ihi, sr.ihi, 3 * sizeof(*re->ihi));
	    ri[sr.nei_rank].n_recv += sr.len;
	    if (sr.nei_rank != sub->mpi_rank) {
	      patt2->n_recv += sr.len;
	    }
	  }
	}
      }
    }
  }

  // reorganize recv_entries into right order (same as send on the sending rank)
  for (int r = 0; r < sub->mpi_size; r++) {
    for (int p = 0, cnt = 0; cnt < ri[r].n_recv_entries; p++) {
      for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
	for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	  for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	    int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });

	    for (int i = 0; i < ri[r].n_recv_entries; i++) {
	      struct mrc_ddc_sendrecv_entry *re = &ri[r].recv_entry_[i];
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
  for (int i = 0; i < ri[sub->mpi_rank].n_recv_entries; i++) {
    struct mrc_ddc_sendrecv_entry *re = &ri[sub->mpi_rank].recv_entry[i];
    local_buf_size = MAX(local_buf_size,
			 (re->ihi[0] - re->ilo[0]) * 
			 (re->ihi[1] - re->ilo[1]) * 
			 (re->ihi[2] - re->ilo[2]));
  }

  patt2->local_buf_size = local_buf_size;
  patt2->send_req = malloc(patt2->n_send_ranks * sizeof(*patt2->send_req));
  patt2->recv_req = malloc(patt2->n_recv_ranks * sizeof(*patt2->recv_req));
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_destroy_rank_info

static void
mrc_ddc_multi_destroy_pattern2(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  free(patt2->send_req);
  free(patt2->recv_req);

  for (int r = 0; r < sub->mpi_size; r++) {
    free(patt2->ri[r].send_entry);
    free(patt2->ri[r].recv_entry);
  }
  free(patt2->ri);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_set_mpi_type
//
// FIXME, this could be done more nicely for mrc_fld interface by having
// mrc_fld know what MPI type to use, or at least use MRC_NT_*

static void
mrc_ddc_multi_set_mpi_type(struct mrc_ddc *ddc)
{
  if (ddc->size_of_type == sizeof(float)) {
    ddc->mpi_type = MPI_FLOAT;
  } else if (ddc->size_of_type == sizeof(double)) {
    ddc->mpi_type = MPI_DOUBLE;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_setup

static void
mrc_ddc_multi_setup(struct mrc_ddc *ddc)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  assert(sub->domain);
  mrc_domain_get_nr_procs(sub->domain, sub->np);
  mrc_domain_get_bc(sub->domain, sub->bc);

  sub->patches = mrc_domain_get_patches(sub->domain,
					&sub->nr_patches);

  MPI_Comm_rank(mrc_ddc_comm(ddc), &sub->mpi_rank);
  MPI_Comm_size(mrc_ddc_comm(ddc), &sub->mpi_size);

  mrc_ddc_multi_setup_pattern2(ddc, &sub->fill_ghosts2,
			       ddc_init_inside, ddc_init_outside, ddc->ibn);
  mrc_ddc_multi_setup_pattern2(ddc, &sub->add_ghosts2,
			       ddc_init_outside, ddc_init_inside, ddc->ibn);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_destroy

static void
mrc_ddc_multi_destroy(struct mrc_ddc *ddc)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  mrc_ddc_multi_free_buffers(ddc, &sub->fill_ghosts2);
  mrc_ddc_multi_free_buffers(ddc, &sub->add_ghosts2);

  mrc_ddc_multi_destroy_pattern2(ddc, &sub->fill_ghosts2);
  mrc_ddc_multi_destroy_pattern2(ddc, &sub->add_ghosts2);

  for (int nr_ghosts = 1; nr_ghosts <= MAX_NR_GHOSTS; nr_ghosts++) {
    if (sub->fill_ghosts[nr_ghosts]) {
      mrc_ddc_multi_free_buffers(ddc, sub->fill_ghosts[nr_ghosts]);
      mrc_ddc_multi_destroy_pattern2(ddc, sub->fill_ghosts[nr_ghosts]);
    }
  }
}

// ----------------------------------------------------------------------
// ddc_run_begin

static void
ddc_run_begin(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
	      int mb, int me, void *ctx,
	      void (*to_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
  struct mrc_ddc_rank_info *ri = patt2->ri;

  // communicate aggregated buffers
  // post receives
  patt2->recv_cnt = 0;
  void *p = patt2->recv_buf;
  for (int r = 0; r < sub->mpi_size; r++) {
    if (r != sub->mpi_rank && ri[r].n_recv_entries) {
      MPI_Irecv(p, ri[r].n_recv * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &patt2->recv_req[patt2->recv_cnt++]);
      p += ri[r].n_recv * (me - mb) * ddc->size_of_type;
    }
  }  
  assert(p == patt2->recv_buf + patt2->n_recv * (me - mb) * ddc->size_of_type);

  // post sends
  patt2->send_cnt = 0;
  p = patt2->send_buf;
  for (int r = 0; r < sub->mpi_size; r++) {
    if (r != sub->mpi_rank && ri[r].n_send_entries) {
      void *p0 = p;
      for (int i = 0; i < ri[r].n_send_entries; i++) {
	struct mrc_ddc_sendrecv_entry *se = &ri[r].send_entry[i];
	to_buf(mb, me, se->patch, se->ilo, se->ihi, p, ctx);
	p += se->len * (me - mb) * ddc->size_of_type;
      }
      MPI_Isend(p0, ri[r].n_send * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &patt2->send_req[patt2->send_cnt++]);
    }
  }  
  assert(p == patt2->send_buf + patt2->n_send * (me - mb) * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// ddc_run_end

static void
ddc_run_end(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
	    int mb, int me, void *ctx,
	    void (*from_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
  struct mrc_ddc_rank_info *ri = patt2->ri;

  MPI_Waitall(patt2->recv_cnt, patt2->recv_req, MPI_STATUSES_IGNORE);

  void *p = patt2->recv_buf;
  for (int r = 0; r < sub->mpi_size; r++) {
    if (r != sub->mpi_rank) {
      for (int i = 0; i < ri[r].n_recv_entries; i++) {
	struct mrc_ddc_sendrecv_entry *re = &ri[r].recv_entry[i];
	from_buf(mb, me, re->patch, re->ilo, re->ihi, p, ctx);
	p += re->len * (me - mb) * ddc->size_of_type;
      }
    }
  }

  MPI_Waitall(patt2->send_cnt, patt2->send_req, MPI_STATUSES_IGNORE);
}

// ----------------------------------------------------------------------
// ddc_run_local

static void
ddc_run_local(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
	      int mb, int me, void *ctx,
	      void (*to_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx),
	      void (*from_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
  struct mrc_ddc_rank_info *ri = patt2->ri;

  // overlap: local exchange
  for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
    struct mrc_ddc_sendrecv_entry *se = &ri[sub->mpi_rank].send_entry[i];
    struct mrc_ddc_sendrecv_entry *re = &ri[sub->mpi_rank].recv_entry[i];
    if (se->ilo[0] == se->ihi[0] ||
	se->ilo[1] == se->ihi[1] ||
	se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
      continue;
    }
    to_buf(mb, me, se->patch, se->ilo, se->ihi, patt2->local_buf, ctx);
    from_buf(mb, me, se->nei_patch, re->ilo, re->ihi, patt2->local_buf, ctx);
  }
}

// ----------------------------------------------------------------------
// ddc_run

static void
ddc_run(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
	int mb, int me, void *ctx,
	void (*to_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx),
	void (*from_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  ddc_run_begin(ddc, patt2, mb, me, ctx, to_buf);
  ddc_run_local(ddc, patt2, mb, me, ctx, to_buf, from_buf);
  ddc_run_end(ddc, patt2, mb, me, ctx, from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_add_ghosts

static void
mrc_ddc_multi_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  mrc_ddc_multi_set_mpi_type(ddc);
  mrc_ddc_multi_alloc_buffers(ddc, &sub->add_ghosts2, me - mb);
  ddc_run(ddc, &sub->add_ghosts2, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->add_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_fill_ghosts

static void
mrc_ddc_multi_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  mrc_ddc_multi_set_mpi_type(ddc);
  mrc_ddc_multi_alloc_buffers(ddc, &sub->fill_ghosts2, me - mb);
  ddc_run(ddc, &sub->fill_ghosts2, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->copy_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_fill_ghosts_begin

static void
mrc_ddc_multi_fill_ghosts_begin(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
 
  mrc_ddc_multi_set_mpi_type(ddc);
  mrc_ddc_multi_alloc_buffers(ddc, &sub->fill_ghosts2, me - mb);
  ddc_run_begin(ddc, &sub->fill_ghosts2, mb, me, ctx,
		ddc->funcs->copy_to_buf);
}

static void
mrc_ddc_multi_fill_ghosts_end(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  ddc_run_end(ddc, &sub->fill_ghosts2, mb, me, ctx,
	      ddc->funcs->copy_from_buf);
}

static void
mrc_ddc_multi_fill_ghosts_local(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  ddc_run_local(ddc, &sub->fill_ghosts2, mb, me, ctx,
		ddc->funcs->copy_to_buf, ddc->funcs->copy_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_fill_ghosts_fld

static void
mrc_ddc_multi_fill_ghosts_fld(struct mrc_ddc *ddc, int mb, int me,
			      struct mrc_fld *fld)
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);

  int nr_ghosts = fld->_nr_ghosts;
  assert(nr_ghosts > 0 && nr_ghosts <= MAX_NR_GHOSTS);
  if (!sub->fill_ghosts[nr_ghosts]) {
    sub->fill_ghosts[nr_ghosts] = calloc(1, sizeof(*sub->fill_ghosts[0]));
    mrc_ddc_multi_setup_pattern2(ddc, sub->fill_ghosts[nr_ghosts],
				 ddc_init_inside, ddc_init_outside,
				 fld->_sw.vals);
  }

  ddc->size_of_type = fld->_nd->size_of_type;
  mrc_ddc_multi_set_mpi_type(ddc);
  mrc_ddc_multi_alloc_buffers(ddc, sub->fill_ghosts[nr_ghosts], me - mb);
  ddc_run(ddc, sub->fill_ghosts[nr_ghosts], mb, me, fld,
	  mrc_fld_ddc_copy_to_buf, mrc_fld_ddc_copy_from_buf);
}

// ======================================================================
// mrc_ddc_multi_ops

#define VAR(x) (void *)offsetof(struct mrc_ddc_multi, x)
static struct param mrc_ddc_multi_descr[] = {
  { "domain"          , VAR(domain)       , PARAM_OBJ(mrc_domain) },
  {},
};
#undef VAR

struct mrc_ddc_ops mrc_ddc_multi_ops = {
  .name                  = "multi",
  .size                  = sizeof(struct mrc_ddc_multi),
  .param_descr           = mrc_ddc_multi_descr,
  .setup                 = mrc_ddc_multi_setup,
  .destroy               = mrc_ddc_multi_destroy,
  .set_domain            = mrc_ddc_multi_set_domain,
  .get_domain            = mrc_ddc_multi_get_domain,
  .fill_ghosts_fld       = mrc_ddc_multi_fill_ghosts_fld,
  .fill_ghosts           = mrc_ddc_multi_fill_ghosts,
  .fill_ghosts_begin     = mrc_ddc_multi_fill_ghosts_begin,
  .fill_ghosts_end       = mrc_ddc_multi_fill_ghosts_end,
  .fill_ghosts_local     = mrc_ddc_multi_fill_ghosts_local,
  .add_ghosts            = mrc_ddc_multi_add_ghosts,
};

