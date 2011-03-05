
#include "mrc_ddc_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>

#include <stdlib.h>
#include <assert.h>

#define to_mrc_ddc_multi(ddc) ((struct mrc_ddc_multi *) (ddc)->obj.subctx)

// ----------------------------------------------------------------------
// mrc_ddc_multi_get_nei_patch_info

static void
mrc_ddc_multi_nei_get_patch_info(struct mrc_ddc *ddc, int p, int dir[3],
				 struct mrc_patch_info *info)
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
      info->rank = -1;
      info->patch = -1;
      return;
    }
  }
  mrc_domain_get_idx3_patch_info(multi->domain, patch_idx_nei, info);
  //  get_rank_patch(ddc, patch_idx_nei, rank, patch);
}

// ----------------------------------------------------------------------
// ddc_init_outside

static void
ddc_init_outside(struct mrc_ddc *ddc, int p, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  struct mrc_patch_info info;
  mrc_ddc_multi_nei_get_patch_info(ddc, p, dir, &info);
  sr->nei_rank = info.rank;
  sr->nei_patch = info.patch;
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

  struct mrc_patch_info info;
  mrc_ddc_multi_nei_get_patch_info(ddc, p, dir, &info);
  sr->nei_rank = info.rank;
  sr->nei_patch = info.patch;
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
// mrc_ddc_multi_setup

static void
mrc_ddc_multi_setup(struct mrc_obj *obj)
{
  struct mrc_ddc *ddc = to_mrc_ddc(obj);
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
    mrc_domain_get_patch_idx3(multi->domain, p, ddc_patch->patch_idx);

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
}

// ----------------------------------------------------------------------
// ddc_run

static void
ddc_run(struct mrc_ddc *ddc, struct mrc_ddc_pattern *patt, int mb, int me,
	void *ctx,
	void (*to_buf)(int mb, int me, int ilo[3], int ihi[3], void *buf, void *ctx),
	void (*from_buf)(int mb, int me, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);
  int dir[3];

  // post all receives
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = mrc_ddc_dir2idx(dir);
	int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	struct mrc_ddc_sendrecv *r = &patt->recv[dir1];
	if (r->len > 0) {
#if 0
	  printf("[%d] recv from %d [%d,%d] x [%d,%d] x [%d,%d] len %d\n", ddc->rank,
		 r->rank_nei,
		 r->ilo[0], r->ihi[0], r->ilo[1], r->ihi[1], r->ilo[2], r->ihi[2],
		 r->len);
#endif
	  MPI_Irecv(r->buf, r->len * (me - mb), ddc->mpi_type, r->nei_rank,
		    0x1000 + dir1neg, ddc->obj.comm, &multi->recv_reqs[dir1]);
	} else {
	  multi->recv_reqs[dir1] = MPI_REQUEST_NULL;
	}
      }
    }
  }

  // post all sends
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = mrc_ddc_dir2idx(dir);
	struct mrc_ddc_sendrecv *s = &patt->send[dir1];
	if (s->len > 0) {
	  to_buf(mb, me, s->ilo, s->ihi, s->buf, ctx);
#if 0
	  printf("[%d] send to %d [%d,%d] x [%d,%d] x [%d,%d] len %d\n", ddc->rank,
		 s->rank_nei,
		 s->ilo[0], s->ihi[0], s->ilo[1], s->ihi[1], s->ilo[2], s->ihi[2],
		 s->len);
#endif
	  MPI_Isend(s->buf, s->len * (me - mb), ddc->mpi_type, s->nei_rank,
		    0x1000 + dir1, ddc->obj.comm, &multi->send_reqs[dir1]);
	} else {
	  multi->send_reqs[dir1] = MPI_REQUEST_NULL;
	}
      }
    }
  }

  // finish receives
  MPI_Waitall(27, multi->recv_reqs, MPI_STATUSES_IGNORE);
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = mrc_ddc_dir2idx(dir);
	struct mrc_ddc_sendrecv *r = &patt->recv[dir1];
	if (r->len > 0) {
	  from_buf(mb, me, r->ilo, r->ihi, r->buf, ctx);
	}
      }
    }
  }

  // finish sends
  MPI_Waitall(27, multi->send_reqs, MPI_STATUSES_IGNORE);
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

static struct mrc_ddc_ops mrc_ddc_multi_ops = {
  .name                  = "multi",
  .size                  = sizeof(struct mrc_ddc_multi),
  .setup                 = mrc_ddc_multi_setup,
  .set_domain            = mrc_ddc_multi_set_domain,
  .fill_ghosts           = mrc_ddc_multi_fill_ghosts,
  .add_ghosts            = mrc_ddc_multi_add_ghosts,
  .get_rank_nei          = mrc_ddc_get_rank_nei,
};

void
libmrc_ddc_register_multi()
{
  libmrc_ddc_register(&mrc_ddc_multi_ops);
}

