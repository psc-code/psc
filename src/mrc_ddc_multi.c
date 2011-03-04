
#include "mrc_ddc_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>

#include <stdlib.h>
#include <assert.h>

#define to_mrc_ddc_multi(ddc) ((struct mrc_ddc_multi *) (ddc)->obj.subctx)

// ----------------------------------------------------------------------
// mrc_ddc_multi_get_rank_nei

static int
get_rank(struct mrc_ddc *ddc, const int patch_idx[3])
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  for (int d = 0; d < 3; d++) {
    assert(patch_idx[d] >= 0 && patch_idx[d] < multi->np[d]);
  }
  return (patch_idx[2] * multi->np[1] + patch_idx[1]) * multi->np[0] + patch_idx[0];
}

static int
mrc_ddc_multi_get_rank_nei(struct mrc_ddc *ddc, int dir[3])
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  int proc_nei[3];
  for (int d = 0; d < 3; d++) {
    proc_nei[d] = ddc->proc[d] + dir[d];
    if (multi->bc[d] == BC_PERIODIC) {
      if (proc_nei[d] < 0) {
	proc_nei[d] += multi->np[d];
      }
      if (proc_nei[d] >= multi->np[d]) {
	proc_nei[d] -= multi->np[d];
      }
    }
    if (proc_nei[d] < 0 || proc_nei[d] >= multi->np[d])
      return -1;
  }
  return get_rank(ddc, proc_nei);
}

// ----------------------------------------------------------------------
// ddc_init_outside

static void
ddc_init_outside(struct mrc_ddc *ddc, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  sr->rank_nei = mrc_ddc_multi_get_rank_nei(ddc, dir);
  if (sr->rank_nei < 0)
    return;

  sr->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = multi->ilo[d] - ddc->ibn[d];
      sr->ihi[d] = multi->ilo[d];
      break;
    case 0:
      sr->ilo[d] = multi->ilo[d];
      sr->ihi[d] = multi->ihi[d];
      break;
    case 1:
      sr->ilo[d] = multi->ihi[d];
      sr->ihi[d] = multi->ihi[d] + ddc->ibn[d];
      break;
    }
    sr->len *= (sr->ihi[d] - sr->ilo[d]);
  }
  sr->buf = malloc(sr->len * ddc->max_n_fields * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// ddc_init_inside

static void
ddc_init_inside(struct mrc_ddc *ddc, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  struct mrc_ddc_multi *multi = to_mrc_ddc_multi(ddc);

  sr->rank_nei = mrc_ddc_multi_get_rank_nei(ddc, dir);
  if (sr->rank_nei < 0)
    return;

  sr->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = multi->ilo[d];
      sr->ihi[d] = multi->ilo[d] + ddc->ibn[d];
      break;
    case 0:
      sr->ilo[d] = multi->ilo[d];
      sr->ihi[d] = multi->ihi[d];
      break;
    case 1:
      sr->ilo[d] = multi->ihi[d] - ddc->ibn[d];
      sr->ihi[d] = multi->ihi[d];
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

  int rr = ddc->rank;
  ddc->proc[0] = rr % multi->np[0]; rr /= multi->np[0];
  ddc->proc[1] = rr % multi->np[1]; rr /= multi->np[1];
  ddc->proc[2] = rr;

  int dir[3];

  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
	  continue;

	ddc_init_outside(ddc, &ddc->add_ghosts.send[mrc_ddc_dir2idx(dir)], dir);
	ddc_init_inside(ddc, &ddc->add_ghosts.recv[mrc_ddc_dir2idx(dir)], dir);

	ddc_init_inside(ddc, &ddc->fill_ghosts.send[mrc_ddc_dir2idx(dir)], dir);
	ddc_init_outside(ddc, &ddc->fill_ghosts.recv[mrc_ddc_dir2idx(dir)], dir);
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
	  MPI_Irecv(r->buf, r->len * (me - mb), ddc->mpi_type, r->rank_nei,
		    0x1000 + dir1neg, ddc->obj.comm, &ddc->recv_reqs[dir1]);
	} else {
	  ddc->recv_reqs[dir1] = MPI_REQUEST_NULL;
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
	  MPI_Isend(s->buf, s->len * (me - mb), ddc->mpi_type, s->rank_nei,
		    0x1000 + dir1, ddc->obj.comm, &ddc->send_reqs[dir1]);
	} else {
	  ddc->send_reqs[dir1] = MPI_REQUEST_NULL;
	}
      }
    }
  }

  // finish receives
  MPI_Waitall(27, ddc->recv_reqs, MPI_STATUSES_IGNORE);
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
  MPI_Waitall(27, ddc->send_reqs, MPI_STATUSES_IGNORE);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_add_ghosts

static void
mrc_ddc_multi_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  ddc_run(ddc, &ddc->add_ghosts, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->add_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_fill_ghosts

static void
mrc_ddc_multi_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  ddc_run(ddc, &ddc->fill_ghosts, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->copy_from_buf);
}

// ======================================================================
// mrc_ddc_multi_ops

#define VAR(x) (void *)offsetof(struct mrc_ddc_multi, x)
static struct param mrc_ddc_multi_params_descr[] = {
  { "ilo"             , VAR(ilo)          , PARAM_INT3(0, 0, 0)    },
  { "ihi"             , VAR(ihi)          , PARAM_INT3(0, 0, 0)    },
  { "bc"              , VAR(bc)           , PARAM_INT3(0, 0, 0)    }, // FIXME, select
  {}
};

static struct mrc_ddc_ops mrc_ddc_multi_ops = {
  .name                  = "multi",
  .size                  = sizeof(struct mrc_ddc_multi),
  .param_descr           = mrc_ddc_multi_params_descr,
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

