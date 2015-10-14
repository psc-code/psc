
#include "mrc_ddc_private.h"

#include <mrc_params.h>

#include <stdlib.h>
#include <assert.h>

// This entire file may be obsolete, since "multi" will do the same and more
// However, "multi" will do aggregation, while "simple" does not, so there may be
// in point in comparing performance...

// TODO: support dynamically adapting "max_n_fields" to whatever is needed
// TODO: add a fill_ghosts_fld() method and have it handle different data types
// dynamically

// ======================================================================
// mrc_ddc_simple

struct mrc_ddc_pattern {
  struct mrc_ddc_sendrecv send[N_DIR];
  struct mrc_ddc_sendrecv recv[N_DIR];
};

struct mrc_ddc_simple {
  // parameters
  int n_proc[3]; // # procs in 3D grid
  int ilo[3], ihi[3]; // local domain (no ghosts)
  int bc[3]; // boundary condition

  int proc[3]; // this proc's position in the 3D proc grid
  struct mrc_ddc_pattern add_ghosts;
  struct mrc_ddc_pattern fill_ghosts;
  MPI_Request send_reqs[N_DIR];
  MPI_Request recv_reqs[N_DIR];
};

#define mrc_ddc_simple(ddc) mrc_to_subobj(ddc, struct mrc_ddc_simple)

// ----------------------------------------------------------------------
// mrc_ddc_simple_get_rank_nei

static int
get_rank(struct mrc_ddc *ddc, const int proc[3])
{
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);
  for (int d = 0; d < 3; d++) {
    assert(proc[d] >= 0 && proc[d] < sub->n_proc[d]);
  }
  return (proc[2] * sub->n_proc[1] + proc[1]) * sub->n_proc[0] + proc[0];
}

static void
mrc_ddc_simple_get_nei_rank_patch(struct mrc_ddc *ddc, int p, int dir[3],
				  int *nei_rank, int *nei_patch)
{
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);
  assert(p == 0);

  int proc_nei[3];
  for (int d = 0; d < 3; d++) {
    proc_nei[d] = sub->proc[d] + dir[d];
    if (sub->bc[d] == BC_PERIODIC) {
      if (proc_nei[d] < 0) {
	proc_nei[d] += sub->n_proc[d];
      }
      if (proc_nei[d] >= sub->n_proc[d]) {
	proc_nei[d] -= sub->n_proc[d];
      }
    }
    if (proc_nei[d] < 0 || proc_nei[d] >= sub->n_proc[d]) {
      *nei_rank = -1;
      *nei_patch = -1;
      return;
    }
  }
  *nei_rank = get_rank(ddc, proc_nei);
  *nei_patch = 0;
}

// ----------------------------------------------------------------------
// ddc_init_outside

static void
ddc_init_outside(struct mrc_ddc *ddc, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);

  int dummy;
  mrc_ddc_simple_get_nei_rank_patch(ddc, 0, dir, &sr->nei_rank, &dummy);
  if (sr->nei_rank < 0)
    return;

  sr->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = sub->ilo[d] - ddc->ibn[d];
      sr->ihi[d] = sub->ilo[d];
      break;
    case 0:
      sr->ilo[d] = sub->ilo[d];
      sr->ihi[d] = sub->ihi[d];
      break;
    case 1:
      sr->ilo[d] = sub->ihi[d];
      sr->ihi[d] = sub->ihi[d] + ddc->ibn[d];
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
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);

  int dummy;
  mrc_ddc_simple_get_nei_rank_patch(ddc, 0, dir, &sr->nei_rank, &dummy);
  if (sr->nei_rank < 0)
    return;

  sr->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = sub->ilo[d];
      sr->ihi[d] = sub->ilo[d] + ddc->ibn[d];
      break;
    case 0:
      sr->ilo[d] = sub->ilo[d];
      sr->ihi[d] = sub->ihi[d];
      break;
    case 1:
      sr->ilo[d] = sub->ihi[d] - ddc->ibn[d];
      sr->ihi[d] = sub->ihi[d];
      break;
    }
    sr->len *= (sr->ihi[d] - sr->ilo[d]);
  }
  sr->buf = malloc(sr->len * ddc->max_n_fields * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// mrc_ddc_simple_setup

static void
mrc_ddc_simple_setup(struct mrc_ddc *ddc)
{
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);

  if (ddc->size_of_type == sizeof(float)) {
    ddc->mpi_type = MPI_FLOAT;
  } else if (ddc->size_of_type == sizeof(double)) {
    ddc->mpi_type = MPI_DOUBLE;
  } else {
    assert(0);
  }

  assert(sub->n_proc[0] * sub->n_proc[1] * sub->n_proc[2] == ddc->size);
  assert(ddc->max_n_fields > 0);

  int rr = ddc->rank;
  sub->proc[0] = rr % sub->n_proc[0]; rr /= sub->n_proc[0];
  sub->proc[1] = rr % sub->n_proc[1]; rr /= sub->n_proc[1];
  sub->proc[2] = rr;

  int dir[3];

  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
	  continue;

	ddc_init_outside(ddc, &sub->add_ghosts.send[mrc_ddc_dir2idx(dir)], dir);
	ddc_init_inside(ddc, &sub->add_ghosts.recv[mrc_ddc_dir2idx(dir)], dir);

	ddc_init_inside(ddc, &sub->fill_ghosts.send[mrc_ddc_dir2idx(dir)], dir);
	ddc_init_outside(ddc, &sub->fill_ghosts.recv[mrc_ddc_dir2idx(dir)], dir);
      }
    }
  }
}

// ----------------------------------------------------------------------
// ddc_run

static void
ddc_run(struct mrc_ddc *ddc, struct mrc_ddc_pattern *patt, int mb, int me,
	void *ctx,
	void (*to_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx),
	void (*from_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);
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
		    0x1000 + dir1neg, ddc->obj.comm, &sub->recv_reqs[dir1]);
	} else {
	  sub->recv_reqs[dir1] = MPI_REQUEST_NULL;
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
	  to_buf(mb, me, 0, s->ilo, s->ihi, s->buf, ctx);
#if 0
	  printf("[%d] send to %d [%d,%d] x [%d,%d] x [%d,%d] len %d\n", ddc->rank,
		 s->rank_nei,
		 s->ilo[0], s->ihi[0], s->ilo[1], s->ihi[1], s->ilo[2], s->ihi[2],
		 s->len);
#endif
	  MPI_Isend(s->buf, s->len * (me - mb), ddc->mpi_type, s->nei_rank,
		    0x1000 + dir1, ddc->obj.comm, &sub->send_reqs[dir1]);
	} else {
	  sub->send_reqs[dir1] = MPI_REQUEST_NULL;
	}
      }
    }
  }

  // finish receives
  MPI_Waitall(27, sub->recv_reqs, MPI_STATUSES_IGNORE);
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = mrc_ddc_dir2idx(dir);
	struct mrc_ddc_sendrecv *r = &patt->recv[dir1];
	if (r->len > 0) {
	  from_buf(mb, me, 0, r->ilo, r->ihi, r->buf, ctx);
	}
      }
    }
  }

  // finish sends
  MPI_Waitall(27, sub->send_reqs, MPI_STATUSES_IGNORE);
}

// ----------------------------------------------------------------------
// mrc_ddc_simple_add_ghosts

static void
mrc_ddc_simple_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);

  assert(me - mb <= ddc->max_n_fields);
  ddc_run(ddc, &sub->add_ghosts, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->add_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_simple_fill_ghosts

static void
mrc_ddc_simple_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_simple *sub = mrc_ddc_simple(ddc);

  assert(me - mb <= ddc->max_n_fields);
  ddc_run(ddc, &sub->fill_ghosts, mb, me, ctx,
	  ddc->funcs->copy_to_buf, ddc->funcs->copy_from_buf);
}

// ======================================================================
// mrc_ddc_simple_ops

#define VAR(x) (void *)offsetof(struct mrc_ddc_simple, x)
static struct param mrc_ddc_simple_params_descr[] = {
  { "n_proc"          , VAR(n_proc)       , PARAM_INT3(0, 0, 0)    },
  { "ilo"             , VAR(ilo)          , PARAM_INT3(0, 0, 0)    },
  { "ihi"             , VAR(ihi)          , PARAM_INT3(0, 0, 0)    },
  { "bc"              , VAR(bc)           , PARAM_INT3(0, 0, 0)    }, // FIXME, select
  {}
};

struct mrc_ddc_ops mrc_ddc_simple_ops = {
  .name                  = "simple",
  .size                  = sizeof(struct mrc_ddc_simple),
  .param_descr           = mrc_ddc_simple_params_descr,
  .setup                 = mrc_ddc_simple_setup,
  .fill_ghosts           = mrc_ddc_simple_fill_ghosts,
  .add_ghosts            = mrc_ddc_simple_add_ghosts,
};

// ======================================================================
// mrc_ddc_funcs_fld for mrc_fld

#include <mrc_fld.h>

static void
mrc_fld_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  //  mprintf("to buf %d:%d x %d:%d x %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
  struct mrc_fld *fld = ctx;
  float *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf,m - mb, ix,iy,iz) = MRC_F3(fld,m, ix,iy,iz);
	}
      }
    }
  }
}

static void
mrc_fld_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  //  mprintf("from buf %d:%d x %d:%d x %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
  struct mrc_fld *fld = ctx;
  float *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_F3(fld,m, ix,iy,iz) = MRC_DDC_BUF3(buf,m - mb, ix,iy,iz);
	}
      }
    }
  }
}

struct mrc_ddc_funcs mrc_ddc_funcs_fld = {
  .copy_to_buf   = mrc_fld_copy_to_buf,
  .copy_from_buf = mrc_fld_copy_from_buf,
};

