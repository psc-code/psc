
#include "mrc_ddc.h"
#include "mrc_ddc_private.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

static int
get_rank(struct mrc_ddc *ddc, const int proc[3])
{
  for (int d = 0; d < 3; d++) {
    assert(proc[d] >= 0 && proc[d] < ddc->prm.n_proc[d]);
  }
  return (proc[2] * ddc->prm.n_proc[1] + proc[1]) * ddc->prm.n_proc[0] + proc[0];
}

int
mrc_ddc_get_rank_nei(struct mrc_ddc *ddc, int dir[3])
{
  int proc_nei[3];
  for (int d = 0; d < 3; d++) {
    proc_nei[d] = ddc->proc[d] + dir[d];
    if (ddc->prm.bc[d] == BC_PERIODIC) {
      if (proc_nei[d] < 0) {
	proc_nei[d] += ddc->prm.n_proc[d];
      }
      if (proc_nei[d] >= ddc->prm.n_proc[d]) {
	proc_nei[d] -= ddc->prm.n_proc[d];
      }
    }
    if (proc_nei[d] < 0 || proc_nei[d] >= ddc->prm.n_proc[d])
      return - 1;
  }
  return get_rank(ddc, proc_nei);
}

// ----------------------------------------------------------------------
// ddc_init_outside

static void
ddc_init_outside(struct mrc_ddc *ddc, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  sr->rank_nei = mrc_ddc_get_rank_nei(ddc, dir);
  if (sr->rank_nei < 0)
    return;

  sr->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = ddc->prm.ilo[d] - ddc->prm.ibn[d];
      sr->ihi[d] = ddc->prm.ilo[d];
      break;
    case 0:
      sr->ilo[d] = ddc->prm.ilo[d];
      sr->ihi[d] = ddc->prm.ihi[d];
      break;
    case 1:
      sr->ilo[d] = ddc->prm.ihi[d];
      sr->ihi[d] = ddc->prm.ihi[d] + ddc->prm.ibn[d];
      break;
    }
    sr->len *= (sr->ihi[d] - sr->ilo[d]);
  }
  sr->buf = malloc(sr->len * ddc->prm.max_n_fields * ddc->prm.size_of_type);
}

// ----------------------------------------------------------------------
// ddc_init_inside

static void
ddc_init_inside(struct mrc_ddc *ddc, struct mrc_ddc_sendrecv *sr, int dir[3])
{
  sr->rank_nei = mrc_ddc_get_rank_nei(ddc, dir);
  if (sr->rank_nei < 0)
    return;

  sr->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      sr->ilo[d] = ddc->prm.ilo[d];
      sr->ihi[d] = ddc->prm.ilo[d] + ddc->prm.ibn[d];
      break;
    case 0:
      sr->ilo[d] = ddc->prm.ilo[d];
      sr->ihi[d] = ddc->prm.ihi[d];
      break;
    case 1:
      sr->ilo[d] = ddc->prm.ihi[d] - ddc->prm.ibn[d];
      sr->ihi[d] = ddc->prm.ihi[d];
      break;
    }
    sr->len *= (sr->ihi[d] - sr->ilo[d]);
  }
  sr->buf = malloc(sr->len * ddc->prm.max_n_fields * ddc->prm.size_of_type);
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
		    0x1000 + dir1neg, ddc->comm, &ddc->recv_reqs[dir1]);
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
		    0x1000 + dir1, ddc->comm, &ddc->send_reqs[dir1]);
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
// mrc_ddc_add_ghosts

void
mrc_ddc_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  ddc_run(ddc, &ddc->add_ghosts, mb, me, ctx,
	  ddc->ops->copy_to_buf, ddc->ops->add_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_fill_ghosts

void
mrc_ddc_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  ddc_run(ddc, &ddc->fill_ghosts, mb, me, ctx,
	  ddc->ops->copy_to_buf, ddc->ops->copy_from_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_create

struct mrc_ddc *
mrc_ddc_create(MPI_Comm comm)
{
  struct mrc_ddc *ddc = malloc(sizeof(*ddc));
  memset(ddc, 0, sizeof(*ddc));

  MPI_Comm ncomm;
  MPI_Comm_dup(comm, &ncomm);
  ddc->comm = ncomm;
  MPI_Comm_rank(comm, &ddc->rank);
  MPI_Comm_size(comm, &ddc->size);

  return ddc;
}

// ----------------------------------------------------------------------
// mrc_ddc_set_ops

void
mrc_ddc_set_ops(struct mrc_ddc *ddc, struct mrc_ddc_ops *ops)
{
  ddc->ops = ops;
}

// ----------------------------------------------------------------------
// mrc_ddc_set_params

void
mrc_ddc_set_params(struct mrc_ddc *ddc, struct mrc_ddc_params *prm)
{
  ddc->prm = *prm;
}

// ----------------------------------------------------------------------
// mrc_ddc_setup

void
mrc_ddc_setup(struct mrc_ddc *ddc)
{
  if (ddc->prm.size_of_type == sizeof(float)) {
    ddc->mpi_type = MPI_FLOAT;
  } else if (ddc->prm.size_of_type == sizeof(double)) {
    ddc->mpi_type = MPI_DOUBLE;
  } else {
    assert(0);
  }

  assert(ddc->prm.n_proc[0] * ddc->prm.n_proc[1] * ddc->prm.n_proc[2] == ddc->size);
  assert(ddc->prm.max_n_fields > 0);

  int rr = ddc->rank;
  ddc->proc[0] = rr % ddc->prm.n_proc[0]; rr /= ddc->prm.n_proc[0];
  ddc->proc[1] = rr % ddc->prm.n_proc[1]; rr /= ddc->prm.n_proc[1];
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
// mrc_ddc_destroy

void
mrc_ddc_destroy(struct mrc_ddc *ddc)
{
  MPI_Comm_free(&ddc->comm);
  free(ddc);
}

// ======================================================================
// mrc_ddc_ops_f3 for mrc_f3

#include <mrc_fld.h>

// FIXME, 0-based offsets and ghost points don't match well (not pretty anyway)

static void
mrc_f3_copy_to_buf(int mb, int me, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct mrc_f3 *fld = ctx;
  float *buf = _buf;
  int bnd = fld->sw;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf,m - mb, ix,iy,iz) = MRC_F3(fld,m, ix+bnd,iy+bnd,iz+bnd);
	}
      }
    }
  }
}

static void
mrc_f3_copy_from_buf(int mb, int me, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct mrc_f3 *fld = ctx;
  float *buf = _buf;
  int bnd = fld->sw;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_F3(fld,m, ix+bnd,iy+bnd,iz+bnd) = MRC_DDC_BUF3(buf,m - mb, ix,iy,iz);
	}
      }
    }
  }
}

struct mrc_ddc_ops mrc_ddc_ops_f3 = {
  .copy_to_buf   = mrc_f3_copy_to_buf,
  .copy_from_buf = mrc_f3_copy_from_buf,
};

