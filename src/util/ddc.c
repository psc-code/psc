
#include "ddc.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

static int
get_rank(struct ddc_subdomain *ddc, const int proc[3])
{
  for (int d = 0; d < 3; d++) {
    assert(proc[d] >= 0 && proc[d] < ddc->prm.n_proc[d]);
  }
  return (proc[2] * ddc->prm.n_proc[1] + proc[1]) * ddc->prm.n_proc[0] + proc[0];
}

int
ddc_get_rank_nei(struct ddc_subdomain *ddc, int dir[3])
{
  int proc_nei[3];
  for (int d = 0; d < 3; d++) {
    proc_nei[d] = ddc->proc[d] + dir[d];
    if (ddc->prm.bc[d] == DDC_BC_PERIODIC) {
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
ddc_init_outside(struct ddc_subdomain *ddc, struct ddc_sendrecv *sr, int dir[3])
{
  sr->rank_nei = ddc_get_rank_nei(ddc, dir);
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
  sr->buf = calloc(sr->len, sizeof(*sr->buf));
}

// ----------------------------------------------------------------------
// ddc_init_inside

static void
ddc_init_inside(struct ddc_subdomain *ddc, struct ddc_sendrecv *sr, int dir[3])
{
  sr->rank_nei = ddc_get_rank_nei(ddc, dir);
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
  sr->buf = calloc(sr->len, sizeof(*sr->buf));
}

// ----------------------------------------------------------------------
// ddc_run

static void
ddc_run(struct ddc_subdomain *ddc, struct ddc_pattern *patt, int m,
	void (*to_buf)(int m, int ilo[3], int ihi[3], ddc_real *buf),
	void (*from_buf)(int m, int ilo[3], int ihi[3], ddc_real *buf))
{
  int dir[3];

  // post all receives
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = dir2idx(dir);
	int dir1neg = dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	struct ddc_sendrecv *r = &patt->recv[dir1];
	if (r->len > 0) {
#if 0
	  printf("[%d] recv from %d [%d,%d] x [%d,%d] x [%d,%d] len %d\n", ddc->rank,
		 r->rank_nei,
		 r->ilo[0], r->ihi[0], r->ilo[1], r->ihi[1], r->ilo[2], r->ihi[2],
		 r->len);
#endif
	  MPI_Irecv(r->buf, r->len, MPI_DOUBLE, r->rank_nei, 0x1000 + dir1neg,
		    ddc->prm.comm, &ddc->recv_reqs[dir1]);
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
	int dir1 = dir2idx(dir);
	struct ddc_sendrecv *s = &patt->send[dir1];
	if (s->len > 0) {
	  to_buf(m, s->ilo, s->ihi, s->buf);
#if 0
	  printf("[%d] send to %d [%d,%d] x [%d,%d] x [%d,%d] len %d\n", ddc->rank,
		 s->rank_nei,
		 s->ilo[0], s->ihi[0], s->ilo[1], s->ihi[1], s->ilo[2], s->ihi[2],
		 s->len);
#endif
	  MPI_Isend(s->buf, s->len, MPI_DOUBLE, s->rank_nei, 0x1000 + dir1,
		    ddc->prm.comm, &ddc->send_reqs[dir1]);
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
	int dir1 = dir2idx(dir);
	struct ddc_sendrecv *r = &patt->recv[dir1];
	if (r->len > 0) {
	  from_buf(m, r->ilo, r->ihi, r->buf);
	}
      }
    }
  }

  // finish sends
  MPI_Waitall(27, ddc->send_reqs, MPI_STATUSES_IGNORE);
}

// ----------------------------------------------------------------------
// ddc_add_ghosts

void
ddc_add_ghosts(struct ddc_subdomain *ddc, int m)
{
  ddc_run(ddc, &ddc->add_ghosts, m, ddc->prm.copy_to_buf, ddc->prm.add_from_buf);
}

// ----------------------------------------------------------------------
// ddc_fill_ghosts

void
ddc_fill_ghosts(struct ddc_subdomain *ddc, int m)
{
  ddc_run(ddc, &ddc->fill_ghosts, m, ddc->prm.copy_to_buf, ddc->prm.copy_from_buf);
}

// ----------------------------------------------------------------------
// ddc_create

struct ddc_subdomain *
ddc_create(struct ddc_params *prm)
{
  struct ddc_subdomain *ddc = malloc(sizeof(*ddc));
  memset(ddc, 0, sizeof(*ddc));
  
  ddc->prm = *prm;
  MPI_Comm_rank(prm->comm, &ddc->rank);
  MPI_Comm_size(prm->comm, &ddc->size);

  assert(prm->n_proc[0] * prm->n_proc[1] * prm->n_proc[2] == ddc->size);

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

	ddc_init_outside(ddc, &ddc->add_ghosts.send[dir2idx(dir)], dir);
	ddc_init_inside(ddc, &ddc->add_ghosts.recv[dir2idx(dir)], dir);

	ddc_init_inside(ddc, &ddc->fill_ghosts.send[dir2idx(dir)], dir);
	ddc_init_outside(ddc, &ddc->fill_ghosts.recv[dir2idx(dir)], dir);
      }
    }
  }

  return ddc;
}

