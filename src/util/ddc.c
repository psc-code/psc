
#include "ddc.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

static inline int
dir2idx(int dir[3])
{
  return ((dir[2] + 1) * 3 + dir[1] + 1) * 3 + dir[0] + 1;
}

static int
get_rank(struct ddc_subdomain *ddc, const int proc[3])
{
  for (int d = 0; d < 3; d++) {
    assert(proc[d] >= 0 && proc[d] < ddc->prm.n_proc[d]);
  }
  return (proc[2] * ddc->prm.n_proc[1] + proc[1]) * ddc->prm.n_proc[0] + proc[0];
}

static int
get_rank_nei(struct ddc_subdomain *ddc, int dir[3])
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
// ddc_send_init

static void
ddc_send_init(struct ddc_subdomain *ddc, int dir[3])
{
  struct ddc_send *s = &ddc->send[dir2idx(dir)];

  s->rank_nei = get_rank_nei(ddc, dir);
  if (s->rank_nei < 0)
    return;

  s->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      s->ilo[d] = ddc->prm.ilo[d] - ddc->prm.ibn[d];
      s->ihi[d] = ddc->prm.ilo[d];
      break;
    case 0:
      s->ilo[d] = ddc->prm.ilo[d];
      s->ihi[d] = ddc->prm.ihi[d];
      break;
    case 1:
      s->ilo[d] = ddc->prm.ihi[d];
      s->ihi[d] = ddc->prm.ihi[d] + ddc->prm.ibn[d];
      break;
    }
    s->len *= (s->ihi[d] - s->ilo[d]);
  }
  s->buf = calloc(s->len, sizeof(*s->buf));
}

// ----------------------------------------------------------------------
// ddc_recv_init

static void
ddc_recv_init(struct ddc_subdomain *ddc, int dir[3])
{
  struct ddc_recv *r = &ddc->recv[dir2idx(dir)];

  r->rank_nei = get_rank_nei(ddc, dir);
  if (r->rank_nei < 0)
    return;

  r->len = 1;
  for (int d = 0; d < 3; d++) {
    switch (dir[d]) {
    case -1:
      r->ilo[d] = ddc->prm.ilo[d];
      r->ihi[d] = ddc->prm.ilo[d] + ddc->prm.ibn[d];
      break;
    case 0:
      r->ilo[d] = ddc->prm.ilo[d];
      r->ihi[d] = ddc->prm.ihi[d];
      break;
    case 1:
      r->ilo[d] = ddc->prm.ihi[d] - ddc->prm.ibn[d];
      r->ihi[d] = ddc->prm.ihi[d];
      break;
    }
    r->len *= (r->ihi[d] - r->ilo[d]);
  }
  r->buf = calloc(r->len, sizeof(*r->buf));
}

// ----------------------------------------------------------------------
// ddc_add_ghosts

void
ddc_add_ghosts(struct ddc_subdomain *ddc, int m)
{
  int dir[3];

  // post all receives
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = dir2idx(dir);
	struct ddc_recv *r = &ddc->recv[dir1];
	if (r->len > 0) {
#if 0
	  printf("[%d] recv from %d [%d,%d] x [%d,%d] x [%d,%d] len %d\n", ddc->rank,
		 r->rank_nei,
		 r->ilo[0], r->ihi[0], r->ilo[1], r->ihi[1], r->ilo[2], r->ihi[2],
		 r->len);
#endif
	  MPI_Irecv(r->buf, r->len, MPI_DOUBLE, r->rank_nei, 0x1000, ddc->prm.comm,
		   &ddc->recv_reqs[dir1]);
	  ddc->prm.add_from_buf(m, r->ilo, r->ihi, r->buf);
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
	struct ddc_send *s = &ddc->send[dir1];
	if (s->len > 0) {
	  ddc->prm.copy_to_buf(m, s->ilo, s->ihi, s->buf);
#if 0
	  printf("[%d] send to %d [%d,%d] x [%d,%d] x [%d,%d] len %d\n", ddc->rank,
		 s->rank_nei,
		 s->ilo[0], s->ihi[0], s->ilo[1], s->ihi[1], s->ilo[2], s->ihi[2],
		 s->len);
#endif
	  MPI_Isend(s->buf, s->len, MPI_DOUBLE, s->rank_nei, 0x1000, ddc->prm.comm,
		    &ddc->send_reqs[dir1]);
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
	struct ddc_recv *r = &ddc->recv[dir1];
	if (r->len > 0) {
	  ddc->prm.add_from_buf(m, r->ilo, r->ihi, r->buf);
	}
      }
    }
  }

  // finish sends
  MPI_Waitall(27, ddc->send_reqs, MPI_STATUSES_IGNORE);
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

  // send to right
  ddc_send_init(ddc, (int[3]) { 1, 0, 0 });

  // recv from left
  ddc_recv_init(ddc, (int[3]) { -1, 0, 0 });

  // send to right
  ddc_send_init(ddc, (int[3]) { -1, 0, 0 });

  // recv from right
  ddc_recv_init(ddc, (int[3]) { 1, 0, 0 });

  return ddc;
}

