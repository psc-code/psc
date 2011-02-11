
#include "mrc_ddc.h"
#include "mrc_ddc_private.h"

#ifdef PROF_DDC
#include "mrc_profile.h"
#else
#include "mrc_profile_nop.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define ID_FILLGHOSTS (0x1000)

static int PR_ddc_fill_ghosts_begin;
static int PR_ddc_fill_ghosts_end;
static int PR_ddc_copy_to_buf;
static int PR_ddc_copy_from_buf;

static inline int
to_rank(struct mrc_ddc *ddc, const int ip[3])
{
  return ((ip[2] * ddc->params.np[1]) + ip[1]) * ddc->params.np[0] + ip[0];
}

struct mrc_ddc *
mrc_ddc_create(MPI_Comm comm, struct mrc_ddc_params *params, struct mrc_ddc_ops *ops)
{
  struct mrc_ddc *ddc = malloc(sizeof(*ddc));
  memset(ddc, 0, sizeof(*ddc));
  
  int *n = params->n;
  int *np = params->np;
  int bnd = params->bnd;

  ddc->ops = ops;
  ddc->comm = comm;
  ddc->params = *params;

  MPI_Comm_rank(comm, &ddc->rank);
  MPI_Comm_size(comm, &ddc->size);
  assert(ddc->size == np[0] * np[1] * np[2]);

  int m = ddc->rank;
  ddc->ip[0] = m % np[0]; m = m / np[0];
  ddc->ip[1] = m % np[1]; m = m / np[1];
  ddc->ip[2] = m;
  //  printf("[%d] ip = %d,%d,%d\n", ddc->rank, ddc->ipx, ddc->ipy, ddc->ipz);
  assert(to_rank(ddc, ddc->ip) == ddc->rank);

  for (int dir3 = 0; dir3 < N_DIR; dir3++) {
    struct mrc_ddc_dir *dir = &ddc->dir[dir3];
    int m = dir3;
    dir->dir[0] = m % 3 - 1; m /= 3;
    dir->dir[1] = m % 3 - 1; m /= 3;
    dir->dir[2] = m     - 1;
    
    dir->recv_rank = -1;
    dir->send_rank = -1;
    if (dir->dir[0] == 0 && dir->dir[1] == 0 && dir->dir[2] == 0)
      continue;
    
    for (int d = 0; d < 3; d++) {
      switch (dir->dir[d]) {
      case -1: dir->recv_ib[d] = n[d]; dir->recv_ie[d] = n[d] + bnd; break;
      case  0: dir->recv_ib[d] = 0;    dir->recv_ie[d] = n[d];       break;
      case  1: dir->recv_ib[d] = -bnd; dir->recv_ie[d] = 0;          break;  
      }

      switch (dir->dir[d]) {
      case -1: dir->send_ib[d] = 0;          dir->send_ie[d] = bnd;  break;
      case  0: dir->send_ib[d] = 0;          dir->send_ie[d] = n[d]; break;
      case  1: dir->send_ib[d] = n[d] - bnd; dir->send_ie[d] = n[d]; break;  
      }

      dir->in[d] = dir->recv_ie[d] - dir->recv_ib[d];
    }
    
    dir->n_bnd = dir->in[0] * dir->in[1] * dir->in[2];
    
    int ipn_recv[3], ipn_send[3], is_in_domain_recv = 1, is_in_domain_send = 1;
    for (int d = 0; d < 3; d++) {
      ipn_recv[d] = ddc->ip[d] - dir->dir[d];
      if (ddc->params.bc[d] == BC_PERIODIC) {
	ipn_recv[d] = (ipn_recv[d] + np[d]) % np[d];
      } else {
	if (ipn_recv[d] < 0 || ipn_recv[d] >= np[d]) {
	  is_in_domain_recv = 0;
	}
      }
      ipn_send[d] = ddc->ip[d] + dir->dir[d];
      if (ddc->params.bc[d] == BC_PERIODIC) {
	ipn_send[d] = (ipn_send[d] + np[d]) % np[d];
      } else {
	if (ipn_send[d] < 0 || ipn_send[d] >= np[d]) {
	  is_in_domain_send = 0;
	}
      }
    }
    if (is_in_domain_recv) {
      dir->recv_rank = to_rank(ddc, ipn_recv);
      dir->recv_buf = calloc(params->max_n_comp * dir->n_bnd, sizeof(float));
    }
    if (is_in_domain_send) {
      dir->send_rank = to_rank(ddc, ipn_send);
      dir->send_buf = calloc(params->max_n_comp * dir->n_bnd, sizeof(float));
    }
  }

  if (!PR_ddc_fill_ghosts_begin) {
    PR_ddc_fill_ghosts_begin = prof_register("ddc_fill_ghosts_begin", 1., 0, 0);
    PR_ddc_fill_ghosts_end = prof_register("ddc_fill_ghosts_end", 1., 0, 0);
    PR_ddc_copy_to_buf = prof_register("ddc_copy_to_buf", 1., 0, 0);
    PR_ddc_copy_from_buf = prof_register("ddc_copy_from_buf", 1., 0, 0);
  }
  return ddc;
}

void
mrc_ddc_destroy(struct mrc_ddc *ddc)
{
  for (int dir3 = 0; dir3 < N_DIR; dir3++) {
    struct mrc_ddc_dir *dir = &ddc->dir[dir3];

    free(dir->recv_buf);
    free(dir->send_buf);
  }
  free(ddc);
}

void
mrc_ddc_fill_ghosts_begin(struct mrc_ddc *ddc, float *x, int mb, int me)
{
  prof_start(PR_ddc_fill_ghosts_begin);
  int mn = me - mb; // # components
  assert(mn <= ddc->params.max_n_comp);
  // post all receives
  ddc->n_recv_reqs = 0;
  for (int dir3 = 0; dir3 < N_DIR; dir3++) {
    struct mrc_ddc_dir *dir = &ddc->dir[dir3];

    if (dir->recv_rank >= 0) {
      //      printf("[%d] recv from %d\n", ddc->rank, dir->recv_rank);
      MPI_Irecv(dir->recv_buf, mn*dir->n_bnd, MPI_FLOAT, dir->recv_rank,
		ID_FILLGHOSTS+dir3, ddc->comm, &ddc->recv_reqs[ddc->n_recv_reqs++]);
    }
  }

  // send out all boundary points
  ddc->n_send_reqs = 0;
  for (int dir3 = 0; dir3 < N_DIR; dir3++) {
    struct mrc_ddc_dir *dir = &ddc->dir[dir3];

    if (dir->send_rank >= 0) {
      prof_start(PR_ddc_copy_to_buf);
      ddc->ops->copy_to_buf(&ddc->params, x, dir->send_buf, dir->send_ib, dir->send_ie, mb, me);
      //      printf("[%d] send to %d\n", ddc->rank, dir->send_rank);
      prof_stop(PR_ddc_copy_to_buf);
      MPI_Isend(dir->send_buf, mn*dir->n_bnd, MPI_FLOAT, dir->send_rank,
		ID_FILLGHOSTS+dir3, ddc->comm, &ddc->send_reqs[ddc->n_send_reqs++]);
    }
  }
  prof_stop(PR_ddc_fill_ghosts_begin);
}

void
mrc_ddc_fill_ghosts_end(struct mrc_ddc *ddc, float *x, int mb, int me)
{
  prof_start(PR_ddc_fill_ghosts_end);
  // wait for all receives -- FIXME, could do copying once any is done
  MPI_Waitall(ddc->n_recv_reqs, ddc->recv_reqs, MPI_STATUSES_IGNORE);

  // after receive
  for (int dir3 = 0; dir3 < N_DIR; dir3++) {
    struct mrc_ddc_dir *dir = &ddc->dir[dir3];

    if (dir->recv_rank >= 0) {
      prof_start(PR_ddc_copy_from_buf);
      ddc->ops->copy_from_buf(&ddc->params, x, dir->recv_buf, dir->recv_ib, dir->recv_ie, mb, me);
      prof_stop(PR_ddc_copy_from_buf);
    }
  }

  MPI_Waitall(ddc->n_send_reqs, ddc->send_reqs, MPI_STATUSES_IGNORE);
  prof_stop(PR_ddc_fill_ghosts_end);
}

void
mrc_ddc_fill_ghosts(struct mrc_ddc *ddc, float *x, int mb, int me)
{
  mrc_ddc_fill_ghosts_begin(ddc, x, mb, me);
  mrc_ddc_fill_ghosts_end(ddc, x, mb, me);
}

// ======================================================================

#define GFLD3(f,m, ix,iy,iz) ((f)[((((m)*nnz+(iz))*nny)+(iy))*nnx+(ix)])

static void
fld3_copy_to_buf(struct mrc_ddc_params *ddc_params, float *x, float *buf,
		 const int ib[3], const int ie[3], int mb, int me)
{
  int bnd = ddc_params->bnd;
  int nnx = ddc_params->n[0] + 2*bnd;
  int nny = ddc_params->n[1] + 2*bnd;
  int nnz = ddc_params->n[2] + 2*bnd;

  for (int iz = ib[2]; iz < ie[2]; iz++) {
    for (int iy = ib[1]; iy < ie[1]; iy++) {
      for (int ix = ib[0]; ix < ie[0]; ix++) {
	for (int m = mb; m < me; m++) {
	  MRC_DDC_BUF3(buf,m, ix,iy,iz) = GFLD3(x,m, ix+bnd,iy+bnd,iz+bnd);
	}
      }
    }
  }
}

static void
fld3_copy_from_buf(struct mrc_ddc_params *ddc_params, float *x, float *buf,
		   const int ib[3], const int ie[3], int mb, int me)
{
  int bnd = ddc_params->bnd;
  int nnx = ddc_params->n[0] + 2*bnd;
  int nny = ddc_params->n[1] + 2*bnd;
  int nnz = ddc_params->n[2] + 2*bnd;

  for (int iz = ib[2]; iz < ie[2]; iz++) {
    for (int iy = ib[1]; iy < ie[1]; iy++) {
      for (int ix = ib[0]; ix < ie[0]; ix++) {
	for (int m = mb; m < me; m++) {
	  GFLD3(x,m, ix+bnd,iy+bnd,iz+bnd) = MRC_DDC_BUF3(buf,m, ix,iy,iz);
	}
      }
    }
  }
}

// for standard Fortran layout fields (component: slow idx)

struct mrc_ddc_ops mrc_ddc_ops_fortran = {
  .copy_to_buf   = fld3_copy_to_buf,
  .copy_from_buf = fld3_copy_from_buf,
};

