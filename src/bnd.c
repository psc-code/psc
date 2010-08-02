
#include "psc.h"
#include "util/profile.h"
#include "util/ddc.h"

#include <mpi.h>
#include <string.h>

// ======================================================================
// C bnd

struct c_bnd_ctx {
  struct ddc_subdomain *ddc;
  struct ddc_particles *ddcp;
};

static void
copy_to_buf(int fld_nr, int ilo[3], int ihi[3], void *_buf)
{
  fields_base_real_t *buf = _buf;

  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	DDC_BUF(buf, ix,iy,iz) = F3_BASE(fld_nr, ix,iy,iz);
      }
    }
  }
}

static void
add_from_buf(int fld_nr, int ilo[3], int ihi[3], void *_buf)
{
  fields_base_real_t *buf = _buf;
  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	F3_BASE(fld_nr, ix,iy,iz) += DDC_BUF(buf, ix,iy,iz);
      }
    }
  }
}

static void
copy_from_buf(int fld_nr, int ilo[3], int ihi[3], void *_buf)
{
  fields_base_real_t *buf = _buf;
  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	F3_BASE(fld_nr, ix,iy,iz) = DDC_BUF(buf, ix,iy,iz);
      }
    }
  }
}

// ======================================================================

struct ddcp_nei {
  struct f_particle *send_buf;
  int n_send;
  int n_recv;
  int rank;
  int send_buf_size;
};

struct ddc_particles {
  struct ddcp_nei nei[27];
  MPI_Request send_reqs[27];
  MPI_Request sendp_reqs[27];
  MPI_Request recv_reqs[27];
  int head;
};

static struct ddc_particles *
ddc_particles_create(struct ddc_subdomain *ddc)
{
  struct ddc_particles *ddcp = malloc(sizeof(*ddcp));
  memset(ddcp, 0, sizeof(*ddcp));

  int dir[3];

  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = dir2idx(dir);
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	ddcp->send_reqs[dir1] = MPI_REQUEST_NULL;
	ddcp->sendp_reqs[dir1] = MPI_REQUEST_NULL;
	ddcp->recv_reqs[dir1] = MPI_REQUEST_NULL;

	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	  nei->rank = -1;
	  continue;
	}

	nei->send_buf_size = 8;
	nei->send_buf = malloc(nei->send_buf_size * sizeof(*nei->send_buf));
	nei->n_send = 0;
	nei->rank = ddc_get_rank_nei(ddc, dir);
      }
    }
  }

  return ddcp;
}

static void
ddc_particles_queue(struct ddc_particles *ddcp, int dir[3], struct f_particle *p)
{
  struct ddcp_nei *nei = &ddcp->nei[dir2idx(dir)];

  if (nei->n_send == nei->send_buf_size) {
    // reallocate a larger buffer, doubling buffer size each time
    assert(nei->send_buf_size > 0);
    nei->send_buf_size *= 2;
    nei->send_buf = realloc(nei->send_buf, 
			    nei->send_buf_size * sizeof(*nei->send_buf));
  }
  nei->send_buf[nei->n_send++] = *p;
}

static void
ddc_particles_comm(struct ddc_particles *ddcp)
{
  int sz = sizeof(struct f_particle) / sizeof(f_real);
  int dir[3];

  // post receives for # particles we'll receive in the next step
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = dir2idx(dir);
	int dir1neg = dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	if (nei->rank < 0) {
	  ddcp->recv_reqs[dir1] = MPI_REQUEST_NULL;
	  continue;
	}
	MPI_Irecv(&nei->n_recv, 1, MPI_INT, nei->rank, 2000 + dir1neg,
		  MPI_COMM_WORLD, &ddcp->recv_reqs[dir1]);
      }
    }
  }

  // post sends for # particles and then the actual particles
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = dir2idx(dir);
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	if (nei->rank < 0) {
	  continue;
	}

	MPI_Isend(&nei->n_send, 1, MPI_INT, nei->rank, 2000 + dir1,
		  MPI_COMM_WORLD, &ddcp->send_reqs[dir1]);
	MPI_Isend(nei->send_buf, sz * nei->n_send, MPI_DOUBLE, nei->rank,
		  2100 + dir1, MPI_COMM_WORLD, &ddcp->sendp_reqs[dir1]);
      }
    }
  }

  MPI_Waitall(27, ddcp->recv_reqs, MPI_STATUSES_IGNORE);

  // calc total # of particles
  int new_n_particles = ddcp->head; // particles which stayed on this proc

  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = dir2idx(dir);
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	if (nei->rank < 0) {
	  continue;
	}
	// add the ones we're about to receive
	new_n_particles += nei->n_recv;
      }
    }
  }
  psc.f_part = REALLOC_particles(new_n_particles);

  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = dir2idx(dir);
	int dir1neg = dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	if (nei->rank < 0) {
	  continue;
	}
	MPI_Irecv(&psc.f_part[ddcp->head], sz * nei->n_recv, MPI_DOUBLE, nei->rank,
		  2100 + dir1neg, MPI_COMM_WORLD, &ddcp->recv_reqs[dir1]);
	ddcp->head += nei->n_recv;
      }
    }
  }

  MPI_Waitall(27, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(27, ddcp->send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(27, ddcp->sendp_reqs, MPI_STATUSES_IGNORE);
}

static void
create_bnd(void)
{
  struct c_bnd_ctx *c_bnd = malloc(sizeof(*c_bnd));
  memset(c_bnd, 0, sizeof(*c_bnd));

  struct ddc_params prm = {
    .comm          = MPI_COMM_WORLD,
    .mpi_type      = MPI_FIELDS_BASE_REAL,
    .size_of_type  = sizeof(fields_base_real_t),
    .n_proc        = { psc.domain.nproc[0], psc.domain.nproc[1], psc.domain.nproc[2] },
    .ilo           = { psc.ilo[0], psc.ilo[1], psc.ilo[2] },
    .ihi           = { psc.ihi[0], psc.ihi[1], psc.ihi[2] },
    .ibn           = { psc.ibn[0], psc.ibn[1], psc.ibn[2] },
    .copy_to_buf   = copy_to_buf,
    .copy_from_buf = copy_from_buf,
    .add_from_buf  = add_from_buf,
  };
  for (int d = 0; d < 3; d++) {
      if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC &&
	  psc.domain.ihi[d] - psc.domain.ilo[d] > 1) {
      prm.bc[d] = DDC_BC_PERIODIC;
    }
  }
  c_bnd->ddc = ddc_create(&prm);
  c_bnd->ddcp = ddc_particles_create(c_bnd->ddc);

  psc.bnd_data = c_bnd;
}
  
static void
c_add_ghosts(int m)
{
  if (!psc.bnd_data) {
    create_bnd();
  }
  struct c_bnd_ctx *c_bnd = psc.bnd_data;

  static int pr;
  if (!pr) {
    pr = prof_register("c_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  ddc_add_ghosts(c_bnd->ddc, m);

  prof_stop(pr);
}

static void
c_fill_ghosts(int m)
{
  if (!psc.bnd_data) {
    create_bnd();
  }
  struct c_bnd_ctx *c_bnd = psc.bnd_data;

  static int pr;
  if (!pr) {
    pr = prof_register("c_fill_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  // FIXME
  // I don't think we need as many points, and only stencil star
  // rather then box
  ddc_fill_ghosts(c_bnd->ddc, m);

  prof_stop(pr);
}

static void
c_exchange_particles(void)
{
  if (!psc.bnd_data) {
    create_bnd();
  }
  static int pr;
  if (!pr) {
    pr = prof_register("c_xchg_part", 1., 0, 0);
  }
  prof_start(pr);

  struct c_bnd_ctx *c_bnd = psc.bnd_data;
  struct ddc_particles *ddcp = c_bnd->ddcp;

  f_real xb[3], xe[3], xgb[3], xge[3], xgl[3];

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.
  // FIXME, calculate once

  for (int d = 0; d < 3; d++) {
    // ensure periodic b.c.
    if (psc.domain.ihi[d] - psc.domain.ilo[d] > 1) {
      assert(psc.domain.bnd_part[d] == 1);
    }
    xb[d] = (psc.ilo[d]-.5) * psc.dx[d];
    xe[d] = (psc.ihi[d]-.5) * psc.dx[d];
    xgb[d] = (psc.domain.ilo[d]-.5) * psc.dx[d];
    xge[d] = (psc.domain.ihi[d]-.5) * psc.dx[d];
    xgl[d] = (psc.domain.ihi[d]-psc.domain.ilo[d]) * psc.dx[d];
  }

  ddcp->head = 0;
  for (int dir1 = 0; dir1 < 27; dir1++) {
    ddcp->nei[dir1].n_send = 0;
  }
  for (int i = 0; i < psc.n_part; i++) {
    struct f_particle *p = &psc.f_part[i];
    f_real *xi = &p->xi; // slightly hacky relies on xi, yi, zi to be contiguous in
                         // the struct. FIXME
    if (xi[0] >= xb[0] && xi[0] <= xe[0] &&
	xi[1] >= xb[1] && xi[1] <= xe[1] &&
	xi[2] >= xb[2] && xi[2] <= xe[2]) {
      // fast path
      // inside domain: move into right position
      psc.f_part[ddcp->head++] = *p;
    } else {
      // slow path
      int dir[3];
      for (int d = 0; d < 3; d++) {
	if (xi[d] < xb[d]) {
	  if (xi[d] < xgb[d]) {
	    xi[d] += xgl[d];
	  }
	  dir[d] = -1;
	} else if (xi[d] > xe[d]) {
	  if (xi[d] > xge[d]) {
	    xi[d] -= xgl[d];
	  }
	  dir[d] = +1;
	} else {
	  dir[d] = 0;
	}
      }
      ddc_particles_queue(ddcp, dir, p);
    }
  }

  ddc_particles_comm(ddcp);
  psc_set_n_particles(ddcp->head);

  prof_stop(pr);
}

struct psc_bnd_ops psc_bnd_ops_c = {
  .name               = "c",
  .add_ghosts         = c_add_ghosts,
  .fill_ghosts        = c_fill_ghosts,
  .exchange_particles = c_exchange_particles,
};
