
#include "psc.h"
#include <mrc_profile.h>
#include <mrc_ddc.h>
#include <mrc_domain.h>

#include <mpi.h>
#include <string.h>

// ======================================================================
// C bnd

struct c_bnd_ctx {
  struct mrc_ddc *ddc;
  struct ddc_particles *ddcp;
};

static void
copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_base_t *flds = ctx;
  fields_base_t *pf = &flds->f[p];
  fields_base_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf, m - mb, ix,iy,iz) = F3_BASE(pf, m, ix,iy,iz);
	}
      }
    }
  }
}

static void
add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_base_t *flds = ctx;
  fields_base_t *pf = &flds->f[p];
  fields_base_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_BASE(pf, m, ix,iy,iz) += MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

static void
copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_base_t *flds = ctx;
  fields_base_t *pf = &flds->f[p];
  fields_base_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_BASE(pf, m, ix,iy,iz) = MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = copy_to_buf,
  .copy_from_buf = copy_from_buf,
  .add_from_buf  = add_from_buf,
};

// ======================================================================

struct ddcp_nei {
  particle_base_t *send_buf;
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
ddc_particles_create(struct mrc_ddc *ddc)
{
  struct ddc_particles *ddcp = malloc(sizeof(*ddcp));
  memset(ddcp, 0, sizeof(*ddcp));

  int dir[3];

  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = mrc_ddc_dir2idx(dir);
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
	int dummy;
	mrc_ddc_get_nei_rank_patch(ddc, 0, dir, &nei->rank, &dummy);
      }
    }
  }

  return ddcp;
}

static void
ddc_particles_queue(struct ddc_particles *ddcp, int dir[3], particle_base_t *p)
{
  struct ddcp_nei *nei = &ddcp->nei[mrc_ddc_dir2idx(dir)];

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
ddc_particles_comm(struct ddc_particles *ddcp, mparticles_base_t *particles)
{
  particles_base_t *pp = &particles->p[0];

  int sz = sizeof(particle_base_t) / sizeof(particle_base_real_t);
  int dir[3];

  // post receives for # particles we'll receive in the next step
  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = mrc_ddc_dir2idx(dir);
	int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
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
	int dir1 = mrc_ddc_dir2idx(dir);
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	if (nei->rank < 0) {
	  continue;
	}

	MPI_Isend(&nei->n_send, 1, MPI_INT, nei->rank, 2000 + dir1,
		  MPI_COMM_WORLD, &ddcp->send_reqs[dir1]);
	MPI_Isend(nei->send_buf, sz * nei->n_send, MPI_PARTICLES_BASE_REAL, nei->rank,
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
	int dir1 = mrc_ddc_dir2idx(dir);
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	if (nei->rank < 0) {
	  continue;
	}
	// add the ones we're about to receive
	new_n_particles += nei->n_recv;
      }
    }
  }
  particles_base_realloc(pp, new_n_particles);

  for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
    for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
      for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	int dir1 = mrc_ddc_dir2idx(dir);
	int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	struct ddcp_nei *nei = &ddcp->nei[dir1];
	if (nei->rank < 0) {
	  continue;
	}
	MPI_Irecv(&pp->particles[ddcp->head], sz * nei->n_recv,
		  MPI_PARTICLES_BASE_REAL,
		  nei->rank, 2100 + dir1neg, MPI_COMM_WORLD, &ddcp->recv_reqs[dir1]);
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

  c_bnd->ddc = mrc_domain_create_ddc(psc.mrc_domain);
  mrc_ddc_set_funcs(c_bnd->ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(c_bnd->ddc, "ibn", psc.ibn);
  mrc_ddc_set_param_int(c_bnd->ddc, "max_n_fields", 6);
  mrc_ddc_set_param_int(c_bnd->ddc, "size_of_type", sizeof(fields_base_real_t));
  mrc_ddc_setup(c_bnd->ddc);

  c_bnd->ddcp = ddc_particles_create(c_bnd->ddc);

  psc.bnd_data = c_bnd;
}
  
static void
c_add_ghosts(mfields_base_t *flds, int mb, int me)
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

  mrc_ddc_add_ghosts(c_bnd->ddc, mb, me, flds);

  prof_stop(pr);
}

static void
c_fill_ghosts(mfields_base_t *flds, int mb, int me)
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
  mrc_ddc_fill_ghosts(c_bnd->ddc, mb, me, flds);

  prof_stop(pr);
}

static void
c_exchange_particles(mparticles_base_t *particles)
{
  particles_base_t *pp = &particles->p[0];

  if (!psc.bnd_data) {
    create_bnd();
  }
  static int pr, pr_A, pr_B;
  if (!pr) {
    pr = prof_register("c_xchg_part", 1., 0, 0);
    pr_A = prof_register("c_xchg_part_A", 1., 0, 0);
    pr_B = prof_register("c_xchg_part_B", 1., 0, 0);
  }
  prof_start(pr);
  prof_start(pr_A);

  struct c_bnd_ctx *c_bnd = psc.bnd_data;
  struct ddc_particles *ddcp = c_bnd->ddcp;

  f_real xb[3], xe[3], xgb[3], xge[3], xgl[3];

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.
  // FIXME, calculate once

  struct psc_patch *patch = &psc.patch[0];
  for (int d = 0; d < 3; d++) {
    xb[d] = (patch->off[d]-.5) * psc.dx[d];
    if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC) {
      xgb[d] = -.5 * psc.dx[d];
    } else {
      xgb[d] = 0.;
      if (patch->off[d] == 0) {
	xb[d] = xgb[d];
      }
    }

    xe[d] = (patch->off[d] + patch->ldims[d] - .5) * psc.dx[d];
    if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC) {
      xge[d] = (psc.domain.gdims[d]-.5) * psc.dx[d];
    } else {
      xge[d] = (psc.domain.gdims[d]-1) * psc.dx[d];
      if (patch->off[d] + patch->ldims[d] == psc.domain.gdims[d]) {
	xe[d] = xge[d];
      }
    }

    xgl[d] = xge[d] - xgb[d];
  }

  ddcp->head = 0;
  for (int dir1 = 0; dir1 < 27; dir1++) {
    ddcp->nei[dir1].n_send = 0;
  }
  for (int i = 0; i < pp->n_part; i++) {
    particle_base_t *p = particles_base_get_one(pp, i);
    particle_base_real_t *xi = &p->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
    particle_base_real_t *pxi = &p->pxi;
    if (xi[0] >= xb[0] && xi[0] <= xe[0] &&
	xi[1] >= xb[1] && xi[1] <= xe[1] &&
	xi[2] >= xb[2] && xi[2] <= xe[2]) {
      // fast path
      // inside domain: move into right position
      pp->particles[ddcp->head++] = *p;
    } else {
      // slow path
      int dir[3];
      for (int d = 0; d < 3; d++) {
	if (xi[d] < xb[d]) {
	  if (xi[d] < xgb[d]) {
	    switch (psc.domain.bnd_part[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] = 2.f * xgb[d] - xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      break;
	    case BND_PART_PERIODIC:
	      xi[d] += xgl[d];
	      dir[d] = -1;
	      break;
	    default:
	      assert(0);
	    }
	  } else {
	    // computational bnd
	    dir[d] = -1;
	  }
	} else if (xi[d] > xe[d]) {
	  if (xi[d] > xge[d]) {
	    switch (psc.domain.bnd_part[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] = 2.f * xge[d] - xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      break;
	    case BND_PART_PERIODIC:
	      xi[d] -= xgl[d];
	      dir[d] = +1;
	      break;
	    default:
	      assert(0);
	    }
	  } else {
	    dir[d] = +1;
	  }
	} else {
	  // computational bnd
	  dir[d] = 0;
	}
      }
      if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	pp->particles[ddcp->head++] = *p;
      } else {
	ddc_particles_queue(ddcp, dir, p);
      }
    }
  }

  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);
  psc_set_n_particles(pp, ddcp->head);
  prof_stop(pr_B);

  prof_stop(pr);
}

struct psc_bnd_ops psc_bnd_ops_c = {
  .name               = "c",
  .add_ghosts         = c_add_ghosts,
  .fill_ghosts        = c_fill_ghosts,
  .exchange_particles = c_exchange_particles,
};
