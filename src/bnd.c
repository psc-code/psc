
#include "psc.h"
#include "psc_bnd_private.h"
#include <mrc_profile.h>
#include <mrc_ddc.h>
#include <mrc_domain.h>

#include <mpi.h>
#include <string.h>

// ======================================================================

#define N_DIR (27)

struct ddcp_nei {
  particle_base_t *send_buf;
  int n_send;
  int n_recv;
  int rank;
  int patch;
  int send_buf_size;
};

struct ddcp_patch {
  int head;
  struct ddcp_nei nei[N_DIR];
};

struct ddc_particles {
  struct ddcp_patch *patches;
  MPI_Request *send_reqs;
  MPI_Request *sendp_reqs;
  MPI_Request *recv_reqs;
};

static struct ddc_particles *
ddc_particles_create(struct mrc_ddc *ddc)
{
  struct ddc_particles *ddcp = malloc(sizeof(*ddcp));
  memset(ddcp, 0, sizeof(*ddcp));

  ddcp->patches = calloc(psc.nr_patches, sizeof(*ddcp->patches));
  ddcp->send_reqs  = calloc(psc.nr_patches * N_DIR, sizeof(*ddcp->send_reqs));
  ddcp->sendp_reqs = calloc(psc.nr_patches * N_DIR, sizeof(*ddcp->sendp_reqs));
  ddcp->recv_reqs  = calloc(psc.nr_patches * N_DIR, sizeof(*ddcp->recv_reqs));
  foreach_patch(p) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  ddcp->send_reqs[dir1 + p * N_DIR] = MPI_REQUEST_NULL;
	  ddcp->sendp_reqs[dir1 + p * N_DIR] = MPI_REQUEST_NULL;
	  ddcp->recv_reqs[dir1 + p * N_DIR] = MPI_REQUEST_NULL;
	  
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    nei->rank = -1;
	    continue;
	  }
	  
	  nei->send_buf_size = 8;
	  nei->send_buf = malloc(nei->send_buf_size * sizeof(*nei->send_buf));
	  nei->n_send = 0;
	  mrc_ddc_get_nei_rank_patch(ddc, p, dir, &nei->rank, &nei->patch);
	}
      }
    }
  }

  return ddcp;
}

static void
ddc_particles_queue(struct ddcp_patch *patch, int dir[3], particle_base_t *p)
{
  struct ddcp_nei *nei = &patch->nei[mrc_ddc_dir2idx(dir)];

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
  int sz = sizeof(particle_base_t) / sizeof(particle_base_real_t);
  int dir[3];

  // post receives for # particles we'll receive in the next step
  foreach_patch(p) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
#if 0
	  mprintf(":%d irecv # from %d:%d tag %d\n", p, nei->rank, nei->patch,
		  dir1neg + p * N_DIR);
#endif
	  MPI_Irecv(&nei->n_recv, 1, MPI_INT, nei->rank, dir1neg + p * N_DIR,
		    MPI_COMM_WORLD, &ddcp->recv_reqs[dir1 + p * N_DIR]);
	}
      }
    }
  }

  // post sends for # particles and then the actual particles
  foreach_patch(p) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
#if 0
	  mprintf(":%d isend # to   %d:%d tag %d: len %d\n", p, nei->rank, nei->patch,
		  dir1 + nei->patch * N_DIR, nei->n_send);
	  mprintf(":%d isend P to   %d:%d tag %d\n", p, nei->rank, nei->patch,
		  dir1 + nei->patch * N_DIR);
#endif
	  MPI_Isend(&nei->n_send, 1, MPI_INT,
		    nei->rank, dir1 + nei->patch * N_DIR, MPI_COMM_WORLD,
		    &ddcp->send_reqs[dir1 + p * N_DIR]);
	  MPI_Isend(nei->send_buf, sz * nei->n_send, MPI_PARTICLES_BASE_REAL,
		    nei->rank, dir1 + nei->patch * N_DIR, MPI_COMM_WORLD,
		    &ddcp->sendp_reqs[dir1 + p * N_DIR]);
	}
      }
    }
  }

  MPI_Waitall(psc.nr_patches * N_DIR, ddcp->recv_reqs, MPI_STATUSES_IGNORE);

  // calc total # of particles
  foreach_patch(p) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    particles_base_t *pp = &particles->p[p];
    int new_n_particles = patch->head; // particles which stayed on this proc

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
	  // add the ones we're about to receive
	  new_n_particles += nei->n_recv;
	}
      }
    }
    particles_base_realloc(pp, new_n_particles);
  }

  // post particle receives
  foreach_patch(p) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    particles_base_t *pp = &particles->p[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0) {
	    continue;
	  }
#if 0
	  mprintf(":%d irecv P from %d:%d tag %d len %d\n", p, nei->rank, nei->patch,
		  dir1neg + p * N_DIR, nei->n_recv);
#endif
	  MPI_Irecv(&pp->particles[patch->head], sz * nei->n_recv,
		    MPI_PARTICLES_BASE_REAL,
		    nei->rank, dir1neg + p * N_DIR, MPI_COMM_WORLD,
		    &ddcp->recv_reqs[dir1 + p * N_DIR]);
	  patch->head += nei->n_recv;
	}
      }
    }
  }

  MPI_Waitall(psc.nr_patches * N_DIR, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(psc.nr_patches * N_DIR, ddcp->send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(psc.nr_patches * N_DIR, ddcp->sendp_reqs, MPI_STATUSES_IGNORE);
}

#define to_psc_bnd_c(bnd) ((struct psc_bnd_c *)((bnd)->obj.subctx))

void
create_bnd(void)
{
  struct psc_bnd *bnd = psc.bnd;
  struct psc_bnd_c *c_bnd = to_psc_bnd_c(bnd);

  c_bnd->ddcp = ddc_particles_create(c_bnd->ddc);
}
  
static void
c_exchange_particles(mparticles_base_t *particles)
{
  struct psc_bnd *bnd = psc.bnd;
  struct psc_bnd_c *c_bnd = to_psc_bnd_c(bnd);

  static int pr, pr_A, pr_B;
  if (!pr) {
    pr = prof_register("c_xchg_part", 1., 0, 0);
    pr_A = prof_register("c_xchg_part_A", 1., 0, 0);
    pr_B = prof_register("c_xchg_part_B", 1., 0, 0);
  }
  prof_start(pr);
  prof_start(pr_A);

  struct ddc_particles *ddcp = c_bnd->ddcp;

  f_real xb[3], xe[3], xgb[3], xge[3], xgl[3];

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.
  // FIXME, calculate once

  foreach_patch(p) {
    struct psc_patch *psc_patch = &psc.patch[p];
    particles_base_t *pp = &particles->p[p];

    for (int d = 0; d < 3; d++) {
      xb[d] = (psc_patch->off[d]-.5) * psc.dx[d];
      if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC) {
	xgb[d] = -.5 * psc.dx[d];
      } else {
	xgb[d] = 0.;
	if (psc_patch->off[d] == 0) {
	  xb[d] = xgb[d];
	}
      }
      
      xe[d] = (psc_patch->off[d] + psc_patch->ldims[d] - .5) * psc.dx[d];
      if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC) {
	xge[d] = (psc.domain.gdims[d]-.5) * psc.dx[d];
      } else {
	xge[d] = (psc.domain.gdims[d]-1) * psc.dx[d];
	if (psc_patch->off[d] + psc_patch->ldims[d] == psc.domain.gdims[d]) {
	  xe[d] = xge[d];
	}
      }
      
      xgl[d] = xge[d] - xgb[d];
    }

    struct ddcp_patch *patch = &ddcp->patches[p];
    patch->head = 0;
    for (int dir1 = 0; dir1 < N_DIR; dir1++) {
      patch->nei[dir1].n_send = 0;
    }
    for (int i = 0; i < pp->n_part; i++) {
      particle_base_t *part = particles_base_get_one(pp, i);
      particle_base_real_t *xi = &part->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
      particle_base_real_t *pxi = &part->pxi;
      if (xi[0] >= xb[0] && xi[0] <= xe[0] &&
	  xi[1] >= xb[1] && xi[1] <= xe[1] &&
	  xi[2] >= xb[2] && xi[2] <= xe[2]) {
	// fast path
	// inside domain: move into right position
	pp->particles[patch->head++] = *part;
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
	  pp->particles[patch->head++] = *part;
	} else {
	  ddc_particles_queue(patch, dir, part);
	}
      }
    }
  }
  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);

  foreach_patch(p) {
    particles_base_t *pp = &particles->p[p];
    struct ddcp_patch *patch = &ddcp->patches[p];
    psc_set_n_particles(pp, patch->head);
  }
  prof_stop(pr_B);

  prof_stop(pr);
}

struct _psc_bnd_ops psc_bnd_ops_c = {
  .name               = "c",
  .exchange_particles = c_exchange_particles,
};
