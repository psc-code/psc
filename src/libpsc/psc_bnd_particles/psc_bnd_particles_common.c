
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_c.h"

#include <mrc_domain.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <string.h>

static inline bool
at_lo_boundary(int p, int d)
{
  return ppsc->patch[p].off[d] == 0;
}

static inline bool
at_hi_boundary(int p, int d)
{
  return ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d];
}

#include "psc_bnd_particles_open.c"

// ----------------------------------------------------------------------
// ddcp_particles helpers

static void
ddcp_particles_realloc(void *_ctx, int p, int new_n_particles)
{
  struct psc_mparticles *mprts = _ctx;
  mparticles_patch_reserve(mprts, p, new_n_particles);
}

static void *
ddcp_particles_get_addr(void *_ctx, int p, int n)
{
  struct psc_mparticles *mprts = _ctx;
  particle_range_t prts = particle_range_mprts(mprts, p);
  return particle_iter_at(prts.begin, n);
}

#include "ddc_particles_inc.c"

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_setup

static void
psc_bnd_particles_sub_setup(struct psc_bnd_particles *bnd)
{
  bnd->ddcp = ddc_particles_create(bnd->psc->mrc_domain, sizeof(particle_t),
				   sizeof(particle_real_t),
				   MPI_PARTICLES_REAL,
				   ddcp_particles_realloc,
				   ddcp_particles_get_addr);

  psc_bnd_particles_open_setup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
{
  ddc_particles_destroy(bnd->ddcp);

  psc_bnd_particles_open_unsetup(bnd);
}

// ======================================================================
//

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_prep

static void
psc_bnd_particles_sub_exchange_particles_prep(struct psc_bnd_particles *bnd,
					      struct psc_mparticles *mprts, int p)
{
  particle_range_t prts = particle_range_mprts(mprts, p);
  struct ddc_particles *ddcp = bnd->ddcp;
  struct psc *psc = bnd->psc;

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.

  struct psc_patch *ppatch = &psc->patch[p];
  particle_real_t b_dxi[3] = { 1.f / ppatch->dx[0], 1.f / ppatch->dx[1], 1.f / ppatch->dx[2] };
  particle_real_t xm[3];
  int b_mx[3];
  for (int d = 0; d < 3; d++ ) {
    b_mx[d] = ppatch->ldims[d];
    xm[d] = b_mx[d] * ppatch->dx[d];
  }
  
  struct ddcp_patch *patch = &ddcp->patches[p];
  patch->head = 0;
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    patch->nei[dir1].n_send = 0;
  }
  unsigned int n_prts = particle_range_size(prts);
  for (int i = 0; i < n_prts; i++) {
    particle_t *part = particle_iter_at(prts.begin, i);
    particle_real_t *xi = &part->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
    particle_real_t *pxi = &part->pxi;
    
    int b_pos[3];
    particle_xi_get_block_pos(xi, b_dxi, b_pos);
    if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
	b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      // fast path
      // inside domain: move into right position
      *particle_iter_at(prts.begin, patch->head++) = *part;
    } else {
      // slow path
      bool drop = false;
      int dir[3];
      for (int d = 0; d < 3; d++) {
	if (b_pos[d] < 0) {
	  if (!at_lo_boundary(p, d) ||
	      psc->domain.bnd_part_lo[d] == BND_PART_PERIODIC) {
	    xi[d] += xm[d];
	    dir[d] = -1;
	    int bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi >= b_mx[d]) {
	      xi[d] = 0.;
	      dir[d] = 0;
	    }
	  } else {
	    switch (psc->domain.bnd_part_lo[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] =  -xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      break;
	    case BND_PART_ABSORBING:
	    case BND_PART_OPEN:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else if (b_pos[d] >= b_mx[d]) {
	  if (!at_hi_boundary(p, d) ||
	      psc->domain.bnd_part_hi[d] == BND_PART_PERIODIC) {
	    xi[d] -= xm[d];
	    dir[d] = +1;
	    int bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi < 0) {
	      xi[d] = 0.;
	    }
	  } else {
	    switch (psc->domain.bnd_part_hi[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] = 2.f * xm[d] - xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      int bi = particle_real_fint(xi[d] * b_dxi[d]);
	      if (bi >= b_mx[d]) {
		xi[d] *= (1. - 1e-6);
	      }
	      break;
	    case BND_PART_ABSORBING:
	    case BND_PART_OPEN:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else {
	  // computational bnd
	  dir[d] = 0;
	}
	if (!drop) {
	  if (xi[d] < 0.f && xi[d] > -1e-6f) {
	    //	    mprintf("d %d xi %g\n", d, xi[d]);
	    xi[d] = 0.f;
	  }
	  assert(xi[d] >= 0.f);
	  assert(xi[d] <= xm[d]);
	}
      }
      if (!drop) {
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	  *particle_iter_at(prts.begin, patch->head++) = *part;
	} else {
	  ddc_particles_queue(ddcp, patch, dir, part);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_post

static void
psc_bnd_particles_sub_exchange_particles_post(struct psc_bnd_particles *bnd,
					      struct psc_mparticles *mprts, int p)
{
  struct ddc_particles *ddcp = bnd->ddcp;
  struct ddcp_patch *patch = &ddcp->patches[p];
  mparticles_patch_resize(mprts, p, patch->head);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

static void
psc_bnd_particles_sub_exchange_particles(struct psc_bnd_particles *bnd, struct psc_mparticles *particles_base)
{
  struct psc *psc = bnd->psc;

  struct psc_mparticles *particles = psc_mparticles_get_as(particles_base, PARTICLE_TYPE, 0);

  static int pr_A, pr_B, pr_C;
  if (!pr_A) {
    pr_A = prof_register("xchg_prep", 1., 0, 0);
    pr_B = prof_register("xchg_comm", 1., 0, 0);
    pr_C = prof_register("xchg_post", 1., 0, 0);
  }
  
  prof_start(pr_A);

  struct ddc_particles *ddcp = bnd->ddcp;

  // FIXME we should make sure (assert) we don't quietly drop particle which left
  // in the invariant direction

  prof_restart(pr_time_step_no_comm);
#pragma omp parallel for
  psc_foreach_patch(psc, p) {
    psc_balance_comp_time_by_patch[p] -= MPI_Wtime();
    psc_bnd_particles_sub_exchange_particles_prep(bnd, particles, p);
    psc_balance_comp_time_by_patch[p] += MPI_Wtime();
  }
  prof_stop(pr_time_step_no_comm);

  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);
  prof_stop(pr_B);

  prof_start(pr_C);
  for (int p = 0; p < particles->nr_patches; p++) {
    psc_bnd_particles_sub_exchange_particles_post(bnd, particles, p);
  }
  prof_stop(pr_C);

  //struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, "c", JXI, JXI + 3);
  //psc_bnd_particles_open_boundary(bnd, particles, mflds);
  //psc_mfields_put_as(mflds, psc->flds, JXI, JXI + 3);

  psc_mparticles_put_as(particles, particles_base, 0);
}

