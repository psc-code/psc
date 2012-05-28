
#include <mrc_domain.h>
#include <mrc_profile.h>
#include <string.h>

struct psc_bnd_sub {
};

#define to_psc_bnd_sub(bnd) ((struct psc_bnd_sub *)((bnd)->obj.subctx))

// ----------------------------------------------------------------------
// ddcp_particles helpers

static void
ddcp_particles_realloc(void *_ctx, int p, int new_n_particles)
{
  mparticles_t *particles = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
  particles_realloc(prts, new_n_particles);
}

static void *
ddcp_particles_get_addr(void *_ctx, int p, int n)
{
  mparticles_t *particles = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
  return particles_get_one(prts, n);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_setup

static void
psc_bnd_sub_setup(struct psc_bnd *bnd)
{
  psc_bnd_setup_super(bnd);
  bnd->ddcp = ddc_particles_create(bnd->ddc, sizeof(particle_t),
				   sizeof(particle_real_t),
				   MPI_PARTICLES_REAL,
				   ddcp_particles_realloc,
				   ddcp_particles_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_unsetup

static void
psc_bnd_sub_unsetup(struct psc_bnd *bnd)
{
  ddc_particles_destroy(bnd->ddcp);
}

// ======================================================================
//
// ----------------------------------------------------------------------
// find_block_position

static inline void
find_block_position(int b_pos[3], particle_real_t xi[3], particle_real_t b_dxi[3])
{
  for (int d = 0; d < 3; d++) {
    b_pos[d] = particle_real_fint(xi[d] * b_dxi[d]);
  }
}

// ----------------------------------------------------------------------
// find_block_indices_count

static inline void
find_block_indices_count(unsigned int *b_idx, unsigned int *b_cnts, struct psc_particles *prts,
			 int off)
{
  struct psc_particles_single *sngl = psc_particles_single(prts);
  int *b_mx = sngl->b_mx;
  for (int i = off; i < prts->n_part; i++) {
    particle_t *part = particles_get_one(prts, i);
    int b_pos[3];
    find_block_position(b_pos, &part->xi, sngl->b_dxi);
    assert(b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
	   b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	   b_pos[2] >= 0 && b_pos[2] < b_mx[2]);
    b_idx[i] = (b_pos[2] * b_mx[1] + b_pos[1]) * b_mx[0] + b_pos[0];
    b_cnts[b_idx[i]]++;
  }
}

// ----------------------------------------------------------------------
// find_block_indices_count_reorder

static void __unused
find_block_indices_count_reorder(struct psc_particles *prts)
{
  struct psc_particles_single *sngl = psc_particles_single(prts);
  unsigned int cnt = prts->n_part;
  int *b_mx = sngl->b_mx;
  memset(sngl->b_cnt, 0, (sngl->nr_blocks + 1) * sizeof(*sngl->b_cnt));
  for (int i = 0; i < prts->n_part; i++) {
    particle_t *part = particles_get_one(prts, i);
    int b_pos[3];
    find_block_position(b_pos, &part->xi, sngl->b_dxi);
    if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
	b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      sngl->b_idx[i] = (b_pos[2] * b_mx[1] + b_pos[1]) * b_mx[0] + b_pos[0];
    } else { // out of bounds
      sngl->b_idx[i] = sngl->nr_blocks;
      assert(cnt < sngl->n_alloced);
      *particles_get_one(prts, cnt) = *part;
      cnt++;
    }
    sngl->b_cnt[sngl->b_idx[i]]++;
  }
}

static void __unused
count_and_reorder_to_back(struct psc_particles *prts)
{
  struct psc_particles_single *sngl = psc_particles_single(prts);

  memset(sngl->b_cnt, 0, (sngl->nr_blocks + 1) * sizeof(*sngl->b_cnt));
  unsigned int cnt = prts->n_part;
  for (int i = 0; i < prts->n_part; i++) {
    if (sngl->b_idx[i] == sngl->nr_blocks) {
      assert(cnt < sngl->n_alloced);
      *particles_get_one(prts, cnt) = *particles_get_one(prts, i);
      cnt++;
    }
    sngl->b_cnt[sngl->b_idx[i]]++;
  }
}

static void __unused
reorder_to_back(struct psc_particles *prts)
{
  struct psc_particles_single *sngl = psc_particles_single(prts);
  unsigned int cnt = prts->n_part;
  for (int i = 0; i < prts->n_part; i++) {
    if (sngl->b_idx[i] == sngl->nr_blocks) {
      assert(cnt < sngl->n_alloced);
      *particles_get_one(prts, cnt) = *particles_get_one(prts, i);
      cnt++;
    }
  }
}

// ----------------------------------------------------------------------
// count_block_indices

static inline void
count_block_indices(unsigned int *b_cnts, unsigned int *b_idx, int n_part, int off)
{
  for (int i = off; i < n_part; i++) {
    b_cnts[b_idx[i]]++;
  }
}

// ----------------------------------------------------------------------
// exclusive_scan

static inline void
exclusive_scan(unsigned int *b_cnts, int n)
{
  unsigned int sum = 0;
  for (int i = 0; i < n; i++) {
    unsigned int cnt = b_cnts[i];
    b_cnts[i] = sum;
    sum += cnt;
  }
}

static inline void
xchg_append(struct psc_particles *prts, void *patch_ctx, particle_t *prt)
{
  struct ddcp_patch *ddcp_patch = patch_ctx;
  *particles_get_one(prts, ddcp_patch->head++) = *prt;
}

static inline int *
get_b_mx(struct psc_particles *prts)
{
  return ppsc->patch[prts->p].ldims;
}

static inline particle_real_t *
get_b_dxi(struct psc_particles *prts)
{
  static particle_real_t b_dxi[3];
  if (!b_dxi[0]) {
    for (int d = 0; d < 3; d++) {
      b_dxi[d] = 1.f / ppsc->dx[d];
    }
  }
  return b_dxi;
}

// ----------------------------------------------------------------------
// exchange particles

static void
exchange_particles(struct psc_bnd *bnd, struct psc_mparticles *particles)
{
  struct psc *psc = bnd->psc;
  struct ddc_particles *ddcp = bnd->ddcp;

  static int pr_A, pr_B, pr_D;
  if (!pr_A) {
    pr_A = prof_register("xchg_prep_" PARTICLE_TYPE, 1., 0, 0);
    pr_B = prof_register("xchg_comm_" PARTICLE_TYPE, 1., 0, 0);
    pr_D = prof_register("xchg_post_" PARTICLE_TYPE, 1., 0, 0);
  }

  prof_start(pr_A);

  // FIXME we should make sure (assert) we don't quietly drop particle which left
  // in the invariant direction

  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    int n_send = psc_particles_single(prts)->n_send;

    struct psc_patch *patch = &psc->patch[p];
    particle_real_t xm[3];
    for (int d = 0; d < 3; d++) {
      xm[d] = patch->ldims[d] * psc->dx[d];
    }
    particle_real_t *b_dxi = get_b_dxi(prts);
    int *b_mx = get_b_mx(prts);

    struct ddcp_patch *ddcp_patch = &ddcp->patches[p];
    ddcp_patch->head = prts->n_part - n_send;
    for (int dir1 = 0; dir1 < N_DIR; dir1++) {
      ddcp_patch->nei[dir1].n_send = 0;
    }
    for (int n = ddcp_patch->head; n < prts->n_part; n++) {
      particle_t *prt = particles_get_one(prts, n);
      particle_real_t *xi = &prt->xi;
      particle_real_t *pxi = &prt->pxi;
      
      bool drop = false;
      int dir[3];
      for (int d = 0; d < 3; d++) {
	int bi = particle_real_fint(xi[d] * b_dxi[d]);
	if (bi < 0) {
	  // FIXME, assumes every patch has same dimensions
	  if (patch->off[d] != 0 || psc->domain.bnd_part_lo[d] == BND_PART_PERIODIC) {
	    xi[d] += xm[d];
	    dir[d] = -1;
	    bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi >= b_mx[d]) {
	      xi[d] = 0.;
	      dir[d] = 0;
	    }
	  } else {
	    switch (psc->domain.bnd_part_lo[d]) {
	    case BND_PART_REFLECTING:
	      xi[d]  = -xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      break;
	    case BND_PART_ABSORBING:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else if (bi >= b_mx[d]) {
	  if (patch->off[d] + patch->ldims[d] != psc->domain.gdims[d] ||
	      psc->domain.bnd_part_hi[d] == BND_PART_PERIODIC) {
	    xi[d] -= xm[d];
	    dir[d] = +1;
	    bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi < 0) {
	      xi[d] = 0.;
	    }
	  } else {
	    switch (psc->domain.bnd_part_hi[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] = 2.f * xm[d] - xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      bi = particle_real_fint(xi[d] * b_dxi[d]);
	      if (bi >= b_mx[d]) {
		xi[d] *= (1. - 1e-6);
	      }
	      break;
	    case BND_PART_ABSORBING:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else {
	  dir[d] = 0;
	}
	assert(xi[d] >= 0.f);
	assert(xi[d] <= xm[d]);
      }
      if (!drop) {
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	  xchg_append(prts, ddcp_patch, prt);
	} else {
	  ddc_particles_queue(ddcp, ddcp_patch, dir, prt);
	}
      }
    }
  }
  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);
  prof_stop(pr_B);

  prof_start(pr_D);
  psc_foreach_patch(psc, p) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct ddcp_patch *patch = &ddcp->patches[p];
    prts->n_part = patch->head;
  }
  prof_stop(pr_D);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_exchange_particles

static void
psc_bnd_sub_exchange_particles(struct psc_bnd *bnd, mparticles_base_t *particles_base)
{
  static int pr, pr_A, pr_B, pr_C, pr_D;
  if (!pr) {
    pr   = prof_register("xchg_parts", 1., 0, 0);
    pr_A = prof_register("xchg_bidx_mv", 1., 0, 0);
    pr_B = prof_register("xchg_actual", 1., 0, 0);
    pr_C = prof_register("xchg_upd_bidx", 1., 0, 0);
    pr_D = prof_register("xchg_reorder", 1., 0, 0);
  }
  
  prof_start(pr);

  struct psc *psc = bnd->psc;
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  prof_start(pr_A);
  psc_foreach_patch(psc, p) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_single *sngl = psc_particles_single(prts);

    if (1) {
      //      find_block_indices_count_reorderx(prts);
      count_and_reorder_to_back(prts);
    }
    sngl->n_part_save = prts->n_part;
    sngl->n_send = sngl->b_cnt[sngl->nr_blocks];
    prts->n_part += sngl->n_send;
  }
  prof_stop(pr_A);

  prof_start(pr_B);
  exchange_particles(bnd, particles);
  prof_stop(pr_B);

  prof_start(pr_C);
  psc_foreach_patch(psc, p) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_single *sngl = psc_particles_single(prts);

    find_block_indices_count(sngl->b_idx, sngl->b_cnt, prts, sngl->n_part_save);
    exclusive_scan(sngl->b_cnt, sngl->nr_blocks + 1);
  }
  prof_stop(pr_C);

  prof_start(pr_D);
  psc_foreach_patch(psc, p) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_single *sngl = psc_particles_single(prts);

    psc_particles_reorder(prts, sngl->b_idx, sngl->b_cnt);
    prts->n_part = sngl->b_cnt[sngl->nr_blocks - 1];
  }
  prof_stop(pr_D);

  psc_mparticles_put_cf(particles, particles_base, 0);
  prof_stop(pr);
}

