
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_c.h"

#include <mrc_domain.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <string.h>

static void
copy_to_mrc_fld(struct mrc_fld *m3, struct psc_mfields *flds)
{
  psc_foreach_patch(ppsc, p) {
    struct psc_fields *pf = psc_mfields_get_patch(flds, p);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < flds->nr_fields; m++) {
	MRC_M3(m3p,m, ix,iy,iz) = F3_C(pf,m, ix,iy,iz);
      }
    } mrc_fld_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

static void
average_9_point(struct psc_mfields *mflds_av, struct psc_mfields *mflds)
{
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields *flds_av = psc_mfields_get_patch(mflds_av, p);
    
    // at lower z bnd only FIXME
    for (int m = 0; m < mflds->nr_fields; m++) {
      for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	F3_C(flds_av, m, 0,iy,0) = (F3_C(flds, m, 0,iy-1,0) +
				    F3_C(flds, m, 0,iy-1,1) +
				    F3_C(flds, m, 0,iy-1,2) +
				    F3_C(flds, m, 0,iy  ,0) +
				    F3_C(flds, m, 0,iy  ,1) +
				    F3_C(flds, m, 0,iy  ,2) +
				    F3_C(flds, m, 0,iy+1,0) +
				    F3_C(flds, m, 0,iy+1,1) +
				    F3_C(flds, m, 0,iy+1,2)) / 9.;
      }
    }
  }
}

static void
average_in_time(struct psc_mfields *mflds_av, struct psc_mfields *mflds_last,
		bool first_time)
{
  const double rr = 0.01; // FIXME
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_av = psc_mfields_get_patch(mflds_av, p);
    struct psc_fields *flds_last = psc_mfields_get_patch(mflds_last, p);
    
    if (first_time) {
      for (int m = 0; m < mflds_av->nr_fields; m++) {
	for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	  F3_C(flds_last, m, 0,iy,0) = F3_C(flds_av, m, 0,iy,0);
	}
      }
    }

    for (int m = 0; m < mflds_av->nr_fields; m++) {
      for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	F3_C(flds_av, m, 0,iy,0) = 
	  rr * F3_C(flds_av, m, 0,iy,0)
	  + (1.0 - rr) * F3_C(flds_last, m, 0,iy,0);
	F3_C(flds_last, m, 0,iy,0) = F3_C(flds_av, m, 0,iy,0);
      }
    }
  }
}

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
// psc_bnd_particles_sub_setup

static void
psc_bnd_particles_sub_setup(struct psc_bnd_particles *bnd)
{
  bnd->ddcp = ddc_particles_create(bnd->psc->mrc_domain, sizeof(particle_t),
				   sizeof(particle_real_t),
				   MPI_PARTICLES_REAL,
				   ddcp_particles_realloc,
				   ddcp_particles_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
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
// psc_bnd_particles_sub_exchange_particles_prep

static void
psc_bnd_particles_sub_exchange_particles_prep(struct psc_bnd_particles *bnd, struct psc_particles *prts)
{
  struct ddc_particles *ddcp = bnd->ddcp;
  struct psc *psc = bnd->psc;

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.

  struct psc_patch *ppatch = &psc->patch[prts->p];
  particle_real_t b_dxi[3] = { 1.f / ppatch->dx[0], 1.f / ppatch->dx[1], 1.f / ppatch->dx[2] };
  particle_real_t xm[3];
  int b_mx[3];
  for (int d = 0; d < 3; d++ ) {
    if (psc->domain.bnd_part_hi[d] == BND_PART_REFLECTING &&
	!psc->prm.gdims_in_terms_of_cells &&
	ppatch->off[d] + ppatch->ldims[d] == psc->domain.gdims[d]) {
      b_mx[d] = ppatch->ldims[d] - 1;
    } else {
      b_mx[d] = ppatch->ldims[d];
    }
    xm[d] = b_mx[d] * ppatch->dx[d];
  }
  
  struct ddcp_patch *patch = &ddcp->patches[prts->p];
  patch->head = 0;
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    patch->nei[dir1].n_send = 0;
  }
  for (int i = 0; i < prts->n_part; i++) {
    particle_t *part = particles_get_one(prts, i);
    particle_real_t *xi = &part->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
    
    int b_pos[3];
    find_block_position(b_pos, xi, b_dxi);
    particle_real_t *pxi = &part->pxi;
    if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
	b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      // fast path
      // inside domain: move into right position
      *particles_get_one(prts, patch->head++) = *part;
    } else {
      // slow path
      bool drop = false;
      int dir[3];
      for (int d = 0; d < 3; d++) {
	if (b_pos[d] < 0) {
	  if (ppatch->off[d] > 0 ||
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
	  if (ppatch->off[d] + ppatch->ldims[d] < psc->domain.gdims[d] ||
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
	  *particles_get_one(prts, patch->head++) = *part;
	} else {
	  ddc_particles_queue(ddcp, patch, dir, part);
	}
      }
    }
  }
}

static void
psc_bnd_particles_open_boundary(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct psc_bnd *flds_bnd = psc_bnd_create(psc_mparticles_comm(mprts));
  psc_bnd_set_type(flds_bnd, "c");
  psc_bnd_set_psc(flds_bnd, ppsc);
  psc_bnd_setup(flds_bnd);

  static int itime=0;
  //  int ninco=0;
  //  static int seed;
  if (!itime){
   srandom(rank);
   //   seed=random();
  }

  struct psc_output_fields_item *item_n =
    psc_output_fields_item_create(psc_mparticles_comm(mprts));
  psc_output_fields_item_set_type(item_n, "n_1st_double");
  psc_output_fields_item_set_psc_bnd(item_n, flds_bnd);
  psc_output_fields_item_setup(item_n);

  struct psc_mfields *mflds_n = psc_output_fields_item_create_mfields(item_n);
  struct psc_mfields *mflds_n_av = psc_output_fields_item_create_mfields(item_n);

  struct psc_output_fields_item *item_v =
  psc_output_fields_item_create(psc_mparticles_comm(mprts));
  psc_output_fields_item_set_type(item_v, "v_1st_double");
  psc_output_fields_item_set_psc_bnd(item_v, flds_bnd);
  psc_output_fields_item_setup(item_v);
  struct psc_mfields *mflds_v = psc_output_fields_item_create_mfields(item_v);
  struct psc_mfields *mflds_v_av = psc_output_fields_item_create_mfields(item_v);

  struct psc_output_fields_item *item_t =
  psc_output_fields_item_create(psc_mparticles_comm(mprts));
  psc_output_fields_item_set_type(item_t, "T_1st_double");
  psc_output_fields_item_set_psc_bnd(item_t, flds_bnd);
  psc_output_fields_item_setup(item_t);
  struct psc_mfields *mflds_t = psc_output_fields_item_create_mfields(item_t);
  struct psc_mfields *mflds_t_av = psc_output_fields_item_create_mfields(item_t);
  struct psc_mfields *mflds_ipr = psc_output_fields_item_create_mfields(item_t);

  static bool first_time = true;
  static struct psc_mfields *mflds_n_last, *mflds_v_last, *mflds_t_last, *mflds_n_in;
  if (!mflds_n_last) {
   mflds_n_last = psc_output_fields_item_create_mfields(item_n);
   mflds_v_last = psc_output_fields_item_create_mfields(item_v);
   mflds_t_last = psc_output_fields_item_create_mfields(item_t);
   mflds_n_in = psc_output_fields_item_create_mfields(item_n);
  }

  int nr_kinds = ppsc->nr_kinds;

  psc_output_fields_item_run(item_n, ppsc->flds, mprts, mflds_n);

  psc_output_fields_item_run(item_v, ppsc->flds, mprts, mflds_v);

  psc_output_fields_item_run(item_t, ppsc->flds, mprts, mflds_t);
  psc_bnd_fill_ghosts(flds_bnd, mflds_t, 0, nr_kinds);

  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_n = psc_mfields_get_patch(mflds_n, p);

    // fix up zero density cells
    for (int m = 0; m < nr_kinds; m++) {
      for (int iz = 0; iz < ppatch->ldims[2]; iz++) {
	for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	  if (F3_C(flds_n, m, 0,iy,iz) == 0.0) {
	    F3_C(flds_n, m, 0,iy,iz) = 0.00001;
	  }
	}
      }
    }
  }    
  psc_bnd_fill_ghosts(flds_bnd, mflds_n, 0, nr_kinds);

  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_n = psc_mfields_get_patch(mflds_n, p);
    struct psc_fields *flds_v = psc_mfields_get_patch(mflds_v, p);
    
    // normalized v moments
    for (int m = 0; m < nr_kinds; m++) {
      for (int mm = 0; mm < 3; mm++) {
	for (int iz = 0; iz < ppatch->ldims[2]; iz++) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	    F3_C(flds_v, 3*m + mm, 0,iy,iz) /= F3_C(flds_n, m, 0,iy,iz);
	  }
	}
      }
    }
  }
  psc_bnd_fill_ghosts(flds_bnd, mflds_v, 0, 3 * nr_kinds);

  average_9_point(mflds_n_av, mflds_n);
  average_9_point(mflds_v_av, mflds_v);

  static struct mrc_io *io;
  if (!io) {
    io = mrc_io_create(psc_comm(ppsc));
    mrc_io_set_type(io, "xdmf_collective");
    mrc_io_set_param_string(io, "basename", "open");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
  }

  if (ppsc->timestep % 1 == 0) {
    mrc_io_open(io, "w", ppsc->timestep, ppsc->timestep * ppsc->dt);
    struct mrc_fld *mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
    mrc_fld_set_name(mrc_fld, "v");
    mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
    mrc_fld_set_param_int(mrc_fld, "nr_comps", 3 * nr_kinds);
    mrc_fld_setup(mrc_fld);
    mrc_fld_set_comp_name(mrc_fld, 0, "vx_e");
    mrc_fld_set_comp_name(mrc_fld, 1, "vy_e");
    mrc_fld_set_comp_name(mrc_fld, 2, "vz_e");
    mrc_fld_set_comp_name(mrc_fld, 3, "vx_i");
    mrc_fld_set_comp_name(mrc_fld, 4, "vy_i");
    mrc_fld_set_comp_name(mrc_fld, 5, "vz_i");
    copy_to_mrc_fld(mrc_fld, mflds_v);
    mrc_fld_write(mrc_fld, io);
    mrc_fld_destroy(mrc_fld);

    mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
    mrc_fld_set_name(mrc_fld, "n_av");
    mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
    mrc_fld_set_param_int(mrc_fld, "nr_comps", nr_kinds);
    mrc_fld_setup(mrc_fld);
    mrc_fld_set_comp_name(mrc_fld, 0, "n_av_e");
    mrc_fld_set_comp_name(mrc_fld, 1, "n_av_i");
    copy_to_mrc_fld(mrc_fld, mflds_n_av);
    mrc_fld_write(mrc_fld, io);
    mrc_fld_destroy(mrc_fld);

    mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
    mrc_fld_set_name(mrc_fld, "v_av");
    mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
    mrc_fld_set_param_int(mrc_fld, "nr_comps", 3 * nr_kinds);
    mrc_fld_setup(mrc_fld);
    mrc_fld_set_comp_name(mrc_fld, 0, "vx_av_e");
    mrc_fld_set_comp_name(mrc_fld, 1, "vy_av_e");
    mrc_fld_set_comp_name(mrc_fld, 2, "vz_av_e");
    mrc_fld_set_comp_name(mrc_fld, 3, "vx_av_i");
    mrc_fld_set_comp_name(mrc_fld, 4, "vy_av_i");
    mrc_fld_set_comp_name(mrc_fld, 5, "vz_av_i");
    copy_to_mrc_fld(mrc_fld, mflds_v_av);
    mrc_fld_write(mrc_fld, io);
    mrc_fld_destroy(mrc_fld);
  }

  average_in_time(mflds_n_av, mflds_n_last, first_time);
  average_in_time(mflds_v_av, mflds_v_last, first_time);

  if (ppsc->timestep % 1 == 0) {
    struct mrc_fld *mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
    mrc_fld_set_name(mrc_fld, "n_avt");
    mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
    mrc_fld_set_param_int(mrc_fld, "nr_comps", nr_kinds);
    mrc_fld_setup(mrc_fld);
    mrc_fld_set_comp_name(mrc_fld, 0, "n_avt_e");
    mrc_fld_set_comp_name(mrc_fld, 1, "n_avt_i");
    copy_to_mrc_fld(mrc_fld, mflds_n_av);
    mrc_fld_write(mrc_fld, io);
    mrc_fld_destroy(mrc_fld);

    mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
    mrc_fld_set_name(mrc_fld, "v_avt");
    mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
    mrc_fld_set_param_int(mrc_fld, "nr_comps", 3 * nr_kinds);
    mrc_fld_setup(mrc_fld);
    mrc_fld_set_comp_name(mrc_fld, 0, "vx_avt_e");
    mrc_fld_set_comp_name(mrc_fld, 1, "vy_avt_e");
    mrc_fld_set_comp_name(mrc_fld, 2, "vz_avt_e");
    mrc_fld_set_comp_name(mrc_fld, 3, "vx_avt_i");
    mrc_fld_set_comp_name(mrc_fld, 4, "vy_avt_i");
    mrc_fld_set_comp_name(mrc_fld, 5, "vz_avt_i");
    copy_to_mrc_fld(mrc_fld, mflds_v_av);
    mrc_fld_write(mrc_fld, io);
    mrc_fld_destroy(mrc_fld);

    mrc_io_close(io);
  }


  psc_mfields_destroy(mflds_v_av);
  psc_mfields_destroy(mflds_n_av);
  psc_mfields_destroy(mflds_t_av);
  psc_mfields_destroy(mflds_v);
  psc_mfields_destroy(mflds_n);
  psc_mfields_destroy(mflds_t);
  psc_mfields_destroy(mflds_ipr);
  psc_output_fields_item_destroy(item_n);
  psc_output_fields_item_destroy(item_v);
  psc_output_fields_item_destroy(item_t);
  psc_bnd_destroy(flds_bnd);

  first_time = false;
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_post

static void
psc_bnd_particles_sub_exchange_particles_post(struct psc_bnd_particles *bnd, struct psc_particles *prts)
{
  struct ddc_particles *ddcp = bnd->ddcp;
  struct ddcp_patch *patch = &ddcp->patches[prts->p];
  prts->n_part = patch->head;
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

static void
psc_bnd_particles_sub_exchange_particles(struct psc_bnd_particles *bnd, mparticles_base_t *particles_base)
{
  struct psc *psc = bnd->psc;

  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

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
    psc_bnd_particles_sub_exchange_particles_prep(bnd, psc_mparticles_get_patch(particles, p));
    psc_balance_comp_time_by_patch[p] += MPI_Wtime();
  }
  prof_stop(pr_time_step_no_comm);

  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);
  prof_stop(pr_B);

  prof_start(pr_C);
  for (int p = 0; p < particles->nr_patches; p++) {
    psc_bnd_particles_sub_exchange_particles_post(bnd, psc_mparticles_get_patch(particles, p));
  }
  prof_stop(pr_C);

  psc_bnd_particles_open_boundary(bnd, particles);

  psc_mparticles_put_cf(particles, particles_base, 0);
}

