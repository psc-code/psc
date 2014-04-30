
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_c.h"

#include <mrc_domain.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <string.h>

static const int debug_every_step = 1;

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
debug_dump(struct mrc_io *io, struct psc_mfields *mflds)
{
  if (ppsc->timestep % debug_every_step != 0) {
    return;
  }

  struct mrc_fld *mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
  mrc_fld_set_name(mrc_fld, psc_mfields_name(mflds));
  mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
  mrc_fld_set_param_int(mrc_fld, "nr_comps", mflds->nr_fields);
  mrc_fld_setup(mrc_fld);
  for (int m = 0; m < mflds->nr_fields; m++) {
    mrc_fld_set_comp_name(mrc_fld, m, psc_mfields_comp_name(mflds, m));
  }
  copy_to_mrc_fld(mrc_fld, mflds);
  mrc_fld_write(mrc_fld, io);
  mrc_fld_destroy(mrc_fld);
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
average_in_time(struct psc_bnd_particles *bnd,
		struct psc_mfields *mflds_av, struct psc_mfields *mflds_last)
{
  const double rr = 0.5; // FIXME
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_av = psc_mfields_get_patch(mflds_av, p);
    struct psc_fields *flds_last = psc_mfields_get_patch(mflds_last, p);
    
    for (int m = 0; m < mflds_av->nr_fields; m++) {
      if (!bnd->first_time) {
	for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	  F3_C(flds_av, m, 0,iy,0) = 
	    rr * F3_C(flds_av, m, 0,iy,0)
	    + (1. - rr) * F3_C(flds_last, m, 0,iy,0);
	}
      }

      for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
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

  // open b.c. setup (FIXME: only when necessary)

  bnd->first_time = true;

  // FIXME: not necessary (?)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srandom(rank);

  bnd->flds_bnd = psc_bnd_create(psc_bnd_particles_comm(bnd));
  psc_bnd_set_type(bnd->flds_bnd, "c");
  psc_bnd_set_psc(bnd->flds_bnd, ppsc);
  psc_bnd_setup(bnd->flds_bnd);

  bnd->item_n = psc_output_fields_item_create(psc_bnd_particles_comm(bnd));
  psc_output_fields_item_set_type(bnd->item_n, "n_1st_double");
  psc_output_fields_item_set_psc_bnd(bnd->item_n, bnd->flds_bnd);
  psc_output_fields_item_setup(bnd->item_n);

  bnd->item_v = psc_output_fields_item_create(psc_bnd_particles_comm(bnd));
  psc_output_fields_item_set_type(bnd->item_v, "v_1st_double");
  psc_output_fields_item_set_psc_bnd(bnd->item_v, bnd->flds_bnd);
  psc_output_fields_item_setup(bnd->item_v);

  bnd->item_t = psc_output_fields_item_create(psc_bnd_particles_comm(bnd));
  psc_output_fields_item_set_type(bnd->item_t, "T_1st_double");
  psc_output_fields_item_set_psc_bnd(bnd->item_t, bnd->flds_bnd);
  psc_output_fields_item_setup(bnd->item_t);

  bnd->mflds_n_av = psc_output_fields_item_create_mfields(bnd->item_n);
  bnd->mflds_v_av = psc_output_fields_item_create_mfields(bnd->item_v);
  bnd->mflds_t_av = psc_output_fields_item_create_mfields(bnd->item_t);
  bnd->mflds_n_last = psc_output_fields_item_create_mfields(bnd->item_n);
  bnd->mflds_v_last = psc_output_fields_item_create_mfields(bnd->item_v);
  bnd->mflds_t_last = psc_output_fields_item_create_mfields(bnd->item_t);
  bnd->mflds_n_in = psc_output_fields_item_create_mfields(bnd->item_n);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
{
  ddc_particles_destroy(bnd->ddcp);

  psc_mfields_destroy(bnd->mflds_v_av);
  psc_mfields_destroy(bnd->mflds_n_av);
  psc_mfields_destroy(bnd->mflds_t_av);
  psc_mfields_destroy(bnd->mflds_v_last);
  psc_mfields_destroy(bnd->mflds_n_last);
  psc_mfields_destroy(bnd->mflds_t_last);
  psc_output_fields_item_destroy(bnd->item_n);
  psc_output_fields_item_destroy(bnd->item_v);
  psc_output_fields_item_destroy(bnd->item_t);
  psc_bnd_destroy(bnd->flds_bnd);
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
psc_bnd_particles_sub_open_calc_moments(struct psc_bnd_particles *bnd,
					struct psc_mparticles *mprts_base)
{
  static struct mrc_io *io;
  if (!io) {
    io = mrc_io_create(psc_comm(ppsc));
    mrc_io_set_type(io, "xdmf_collective");
    mrc_io_set_param_string(io, "basename", "open");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
  }

  if (ppsc->timestep % debug_every_step == 0) {
    mrc_io_open(io, "w", ppsc->timestep, ppsc->timestep * ppsc->dt);
  }

  struct psc_mparticles *mprts = psc_mparticles_get_cf(mprts_base, 0);
  int nr_kinds = ppsc->nr_kinds;

  struct psc_mfields *mflds_n = psc_output_fields_item_create_mfields(bnd->item_n);
  struct psc_mfields *mflds_v = psc_output_fields_item_create_mfields(bnd->item_v);
  struct psc_mfields *mflds_t = psc_output_fields_item_create_mfields(bnd->item_t);

  psc_output_fields_item_run(bnd->item_n, ppsc->flds, mprts, mflds_n);
  psc_output_fields_item_run(bnd->item_v, ppsc->flds, mprts, mflds_v);
  psc_output_fields_item_run(bnd->item_t, ppsc->flds, mprts, mflds_t);

  debug_dump(io, mflds_n);

  // fix up zero density cells
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_n = psc_mfields_get_patch(mflds_n, p);

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
  psc_bnd_fill_ghosts(bnd->flds_bnd, mflds_n, 0, nr_kinds);

  // normalize v moments
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_n = psc_mfields_get_patch(mflds_n, p);
    struct psc_fields *flds_v = psc_mfields_get_patch(mflds_v, p);
    
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
  psc_bnd_fill_ghosts(bnd->flds_bnd, mflds_v, 0, 3 * nr_kinds);

  // calculate <(v - U)(v - U)> moments
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_n = psc_mfields_get_patch(mflds_n, p);
    struct psc_fields *flds_v = psc_mfields_get_patch(mflds_v, p);
    struct psc_fields *flds_t = psc_mfields_get_patch(mflds_t, p);
    
    const int mm2mx[6] = { 0, 1, 2, 0, 0, 1 };
    const int mm2my[6] = { 0, 1, 2, 1, 2, 2 };
    for (int m = 0; m < nr_kinds; m++) {
      for (int mm = 0; mm < 6; mm++) {
	int mx = mm2mx[mm], my = mm2my[mm];
	for (int iz = 0; iz < ppatch->ldims[2]; iz++) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	    F3_C(flds_t, 6*m + mm, 0,iy,iz) =
	      F3_C(flds_t, 6*m + mm, 0,iy,iz) / (ppsc->kinds[m].m * F3_C(flds_n, m, 0,iy,iz)) - 
	      F3_C(flds_v, 3*m + mx, 0,iy,iz) * F3_C(flds_v, 3*m + my, 0,iy,iz);
	  }
	}
      }
    }
  }
  psc_bnd_fill_ghosts(bnd->flds_bnd, mflds_t, 0, 6 * nr_kinds);

  debug_dump(io, mflds_n);
  debug_dump(io, mflds_v);
  debug_dump(io, mflds_t);

  average_9_point(bnd->mflds_n_av, mflds_n);
  average_9_point(bnd->mflds_v_av, mflds_v);
  average_9_point(bnd->mflds_t_av, mflds_t);

  average_in_time(bnd, bnd->mflds_n_av, bnd->mflds_n_last);
  average_in_time(bnd, bnd->mflds_v_av, bnd->mflds_v_last);
  average_in_time(bnd, bnd->mflds_t_av, bnd->mflds_t_last);

  psc_mfields_destroy(mflds_v);
  psc_mfields_destroy(mflds_n);
  psc_mfields_destroy(mflds_t);
  psc_mfields_destroy(mflds_nvt);

  bnd->first_time = false;

  psc_mparticles_put_cf(mprts, mprts_base, MP_DONT_COPY);

  if (ppsc->timestep % debug_every_step == 0) {
    mrc_io_close(io);
  }
}

static void
psc_bnd_particles_open_boundary(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
#ifndef NO_OPEN_BC
  int nr_kinds = ppsc->nr_kinds;

  struct psc_mfields *mflds_ipr = psc_output_fields_item_create_mfields(bnd->item_t);
  psc_mfields_set_comp_name(mflds_ipr, 0, "Wxxe");
  psc_mfields_set_comp_name(mflds_ipr, 1, "Wyye");
  psc_mfields_set_comp_name(mflds_ipr, 2, "Wzze");
  psc_mfields_set_comp_name(mflds_ipr, 3, "Wxye");
  psc_mfields_set_comp_name(mflds_ipr, 4, "Wxze");
  psc_mfields_set_comp_name(mflds_ipr, 5, "Wyze");
  psc_mfields_set_comp_name(mflds_ipr, 6, "Wxxi");
  psc_mfields_set_comp_name(mflds_ipr, 7, "Wyyi");
  psc_mfields_set_comp_name(mflds_ipr, 8, "Wzzi");
  psc_mfields_set_comp_name(mflds_ipr, 9, "Wxyi");
  psc_mfields_set_comp_name(mflds_ipr, 10, "Wxzi");
  psc_mfields_set_comp_name(mflds_ipr, 11, "Wyzi");

  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_ipr = psc_mfields_get_patch(mflds_ipr, p);
    struct psc_fields *flds_t_av = psc_mfields_get_patch(bnd->mflds_t_av, p);
    for (int m = 0; m < nr_kinds; m++) {
      for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	double determ = F3_C(flds_t_av, 6*m+0, 0,iy,0)*F3_C(flds_t_av, 6*m+1, 0,iy,0)*F3_C(flds_t_av, 6*m+2, 0,iy,0) 
	  +F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0)*F3_C(flds_t_av, 6*m+4, 0,iy,0)
	  +F3_C(flds_t_av, 6*m+4, 0,iy,0)*F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0)
	  -F3_C(flds_t_av, 6*m+4, 0,iy,0)*F3_C(flds_t_av, 6*m+1, 0,iy,0)*F3_C(flds_t_av, 6*m+4, 0,iy,0)
	  -F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+2, 0,iy,0)
	  -F3_C(flds_t_av, 6*m+0, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0);
	F3_C(flds_ipr, 6*m+0, 0,iy,0)=(F3_C(flds_t_av, 6*m+1, 0,iy,0)*F3_C(flds_t_av, 6*m+2, 0,iy,0)
				   -F3_C(flds_t_av, 6*m+5, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0))
	  /determ*0.5;
	F3_C(flds_ipr, 6*m+1, 0,iy,0)=(F3_C(flds_t_av, 6*m+0, 0,iy,0)*F3_C(flds_t_av, 6*m+2, 0,iy,0)
				   -F3_C(flds_t_av, 6*m+4, 0,iy,0)*F3_C(flds_t_av, 6*m+4, 0,iy,0))
	  /determ*0.5;
	F3_C(flds_ipr, 6*m+2, 0,iy,0)=(F3_C(flds_t_av, 6*m+0, 0,iy,0)*F3_C(flds_t_av, 6*m+1, 0,iy,0)
				   -F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+3, 0,iy,0))
	  /determ*0.5;
	F3_C(flds_ipr, 6*m+3, 0,iy,0)=(F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+2, 0,iy,0)
				   -F3_C(flds_t_av, 6*m+4, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0))
	  /determ*0.5;
	F3_C(flds_ipr, 6*m+4, 0,iy,0)=(F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0)
				   -F3_C(flds_t_av, 6*m+4, 0,iy,0)*F3_C(flds_t_av, 6*m+1, 0,iy,0))
	  /determ*0.5;
	F3_C(flds_ipr, 6*m+5, 0,iy,0)=(F3_C(flds_t_av, 6*m+0, 0,iy,0)*F3_C(flds_t_av, 6*m+5, 0,iy,0)
				   -F3_C(flds_t_av, 6*m+3, 0,iy,0)*F3_C(flds_t_av, 6*m+4, 0,iy,0))
	  /determ*0.5;
      }
    }
  }

  //  debug_dump(io, mflds_ipr);

  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_ipr = psc_mfields_get_patch(mflds_ipr, p);
    struct psc_fields *flds_n_av = psc_mfields_get_patch(bnd->mflds_n_av, p);
    struct psc_fields *flds_v_av = psc_mfields_get_patch(bnd->mflds_v_av, p);
    struct psc_fields *flds_t_av = psc_mfields_get_patch(bnd->mflds_t_av, p);
    struct psc_fields *flds_n_in = psc_mfields_get_patch(bnd->mflds_n_in, p);
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);

    // inject at z = 0
   if (ppatch->off[2]  == 0 ) {
    for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
     double vxea=F3_C(flds_v_av, 0, 0,iy,0);
     double vyea=F3_C(flds_v_av, 1, 0,iy,0);
     double vzea=F3_C(flds_v_av, 2, 0,iy,0);
     double c=1.0;
     int ninjo;
     int ninjn;
     int nvdx=100000;
     double dvz=c/((double) nvdx);

//     F3_C(flds_n_av, 0, 0,iy,0)=0.2;

     double  vsz=sqrt(2.0*(F3_C(flds_t_av, 2, 0,iy,0)));
     double  gs0=exp(-vzea*vzea/vsz/vsz)-exp(-(c-vzea)*(c-vzea)/vsz/vsz)
           +sqrt(M_PI)*vzea/vsz*(erf((c-vzea)/vsz)+erf(vzea/vsz)) ;
     ninjo=F3_C(flds_n_in, 0, 0,iy,0);
     F3_C(flds_n_in, 0, 0,iy,0)=F3_C(flds_n_in, 0, 0,iy,0)+ppsc->dt*gs0*F3_C(flds_n_av, 0, 0,iy,0)
           *vsz/sqrt(M_PI)/2.0/ppsc->patch[p].dx[2]/ppsc->coeff.cori;
     ninjn=F3_C(flds_n_in, 0, 0,iy,0);
     int ninjc=0;
     if(ninjo < ninjn){
      ninjc=ninjn-ninjo;
      //ninco=ninjn-ninjo;
     }
     //     if(ninjo > ninjn){ninco=0;}
/*     if(iy==ppatch->ldims[1]-1){ 
      mprintf("ninjn ele %d\n",ninjn);
      mprintf("ninjo ele %d\n",ninjo);
      mprintf("ninjc ele %d\n",ninjc);
      mprintf("vsz ele %f\n",vsz);
      mprintf("iy %d\n",iy);*/
      mprintf("iy ele %d\n",iy);
      mprintf("n ele %f\n",F3_C(flds_n_av, 0, 0,iy,0));
      mprintf("ninjc %d\n", ninjc);
/*      mprintf("nav ele %f\n",F3_C(flds_n_av, 0, 0,iy,0));
      mprintf("vxav ele %f\n",F3_C(flds_v_av, 0, 0,iy,0));
      mprintf("vyav ele %f\n",F3_C(flds_v_av, 1, 0,iy,0));
      mprintf("vzav ele %f\n",F3_C(flds_v_av, 2, 0,iy,0));
      mprintf("tav ele %f\n",F3_C(flds_t_av, 0, 0,iy,0));
     } */
     if(ninjc!=0) {
      double vzdin=0.0;
      double  fin[nvdx];
      for(int jj=0; jj < nvdx; jj++){
       vzdin=vzdin+dvz;
       fin[jj]=(exp(-vzea*vzea/vsz/vsz)-exp(-(vzdin-vzea)*
         (vzdin-vzea)/vsz/vsz)+sqrt(M_PI)*vzea/vsz*(erf((vzdin-vzea)/vsz)+ erf(vzea/vsz)))/gs0;
      }
      for(int n=0; n< ninjc; n++){
       /* if(!n){ */
       /*  srandom(seed); */
       /* } */
       particle_t *prt = particles_get_one(prts, prts->n_part++); 
       double sr;

          int nnm;
        nnm=0;
       do{
          nnm++;
       do{
	 long seed=random();
        sr=((double) seed)/((double)RAND_MAX);
      
        for(int k=0;k<nvdx-1;k++){
         prt->pzi=0.0;
         if(sr > fin[k] && sr < fin[k+1]){
          prt->pzi=dvz*((double) k+1.0)+(sr-fin[k])*dvz/(fin[k+1]-fin[k]);
	  prt->pzi = .1;
         }
         if(prt->pzi !=0.0) break;
        }
       } while(prt->pzi==0);
        
       long seed=random();
       sr=((double) seed)/((double)RAND_MAX);
       double yya=0.0;
       double yy0;
       int icount=0;
       do {
        icount++;
        yy0=yya;
        yya=yy0-(erf(yy0)-(2.0*sr-1.0))/(2.0/sqrt(M_PI)*exp(-yy0*yy0));
       } while(fabs(yya-yy0) > 1.0E-15 && icount!=100);
       prt->pxi=vxea+yya*sqrt(F3_C(flds_ipr, 1, 0,iy,0)/(F3_C(flds_ipr, 0, 0,iy,0)*F3_C(flds_ipr, 1, 0,iy,0)
          -F3_C(flds_ipr, 3, 0,iy,0)*F3_C(flds_ipr, 3, 0,iy,0)))
          +(prt->pzi-vzea)*F3_C(flds_t_av, 4, 0,iy,0)/F3_C(flds_t_av, 2, 0,iy,0);
   

       seed=random();
       sr=((double) seed)/((double)RAND_MAX);
       yya=0.0;
       icount=0;
       do {
        icount++;
        yy0=yya;
        yya=yy0-(erf(yy0)-(2.0*sr-1.0))/(2.0/sqrt(M_PI)*exp(-yy0*yy0));
       } while(fabs(yya-yy0) > 1.0E-15 && icount!=100);
       prt->pyi=vyea+1.0/F3_C(flds_ipr, 1, 0,iy,0)*(sqrt(F3_C(flds_ipr, 1, 0,iy,0))*yya
             -(prt->pzi-vzea)*F3_C(flds_ipr, 5, 0,iy,0)-(prt->pxi-vxea)*F3_C(flds_ipr, 3, 0,iy,0));

//      mprintf("p2 ele1 %f\n",prt->pxi*prt->pxi+prt->pyi*prt->pyi+prt->pzi*prt->pzi);
//      mprintf("nnm ele1 %d\n",nnm);

          if(prt->pxi*prt->pxi+prt->pyi*prt->pyi+prt->pzi*prt->pzi<1.0) break;
          if(nnm>100) break;
        } while(prt->pxi*prt->pxi+prt->pyi*prt->pyi+prt->pzi*prt->pzi>1.0);

       long seed=random();
       sr=((double) seed)/((double)RAND_MAX);
       sr = .5;
       prt->xi=sr*ppsc->patch[p].dx[0];
//       prt->yi=((double)ppatch->off[1]+(double)iy+sr)*ppsc->patch[p].dx[1];
       seed=random();
       sr=((double) seed)/((double)RAND_MAX);
       sr = .5;
       prt->yi=((double)iy+sr)*ppsc->patch[p].dx[1];
       seed=random();
       /* sr=((double) seed)/((double)RAND_MAX); */
//       prt->zi=sr*ppsc->patch[p].dx[2];
       prt->zi=0.;//0.0001*ppsc->patch[p].dx[2];
       prt->qni_wni=ppsc->kinds[0].q;
       //prt->mni=ppsc->kinds[0].m;
       prt->kind = KIND_ELECTRON;
/*      mprintf("xi ele1 %f\n",prt->xi);
      mprintf("yi ele1 %f\n",prt->yi);
      mprintf("zi ele2 %f\n",prt->zi);
      mprintf("pxi ele1 %f\n",prt->pxi);
      mprintf("pyi ele1 %f\n",prt->pyi);
      mprintf("pzi ele2 %f\n",prt->pzi);*/
       prt->pxi = 0.; prt->pyi = 0.;
       double gamma=1.0/sqrt(1.0-(prt->pxi*prt->pxi+prt->pyi*prt->pyi+prt->pzi*prt->pzi));
//      mprintf("gamma %f\n",gamma);
       if(prt->pxi*prt->pxi+prt->pyi*prt->pyi+prt->pzi*prt->pzi>1.0) gamma=1.0;
       prt->pxi=prt->pxi*gamma;
       prt->pyi=prt->pyi*gamma;
       prt->pzi=prt->pzi*gamma;
/*      mprintf("pxi ele z0%f\n",prt->pxi);
      mprintf("pyi ele z0%f\n",prt->pyi);
      mprintf("pzi ele z0%f\n",prt->pzi);*/
//       mprintf("mass %g\n",ppsc->kinds[1].q);
 
      }
     }
    }
   }
  }

  psc_mfields_destroy(mflds_ipr);
#endif
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

