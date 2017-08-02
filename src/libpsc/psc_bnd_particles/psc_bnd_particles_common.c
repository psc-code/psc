
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_c.h"

#include <mrc_domain.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <string.h>

static const int debug_every_step = 10;

static inline double
random1()
{
  return random() / (double) RAND_MAX;
}

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
  const int b = 1;
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields *flds_av = psc_mfields_get_patch(mflds_av, p);
    
    if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int izb = 0, ize = ppatch->ldims[2];
	if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  izb++;
	}
	if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  ize--;
	}
	for (int iz = izb; iz < ize; iz++) {
	  F3_C(flds_av, m, 0,0,iz) = (F3_C(flds, m, 0,b+0,iz-1) +
				      F3_C(flds, m, 0,b+1,iz-1) +
				      F3_C(flds, m, 0,b+2,iz-1) +
				      F3_C(flds, m, 0,b+0,iz  ) +
				      F3_C(flds, m, 0,b+1,iz  ) +
				      F3_C(flds, m, 0,b+2,iz  ) +
				      F3_C(flds, m, 0,b+0,iz+1) +
				      F3_C(flds, m, 0,b+1,iz+1) +
				      F3_C(flds, m, 0,b+2,iz+1)) / 9.;
	  //F3_C(flds_av, m, 0,0,iz) = F3_C(flds, m, 0,b,iz);
	}
      }
    }

    if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int izb = 0, ize = ppatch->ldims[2];
	if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  izb++;
	}
	if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  ize--;
	}
	int iy = ppatch->ldims[1] - 1;
	for (int iz = izb; iz < ize; iz++) {
	  F3_C(flds_av, m, 0,iy,iz) = (F3_C(flds, m, 0,iy-b-0,iz-1) +
				       F3_C(flds, m, 0,iy-b-1,iz-1) +
				       F3_C(flds, m, 0,iy-b-2,iz-1) +
				       F3_C(flds, m, 0,iy-b-0,iz  ) +
				       F3_C(flds, m, 0,iy-b-1,iz  ) +
				       F3_C(flds, m, 0,iy-b-2,iz  ) +
				       F3_C(flds, m, 0,iy-b-0,iz+1) +
				       F3_C(flds, m, 0,iy-b-1,iz+1) +
				       F3_C(flds, m, 0,iy-b-2,iz+1)) / 9.;
	  //F3_C(flds_av, m, 0,iy,iz) = F3_C(flds, m, 0,iy-b,iz);
	}
      }
    }

    if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iyb = 0, iye = ppatch->ldims[1];
	if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iyb++;
	}
	if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iye--;
	}
	for (int iy = iyb; iy < iye; iy++) {
	  F3_C(flds_av, m, 0,iy,0) = (F3_C(flds, m, 0,iy-1,b+0) +
				      F3_C(flds, m, 0,iy-1,b+1) +
				      F3_C(flds, m, 0,iy-1,b+2) +
				      F3_C(flds, m, 0,iy  ,b+0) +
				      F3_C(flds, m, 0,iy  ,b+1) +
				      F3_C(flds, m, 0,iy  ,b+2) +
				      F3_C(flds, m, 0,iy+1,b+0) +
				      F3_C(flds, m, 0,iy+1,b+1) +
				      F3_C(flds, m, 0,iy+1,b+2)) / 9.;
	  //	  F3_C(flds_av, m, 0,iy,0) = F3_C(flds, m, 0,iy,b);
	}
      }
    }

    if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iz = ppatch->ldims[2] - 1;
	int iyb = 0, iye = ppatch->ldims[1];
	if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iyb++;
	}
	if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iye--;
	}
	for (int iy = iyb; iy < iye; iy++) {
	  F3_C(flds_av, m, 0,iy,iz) = (F3_C(flds, m, 0,iy-1,iz-b-0) +
				       F3_C(flds, m, 0,iy-1,iz-b-1) +
				       F3_C(flds, m, 0,iy-1,iz-b-2) +
				       F3_C(flds, m, 0,iy  ,iz-b-0) +
				       F3_C(flds, m, 0,iy  ,iz-b-1) +
				       F3_C(flds, m, 0,iy  ,iz-b-2) +
				       F3_C(flds, m, 0,iy+1,iz-b-0) +
				       F3_C(flds, m, 0,iy+1,iz-b-1) +
				       F3_C(flds, m, 0,iy+1,iz-b-2)) / 9.;
	  //	  F3_C(flds_av, m, 0,iy,iz) = F3_C(flds, m, 0,iy,iz-b);
	}
      }
    }

    if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN &&
	at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = 0;
	int iz = 0;
	F3_C(flds_av, m, 0,iy,iz) = (F3_C(flds, m, 0,iy+1,iz+1) +
				     F3_C(flds, m, 0,iy+1,iz+2) +
				     F3_C(flds, m, 0,iy+1,iz+3) +
				     F3_C(flds, m, 0,iy+2,iz+1) +
				     F3_C(flds, m, 0,iy+2,iz+2) +
				     F3_C(flds, m, 0,iy+2,iz+3) +
				     F3_C(flds, m, 0,iy+3,iz+1) +
				     F3_C(flds, m, 0,iy+3,iz+2) +
				     F3_C(flds, m, 0,iy+3,iz+3)) / 9.;
      }
    }

    if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN &&
	at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = 0;
	int iz = ppatch->ldims[2] - 1;
	F3_C(flds_av, m, 0,iy,iz) = (F3_C(flds, m, 0,iy+1,iz-1) +
				     F3_C(flds, m, 0,iy+1,iz-2) +
				     F3_C(flds, m, 0,iy+1,iz-3) +
				     F3_C(flds, m, 0,iy+2,iz-1) +
				     F3_C(flds, m, 0,iy+2,iz-2) +
				     F3_C(flds, m, 0,iy+2,iz-3) +
				     F3_C(flds, m, 0,iy+3,iz-1) +
				     F3_C(flds, m, 0,iy+3,iz-2) +
				     F3_C(flds, m, 0,iy+3,iz-3)) / 9.;
      }
    }

    if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN &&
	at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = ppatch->ldims[1] - 1;
	int iz = 0;
	F3_C(flds_av, m, 0,iy,iz) = (F3_C(flds, m, 0,iy-1,iz+1) +
				     F3_C(flds, m, 0,iy-1,iz+2) +
				     F3_C(flds, m, 0,iy-1,iz+3) +
				     F3_C(flds, m, 0,iy-2,iz+1) +
				     F3_C(flds, m, 0,iy-2,iz+2) +
				     F3_C(flds, m, 0,iy-2,iz+3) +
				     F3_C(flds, m, 0,iy-3,iz+1) +
				     F3_C(flds, m, 0,iy-3,iz+2) +
				     F3_C(flds, m, 0,iy-3,iz+3)) / 9.;
      }
    }

    if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN &&
	at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = ppatch->ldims[1] - 1;
	int iz = ppatch->ldims[2] - 1;
	F3_C(flds_av, m, 0,iy,iz) = (F3_C(flds, m, 0,iy-1,iz-1) +
				     F3_C(flds, m, 0,iy-1,iz-2) +
				     F3_C(flds, m, 0,iy-1,iz-3) +
				     F3_C(flds, m, 0,iy-2,iz-1) +
				     F3_C(flds, m, 0,iy-2,iz-2) +
				     F3_C(flds, m, 0,iy-2,iz-3) +
				     F3_C(flds, m, 0,iy-3,iz-1) +
				     F3_C(flds, m, 0,iy-3,iz-2) +
				     F3_C(flds, m, 0,iy-3,iz-3)) / 9.;
      }
    }
  }
}

static void
average_in_time(struct psc_bnd_particles *bnd,
		struct psc_mfields *mflds_av, struct psc_mfields *mflds_last)
{
  const double R = bnd->time_relax;
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_av = psc_mfields_get_patch(mflds_av, p);
    struct psc_fields *flds_last = psc_mfields_get_patch(mflds_last, p);
    
    for (int m = 0; m < mflds_av->nr_fields; m++) {
      if (bnd->first_time) {
	for (int iz = 0; iz < ppatch->ldims[2]; iz++) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	    F3_C(flds_last, m, 0,iy,iz) = F3_C(flds_av, m, 0,iy,iz);
	  }
	}
      } else {
	// at lower and upper z bnd only FIXME
	for (int iz = 0; iz < ppatch->ldims[2]; iz += ppatch->ldims[2] + 1) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	    F3_C(flds_av, m, 0,iy,iz) = 
	      R * F3_C(flds_av, m, 0,iy,iz)
	      + (1. - R) * F3_C(flds_last, m, 0,iy,iz);
	    F3_C(flds_last, m, 0,iy,iz) = F3_C(flds_av, m, 0,iy,iz);
	  }
	}
	for (int iz = 1; iz < ppatch->ldims[2] - 1; iz++) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy += ppatch->ldims[1] + 1) {
	    F3_C(flds_av, m, 0,iy,iz) = 
	      R * F3_C(flds_av, m, 0,iy,iz)
	      + (1. - R) * F3_C(flds_last, m, 0,iy,iz);
	    F3_C(flds_last, m, 0,iy,iz) = F3_C(flds_av, m, 0,iy,iz);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// ddcp_particles helpers

static void
ddcp_particles_realloc(void *_ctx, int p, int new_n_particles)
{
  mparticles_t *mprts = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  particles_realloc(prts, new_n_particles);
}

static void *
ddcp_particles_get_addr(void *_ctx, int p, int n)
{
  mparticles_t *mprts = _ctx;
  particle_range_t prts = particle_range_mprts(mprts, p);
  return particle_iter_at(prts.begin, n);
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

  bnd->item_nvt = psc_output_fields_item_create(psc_bnd_particles_comm(bnd));
  psc_output_fields_item_set_type(bnd->item_nvt, "nvp_1st_double");
  psc_output_fields_item_set_psc_bnd(bnd->item_nvt, bnd->flds_bnd);
  psc_output_fields_item_setup(bnd->item_nvt);

  bnd->mflds_nvt_av = psc_output_fields_item_create_mfields(bnd->item_nvt);
  bnd->mflds_nvt_last = psc_output_fields_item_create_mfields(bnd->item_nvt);
  // A slight FIXME: We're not tracking partial particles at corners separately per direction,
  // which could lead so systematic errors, though it seems unlikely.
  bnd->mflds_n_in = psc_mfields_create(psc_bnd_particles_comm(bnd));
  psc_mfields_set_type(bnd->mflds_n_in, "c");
  psc_mfields_set_domain(bnd->mflds_n_in, ppsc->mrc_domain);
  psc_mfields_set_param_int(bnd->mflds_n_in, "nr_fields", ppsc->nr_kinds);
  psc_mfields_set_param_int3(bnd->mflds_n_in, "ibn", ppsc->ibn);
  psc_mfields_setup(bnd->mflds_n_in);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
{
  ddc_particles_destroy(bnd->ddcp);

  psc_mfields_destroy(bnd->mflds_nvt_av);
  psc_mfields_destroy(bnd->mflds_nvt_last);
  psc_mfields_destroy(bnd->mflds_n_in);
  psc_output_fields_item_destroy(bnd->item_nvt);
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
psc_bnd_particles_sub_exchange_particles_prep(struct psc_bnd_particles *bnd,
					      struct psc_mparticles *mprts, int p)
{
  struct psc_particles *_prts = psc_mparticles_get_patch(mprts, p);
  particle_range_t prts = particle_range_mprts(mprts, p);
  struct ddc_particles *ddcp = bnd->ddcp;
  struct psc *psc = bnd->psc;

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.

  struct psc_patch *ppatch = &psc->patch[_prts->p];
  particle_real_t b_dxi[3] = { 1.f / ppatch->dx[0], 1.f / ppatch->dx[1], 1.f / ppatch->dx[2] };
  particle_real_t xm[3];
  int b_mx[3];
  for (int d = 0; d < 3; d++ ) {
    if (psc->domain.bnd_part_hi[d] == BND_PART_REFLECTING &&
	!psc->prm.gdims_in_terms_of_cells &&
	at_hi_boundary(_prts->p, d)) {
      b_mx[d] = ppatch->ldims[d] - 1;
    } else {
      b_mx[d] = ppatch->ldims[d];
    }
    xm[d] = b_mx[d] * ppatch->dx[d];
  }
  
  struct ddcp_patch *patch = &ddcp->patches[_prts->p];
  patch->head = 0;
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    patch->nei[dir1].n_send = 0;
  }
  for (int i = 0; i < _prts->n_part; i++) {
    particle_t *part = particle_iter_at(prts.begin, i);
    particle_real_t *xi = &part->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
    
    int b_pos[3];
    find_block_position(b_pos, xi, b_dxi);
    particle_real_t *pxi = &part->pxi;
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
	  if (!at_lo_boundary(_prts->p, d) ||
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
	  if (!at_hi_boundary(_prts->p, d) ||
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

static void _mrc_unused
psc_bnd_particles_sub_open_calc_moments(struct psc_bnd_particles *bnd,
					struct psc_mparticles *mprts_base)
{
  static int pr_A, pr_B, pr_C;
  if (!pr_A) {
    pr_A = prof_register("open_calc_moments", 1., 0, 0);
    pr_B = prof_register("open_item_run", 1., 0, 0);
    pr_C = prof_register("open_average", 1., 0, 0);
  }

  prof_start(pr_A);

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
  struct psc_mfields *mflds_nvt = psc_output_fields_item_create_mfields(bnd->item_nvt);

  prof_start(pr_B);
  psc_output_fields_item_run(bnd->item_nvt, ppsc->flds, mprts, mflds_nvt);
  prof_stop(pr_B);

  //debug_dump(io, mflds_nvt);

  prof_start(pr_C);
  average_9_point(bnd->mflds_nvt_av, mflds_nvt);
  debug_dump(io, bnd->mflds_nvt_av);

  average_in_time(bnd, bnd->mflds_nvt_av, bnd->mflds_nvt_last);
  prof_stop(pr_C);

  psc_mfields_destroy(mflds_nvt);

  bnd->first_time = false;

  psc_mparticles_put_cf(mprts, mprts_base, MP_DONT_COPY);

  if (ppsc->timestep % debug_every_step == 0) {
    mrc_io_close(io);
  }
  prof_stop(pr_A);
}

enum {
  XX,
  YY,
  ZZ,
  XY,
  YZ,
  ZX,
};

#ifndef NO_OPEN_BC

static void
calc_W(double W[6], double vv[6])
{
  double determ =
    (vv[XX]*vv[YY]*vv[ZZ] + vv[XY]*vv[YZ]*vv[ZX] + vv[ZX]*vv[XY]*vv[YZ] -
     vv[ZX]*vv[YY]*vv[ZX] - vv[XY]*vv[XY]*vv[ZZ] - vv[XX]*vv[YZ]*vv[YZ]);
  W[XX] = .5 * (vv[YY]*vv[ZZ] - vv[YZ]*vv[YZ]) / determ;
  W[YY] = .5 * (vv[XX]*vv[ZZ] - vv[ZX]*vv[ZX]) / determ;
  W[ZZ] = .5 * (vv[XX]*vv[YY] - vv[XY]*vv[XY]) / determ;
  W[XY] = .5 * (vv[XY]*vv[ZZ] - vv[ZX]*vv[YZ]) / determ;
  W[ZX] = .5 * (vv[XY]*vv[YZ] - vv[ZX]*vv[YY]) / determ;
  W[YZ] = .5 * (vv[XX]*vv[YZ] - vv[XY]*vv[ZX]) / determ;
}

enum {
  NVT_N,
  NVT_VX,
  NVT_VY,
  NVT_VZ,
  NVT_VXVX,
  NVT_VYVY,
  NVT_VZVZ,
  NVT_VXVY,
  NVT_VYVZ,
  NVT_VZVX,
};

static inline double
inject_particles(int p, struct psc_mparticles *mprts, struct psc_fields *flds, 
		 struct psc_fields *flds_nvt_av, int ix, int iy, int iz,
		 double ninjo, int kind, double pos[3], double dir,
		 int X, int Y, int Z)
{
  struct psc_particles *_prts = psc_mparticles_get_patch(mprts, p);
  particle_range_t prts = particle_range_mprts(mprts, p);

  double n     =         F3_C(flds_nvt_av, 10*kind + NVT_N     , ix,iy,iz);
  double v[3]  = {       F3_C(flds_nvt_av, 10*kind + NVT_VX, ix,iy,iz),
		         F3_C(flds_nvt_av, 10*kind + NVT_VY, ix,iy,iz),
			 F3_C(flds_nvt_av, 10*kind + NVT_VZ, ix,iy,iz), };
  double vv[6] = {       F3_C(flds_nvt_av, 10*kind + NVT_VXVX + X, ix,iy,iz),
			 F3_C(flds_nvt_av, 10*kind + NVT_VXVX + Y, ix,iy,iz),
			 F3_C(flds_nvt_av, 10*kind + NVT_VXVX + Z, ix,iy,iz),
			 F3_C(flds_nvt_av, 10*kind + NVT_VXVY + X, ix,iy,iz),
		   dir * F3_C(flds_nvt_av, 10*kind + NVT_VXVY + Y, ix,iy,iz),
		   dir * F3_C(flds_nvt_av, 10*kind + NVT_VXVY + Z, ix,iy,iz), };
  /* n = 1.; */
  /* v[0] = 0.; v[1] = 0.; v[2] = .1; */

  v[Z] *= dir;

  double W[6];
  calc_W(W, vv);

  double c=1.0;
  double vsz = sqrt(2. * vv[ZZ]);
  double gs0 = exp(-sqr(v[Z]) / sqr(vsz)) - exp(-sqr(c - v[Z]) / sqr(vsz))
    +sqrt(M_PI) * v[Z] / vsz * (erf((c - v[Z]) / vsz) + erf(v[Z] / vsz));
  double ninjn = ninjo + ppsc->dt * gs0 * n
    * vsz / sqrt(M_PI) / 2. / ppsc->patch[p].dx[Z] / ppsc->coeff.cori;

  int ninjc = (int) ninjn;
  /* mprintf("n ele %g ninjo %g ninjn %g ninjon %g ninjc %d\n", n, */
  /* 	  ninjo, ninjn, ninjn - ninjc, ninjc); */
  ninjo = ninjn - ninjc;

  int nvdx = 1000;
  double dvz = c / ((double) nvdx);
	  
  if (ninjc != 0) {
    double vzdin = 0.;
    double  fin[nvdx];
    for (int jj = 0; jj < nvdx; jj++){
      vzdin += dvz;
      fin[jj] = (exp(-sqr(v[Z] / vsz)) - exp(-sqr(vzdin - v[Z]) / sqr(vsz)) + 
		 sqrt(M_PI) * v[Z] / vsz * (erf((vzdin - v[Z]) / vsz) + erf(v[Z] / vsz))) / gs0;
    }
    for (int n = 0; n < ninjc; n++) {
      particle_t *prt = particle_iter_at(prts.begin, _prts->n_part++); 
      prt->kind = kind;
      prt->qni_wni = ppsc->kinds[kind].q;

      particle_real_t *pxi = &prt->pxi;
      particle_real_t *xi  = &prt->xi;

      int nnm = 0;
      do {
	nnm++;
	// FIXME, shouldn't have to loop here
	do {
	  double sr = random1();
      
	  pxi[Z] = 0.;
	  for (int k = 0; k < nvdx - 1; k++) {
	    if (sr >= fin[k] && sr < fin[k+1]) {
	      pxi[Z] = dvz * (k + 1) + (sr - fin[k]) * dvz / (fin[k+1] - fin[k]);
	      break;
	    }
	  }
	} while (pxi[Z] == 0);
        
	double sr = random1();
	double yya = 0.;
	double yy0;
	int icount = 0;
	do {
	  icount++;
	  yy0 = yya;
	  yya = yy0 - (erf(yy0) - (2.*sr - 1.)) / (2./sqrt(M_PI) * exp(-sqr(yy0)));
	} while (fabs(yya - yy0) > 1.e-15 && icount != 100);
	pxi[X] = v[X] + yya * sqrt(W[YY] / (W[XX] * W[YY] - sqr(W[XY])))
	  + (pxi[Z] - v[Z]) * vv[ZX] / vv[ZZ];
   
	sr = random1();
	yya = 0.0;
	icount = 0;
	do {
	  icount++;
	  yy0 = yya;
	  yya = yy0 - (erf(yy0) - (2.*sr - 1.)) / (2./sqrt(M_PI) * exp(-sqr(yy0)));
	} while (fabs(yya - yy0) > 1.e-15 && icount != 100);
	pxi[Y] = v[Y] + 1. / W[YY] * (yya * sqrt(W[YY])
				     - (pxi[Z] - v[Z]) * W[YZ] - (pxi[X] - v[X]) * W[XY]);
		
	if (nnm > 100) { assert(0); break; }
      } while (sqr(pxi[X]) + sqr(pxi[Y]) + sqr(pxi[Z]) > 1.);

      double xr = random1();
      xi[X] = pos[X] + xr * ppsc->patch[p].dx[X];
      double yr = random1();
      xi[Y] = pos[Y] + yr * ppsc->patch[p].dx[Y];
      double zr = random1();
      double dz = dir * zr * pxi[Z] * ppsc->dt;
      xi[Z] = pos[Z] + dz;

      double Jz = prt->qni_wni * dz * ppsc->coeff.cori / ppsc->dt;
      if (Z == 2) {
	F3_C(flds, JXI + Z, ix,iy  ,iz) += (1 - yr) * Jz;
	F3_C(flds, JXI + Z, ix,iy+1,iz) += (    yr) * Jz;
      } else if (Z == 1) {
	F3_C(flds, JXI + Z, ix,iy,iz  ) += (1 - xr) * Jz;
	F3_C(flds, JXI + Z, ix,iy,iz+1) += (    xr) * Jz;
      } else {
	assert(0);
      }
      double gamma = 1. / sqrt(1. - (sqr(pxi[X]) + sqr(pxi[Y]) + sqr(pxi[Z])));
      if (sqr(pxi[X]) + sqr(pxi[Y]) + sqr(pxi[Z]) > 1.) {
	gamma = 1.;
	assert(0);
      }
      pxi[X] *= gamma;
      pxi[Y] *= gamma;
      pxi[Z] *= dir * gamma;
    }
  }

  return ninjo;
}

static double
inject_particles_y(int p, struct psc_mparticles *mprts, struct psc_fields *flds, 
		   struct psc_fields *flds_nvt_av, int ix, int iy, int iz,
		   double ninjo, int kind, double pos[3], double dir)
{
  return inject_particles(p, mprts, flds, flds_nvt_av, ix, iy, iz, ninjo, kind, pos, dir,
			  2, 0, 1);
}

static double
inject_particles_z(int p, struct psc_mparticles *mprts, struct psc_fields *flds, 
		   struct psc_fields *flds_nvt_av, int ix, int iy, int iz,
		   double ninjo, int kind, double pos[3], double dir)
{
  return inject_particles(p, mprts, flds, flds_nvt_av, ix, iy, iz, ninjo, kind, pos, dir,
			  0, 1, 2);
}

#endif

static void _mrc_unused
psc_bnd_particles_open_boundary(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts,
				struct psc_mfields *mflds)
{
#ifndef NO_OPEN_BC
  static int pr_A;
  if (!pr_A) {
    pr_A = prof_register("open_boundary", 1., 0, 0);
  }

  prof_start(pr_A);
  int nr_kinds = ppsc->nr_kinds;

  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    struct psc_fields *flds_nvt_av = psc_mfields_get_patch(bnd->mflds_nvt_av, p);
    struct psc_fields *flds_n_in = psc_mfields_get_patch(bnd->mflds_n_in, p);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);

    for (int m = 0; m < nr_kinds; m++) {
      // inject at y lo
      if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN) {
	int iy = 0;
	for (int iz = 0; iz < ppatch->ldims[2]; iz++) {
	  double ninjo = F3_C(flds_n_in, m, 0,iy,iz);
	  double pos[3] = { 0., 0., iz * ppatch->dx[2], };
	  F3_C(flds_n_in, m, 0,iy,iz) =
	    inject_particles_y(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, +1.);
	}
      }
      // inject at y hi
      if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN) {
	int iy = ppatch->ldims[1] - 1;
	for (int iz = 0; iz < ppatch->ldims[2]; iz++) {
	  double ninjo = F3_C(flds_n_in, m, 0,iy,iz);
	  double pos[3] = { 0., (iy + 1) * (1-1e-6) * ppatch->dx[1], iz * ppatch->dx[2] };
	  F3_C(flds_n_in, m, 0,iy,iz) =
	    inject_particles_y(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, -1.);
	}
      }
      // inject at z lo
      if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
	int iz = 0;
	for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	  double ninjo = F3_C(flds_n_in, m, 0,iy,iz);
	  double pos[3] = { 0., iy * ppatch->dx[1], 0. };
	  F3_C(flds_n_in, m, 0,iy,iz) =
	    inject_particles_z(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, +1.);
	}
      }
      // inject at z hi
      if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
	int iz = ppatch->ldims[2] - 1;
	for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	  double ninjo = F3_C(flds_n_in, m, 0,iy,iz);
	  double pos[3] = { 0., iy * ppatch->dx[1], (iz + 1) * (1-1e-6) * ppatch->dx[2] };
	  F3_C(flds_n_in, m, 0,iy,iz) =
	    inject_particles_z(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, -1.);
	}
      }
    }
  }
  prof_stop(pr_A);
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
  //struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, "c", JXI, JXI + 3);

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
    psc_bnd_particles_sub_exchange_particles_post(bnd, psc_mparticles_get_patch(particles, p));
  }
  prof_stop(pr_C);

  //psc_bnd_particles_open_boundary(bnd, particles, mflds);

  psc_mparticles_put_cf(particles, particles_base, 0);
  //psc_mfields_put_as(mflds, psc->flds, JXI, JXI + 3);
}

