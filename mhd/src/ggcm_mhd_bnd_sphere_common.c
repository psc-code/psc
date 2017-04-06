
#include "ggcm_mhd_bnd_private.h"

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>
#include <ggcm_mhd_gkeyll.h>

#include <mrc_domain.h>
#include <mrc_bits.h>
#include <math.h>

#include "pde/pde_defs.h"
#include "pde/pde_mhd_convert.c"

// FIXME, consolidate with ggcm_mhd_iono

// ======================================================================
// ggcm_mhd_bnd subclass "sphere"

struct ggcm_mhd_bnd_sphere {
  // params
  double radius;
  double bnvals[N_PRIMITIVE];  // constant values to set
  int test; // for testing, set to != 0
  int radial_velocity; // 0 : float, 1: reflect, 2: reflect if outflow
  double dr, extra_dr; // parameters for determining ghost points

  bool const_rr;
  bool const_pp;

  // state
  struct ggcm_mhd_bnd_sphere_map map;
};

#define ggcm_mhd_bnd_sphere(bnd) mrc_to_subobj(bnd, struct ggcm_mhd_bnd_sphere)

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_setup

static void
ggcm_mhd_bnd_sphere_setup(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;
  double xxnorm = bnd->mhd->xxnorm;

  pde_mhd_setup(bnd->mhd, mrc_fld_nr_comps(bnd->mhd->fld));

  ggcm_mhd_bnd_sphere_map_setup(map, bnd->mhd, sub->radius / xxnorm,
				sub->dr / xxnorm, sub->extra_dr / xxnorm);
  ggcm_mhd_bnd_sphere_map_setup_flds(map);
  ggcm_mhd_bnd_setup_member_objs_sub(bnd);
  ggcm_mhd_bnd_sphere_map_setup_cc(map);
  ggcm_mhd_bnd_sphere_map_setup_ec(map);
  ggcm_mhd_bnd_sphere_map_setup_fc(map);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_destroy

static void
ggcm_mhd_bnd_sphere_destroy(struct ggcm_mhd_bnd *bnd)
{
  pde_free();
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts

static void
sphere_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;
  struct ggcm_mhd *mhd = bnd->mhd;

  mrc_fld_data_t bnvals[N_PRIMITIVE];
  bnvals[RR] = sub->bnvals[RR] / mhd->rrnorm;
  bnvals[VX] = sub->bnvals[VX] / mhd->vvnorm;
  bnvals[VY] = sub->bnvals[VY] / mhd->vvnorm;
  bnvals[VZ] = sub->bnvals[VZ] / mhd->vvnorm;
  bnvals[PP] = sub->bnvals[PP] / mhd->ppnorm;
  bnvals[BX] = sub->bnvals[BX] / mhd->bbnorm;
  bnvals[BY] = sub->bnvals[BY] / mhd->bbnorm;
  bnvals[BZ] = sub->bnvals[BZ] / mhd->bbnorm;

  for (int i = 0; i < map->cc_n_map; i++) {
    int ix = MRC_I2(map->cc_imap, 0, i);
    int iy = MRC_I2(map->cc_imap, 1, i);
    int iz = MRC_I2(map->cc_imap, 2, i);
    int p  = MRC_I2(map->cc_imap, 3, i);
    
    // FIXME, this is still kinda specific / hacky to ganymede
    // to avoid cutting off the initial perturbation from e.g., the mirror dipole,
    // let's just keep B as-is, rather than using the fixed values above
    if (MT == MT_FCONS_CC) {
      bnvals[BX] = M3(fld, BX, ix,iy,iz, p);
      bnvals[BY] = M3(fld, BY, ix,iy,iz, p);
      bnvals[BZ] = M3(fld, BZ, ix,iy,iz, p);
    }
      
    mrc_fld_data_t state[s_n_state];
      convert_state_from_prim(state, bnvals);
      if (MT_BGRID(MT) == MT_BGRID_CC) {
	convert_put_state_to_3d(state, fld, ix,iy,iz, p);
      } else {
	// If B is not cell-centered, we don't set it at all, just the
	// fluid quantities. We'd need a face-centered map and things
	// are much more complicated, but fortunately, our boundary
	// condition is evolved using E = -v x B to update B anyway.
	convert_put_fluid_state_to_3d(state, fld, ix,iy,iz, p);
      }
  }
}

static void
sphere_fill_ghosts_test_1(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix - (1-bndp),iy,iz, p) = 0.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy - (1-bndp),iz, p) = 0.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz - (1-bndp), p) = 0.;
      }
    }
  }
  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix - (1-bndp),iy,iz, p) += 2.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy - (1-bndp),iz, p) += 2.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz - (1-bndp), p) += 2.;
      }
    }
  }
}

static void
sphere_fill_ghosts_test_2(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  struct mrc_fld *cnt = ggcm_mhd_get_3d_fld(bnd->mhd, 1);

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix-(1-bndp),iy,iz, p) = 0.;
	M3(cnt, 0, ix-(1-bndp),iy,iz, p) = 0.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy-(1-bndp),iz, p) = 0.;
	M3(cnt, 0, ix,iy-(1-bndp),iz, p) = 0.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz-(1-bndp), p) = 0.;
	M3(cnt, 0, ix,iy,iz-(1-bndp), p) = 0.;
      }
    }
  }
  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix-(1-bndp),iy,iz, p) += M3(fld, 0, ix-bndp,iy,iz, p);
	M3(cnt, 0, ix-(1-bndp),iy,iz, p) += 1.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy-(1-bndp),iz, p) += M3(fld, 0, ix,iy-bndp,iz, p);
	M3(cnt, 0, ix,iy-(1-bndp),iz, p) += 1.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz-(1-bndp), p) += M3(fld, 0, ix,iy,iz-bndp, p);
	M3(cnt, 0, ix,iy,iz-(1-bndp), p) += 1.;
      }
    }
  }

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix-(1-bndp),iy,iz, p) /= M3(cnt, 0, ix-(1-bndp),iy,iz, p);
	M3(cnt, 0, ix-(1-bndp),iy,iz, p) = 1.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy-(1-bndp),iz, p) /= M3(cnt, 0, ix,iy-(1-bndp),iz, p);
	M3(cnt, 0, ix,iy-(1-bndp),iz, p) = 1.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz-(1-bndp), p) /= M3(cnt, 0, ix,iy,iz-(1-bndp), p);
	M3(cnt, 0, ix,iy,iz-(1-bndp), p) = 1.;
      }
    }
  }
  ggcm_mhd_put_3d_fld(bnd->mhd, cnt);
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts_3

static void
sphere_fill_ghosts_test_3(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld)
{
#if MT == MT_GKEYLL
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;
  struct ggcm_mhd *mhd = bnd->mhd;

  mrc_fld_data_t bnvals[N_PRIMITIVE] = {};
  bnvals[RR] = sub->bnvals[RR] / mhd->rrnorm;
  bnvals[PP] = sub->bnvals[PP] / mhd->ppnorm;;

  mrc_fld_data_t state[s_n_state];
  convert_state_from_prim(state, bnvals);

  for (int i = 0; i < map->cc_n_map; i++) {
    int ix = MRC_I2(map->cc_imap, 0, i);
    int iy = MRC_I2(map->cc_imap, 1, i);
    int iz = MRC_I2(map->cc_imap, 2, i);
    int p  = MRC_I2(map->cc_imap, 3, i);

    convert_put_state_to_3d(state, fld, ix,iy,iz, p);
  }
#else
  assert(0);
#endif
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts_4
// use neighbor real cell values to float rho and pressure and B field
// and to radially reflect momentum

static void
sphere_fill_ghosts_test_4(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;
  struct ggcm_mhd *mhd = bnd->mhd;

  mrc_fld_data_t gamm = mhd->par.gamm;

  mrc_fld_data_t rrbn = sub->bnvals[RR] / mhd->rrnorm;
  mrc_fld_data_t ppbn = sub->bnvals[PP] / mhd->ppnorm;

  struct mrc_fld *ymask = mrc_fld_get_as(mhd->ymask, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double xc[3];

  if (MT == MT_GKEYLL) {
    int nr_fluids = mhd->par.gk_nr_fluids;
    int nr_moments = mhd->par.gk_nr_moments;
    int nr_comps = nr_fluids * nr_moments + 8;
    assert(nr_comps == fld->_nr_comps);
    
    float bn[nr_comps];

    float *mass = mhd->par.gk_mass.vals;
    float *pressure_ratios = mhd->par.gk_pressure_ratios.vals;

    int idx[nr_fluids];
    ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);
    int idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);

    float mass_ratios[nr_fluids];
    float mass_total = 0.;

    for (int s = 0; s < nr_fluids; s++)
      mass_total += mass[s];
    for (int s = 0; s < nr_fluids; s++)
      mass_ratios[s] = mass[s] / mass_total;

    for (int i = 0; i < map->cc_n_map; i++) {
      int ix = MRC_I2(map->cc_imap, 0, i);
      int iy = MRC_I2(map->cc_imap, 1, i);
      int iz = MRC_I2(map->cc_imap, 2, i);
      int p  = MRC_I2(map->cc_imap, 3, i);

      mrc_dcrds_at_cc(crds, ix,iy,iz, p, xc);
      mrc_fld_data_t x = xc[0];
      mrc_fld_data_t y = xc[1];
      mrc_fld_data_t z = xc[2];

      for (int c = 0; c < nr_comps; c++)
        bn[c] = 0.;
      int n_real = 0;

      for (int jx = ix-1; jx < ix+2; jx++) {
        for (int jy = iy-1; jy < iy+2; jy++) {
          for (int jz = iz-1; jz < iz+2; jz++) {
            if ((jx == ix && jy == iy && jz == iz)
            || jx < fld->_ghost_offs[1] || jx >= fld->_ghost_offs[1] + fld->_ghost_dims[1]
            || jy < fld->_ghost_offs[2] || jy >= fld->_ghost_offs[2] + fld->_ghost_dims[2]
            || jz < fld->_ghost_offs[3] || jz >= fld->_ghost_offs[3] + fld->_ghost_dims[3])
             continue;
            if (M3(ymask, 0, jx,jy,jz,p) > 0.) {
              for (int c = 0; c < nr_comps; c++) {
                bn[c] += M3(fld, c, jx,jy,jz, p);
              }
              n_real += 1;
            }
          }
        }
      }

      if (n_real == 0)
        continue;

      for (int c = 0; c < nr_comps; c++)
        bn[c] /= n_real;

      mrc_fld_data_t ir2 = 1/(x*x + y*y + z*z);
      for (int s = 0; s < nr_fluids; s++) {
        mrc_fld_data_t rr_s  = bn[idx[s] + G5M_RRS ];
        mrc_fld_data_t rvx_s = bn[idx[s] + G5M_RVXS];
        mrc_fld_data_t rvy_s = bn[idx[s] + G5M_RVYS];
        mrc_fld_data_t rvz_s = bn[idx[s] + G5M_RVZS];
        mrc_fld_data_t uu_s  = bn[idx[s] + G5M_UUS ];
        mrc_fld_data_t pp_s = (gamm-1.) * (uu_s 
            - .5 * (rvx_s*rvx_s + rvy_s*rvy_s + rvz_s*rvz_s) / rr_s);

        mrc_fld_data_t rr_s_ = sub->const_rr ? rrbn * mass_ratios[s] :rr_s;
        mrc_fld_data_t pp_s_ = sub->const_pp ? ppbn * pressure_ratios[s] :pp_s;
        mrc_fld_data_t rv_dot_r_s = rvx_s*x + rvy_s*y + rvz_s*z;
        mrc_fld_data_t rvx_s_ = rvx_s - 2*rv_dot_r_s*x*ir2;
        mrc_fld_data_t rvy_s_ = rvy_s - 2*rv_dot_r_s*y*ir2;
        mrc_fld_data_t rvz_s_ = rvz_s - 2*rv_dot_r_s*z*ir2;

        mrc_fld_data_t uu_s_ = pp_s_ / (gamm-1.) 
          + .5 * (rvx_s_*rvx_s_ + rvy_s_*rvy_s_ + rvz_s_*rvz_s_) / rr_s_;

        M3 (fld, idx[s] + G5M_RRS , ix,iy,iz, p) = rr_s_;
        M3 (fld, idx[s] + G5M_RVXS, ix,iy,iz, p) = rvx_s_;
        M3 (fld, idx[s] + G5M_RVYS, ix,iy,iz, p) = rvy_s_;
        M3 (fld, idx[s] + G5M_RVZS, ix,iy,iz, p) = rvz_s_;
        M3 (fld, idx[s] + G5M_UUS , ix,iy,iz, p) = uu_s_;
      }

      M3 (fld, idx_em+GK_EX, ix,iy,iz, p) = 0.;
      M3 (fld, idx_em+GK_EY, ix,iy,iz, p) = 0.;
      M3 (fld, idx_em+GK_EZ, ix,iy,iz, p) = 0.;

      M3 (fld, idx_em+GK_BX, ix,iy,iz, p) = bn[idx_em+GK_BX];
      M3 (fld, idx_em+GK_BY, ix,iy,iz, p) = bn[idx_em+GK_BY];
      M3 (fld, idx_em+GK_BZ, ix,iy,iz, p) = bn[idx_em+GK_BZ];
      
      M3 (fld, idx_em+GK_PHI, ix,iy,iz, p) = bn[idx_em+GK_PHI];
      M3 (fld, idx_em+GK_PSI, ix,iy,iz, p) = bn[idx_em+GK_PSI];
    }
  } else {
    int nr_comps = fld->_nr_comps;
    float bn[nr_comps];

    for (int i = 0; i < map->cc_n_map; i++) {
      int ix = MRC_I2(map->cc_imap, 0, i);
      int iy = MRC_I2(map->cc_imap, 1, i);
      int iz = MRC_I2(map->cc_imap, 2, i);
      int p  = MRC_I2(map->cc_imap, 3, i);

      mrc_dcrds_at_cc(crds, ix,iy,iz, p, xc);
      mrc_fld_data_t x = xc[0];
      mrc_fld_data_t y = xc[1];
      mrc_fld_data_t z = xc[2];

      for (int c = 0; c < nr_comps; c++)
        bn[c] = 0.;
      int n_real = 0;

      for (int jx = ix-1; jx < ix+2; jx++) {
        for (int jy = iy-1; jy < iy+2; jy++) {
          for (int jz = iz-1; jz < iz+2; jz++) {
            if ((jx == ix && jy == iy && jz == iz)
            || jx < fld->_ghost_offs[1] || jx >= fld->_ghost_offs[1] + fld->_ghost_dims[1]
            || jy < fld->_ghost_offs[2] || jy >= fld->_ghost_offs[2] + fld->_ghost_dims[2]
            || jz < fld->_ghost_offs[3] || jz >= fld->_ghost_offs[3] + fld->_ghost_dims[3])
             continue;
            if (M3(ymask, 0, jx,jy,jz,p) > 0.) {
              for (int c = 0; c < nr_comps; c++) {
                bn[c] += M3(fld, c, jx,jy,jz, p);
              }
              n_real += 1;
            }
          }
        }
      }

      if (n_real == 0)
        continue;

      for (int c = 0; c < nr_comps; c++)
        bn[c] /= n_real;

      mrc_fld_data_t rr = bn[RR];
      mrc_fld_data_t rr_ = sub->const_rr? rrbn : rr;

      mrc_fld_data_t rvx = rr * bn[VX];
      mrc_fld_data_t rvy = rr * bn[VY];
      mrc_fld_data_t rvz = rr * bn[VZ];
      
      mrc_fld_data_t ir2 = 1/(x*x + y*y + z*z);
      mrc_fld_data_t rv_dot_r = rvx*x + rvy*y + rvz*z;
      mrc_fld_data_t rvx_ = rvx - 2*rv_dot_r*x*ir2;
      mrc_fld_data_t rvy_ = rvy - 2*rv_dot_r*y*ir2;
      mrc_fld_data_t rvz_ = rvz - 2*rv_dot_r*z*ir2;

      M3(fld, RR,  ix,iy,iz, p) = rr_;
      M3(fld, RVX, ix,iy,iz, p) = rvx_;
      M3(fld, RVY, ix,iy,iz, p) = rvy_;
      M3(fld, RVZ, ix,iy,iz, p) = rvz_;

      mrc_fld_data_t bx = bn[BX];
      mrc_fld_data_t by = bn[BY];
      mrc_fld_data_t bz = bn[BZ];
      
      mrc_fld_data_t bx_ = bx;
      mrc_fld_data_t by_ = by;
      mrc_fld_data_t bz_ = bz;

      if (MT_FORMULATION(MT) == MT_FORMULATION_SCONS) {
        mrc_fld_data_t pp = (gamm -1.) * (bn[UU]
            - .5 * ( sqr(rvx ) + sqr(rvy ) + sqr(rvz ) ) / rr);
        mrc_fld_data_t pp_ = sub->const_pp? pp : ppbn;

        M3(fld, UU, ix,iy,iz, p) = pp_ / (gamm - 1.)
          + .5 * ( sqr(rvx_) + sqr(rvy_) + sqr(rvz_) ) / rr_;

      } else if (MT_FORMULATION(MT) == MT_FORMULATION_FCONS) {
        mrc_fld_data_t pp = (gamm -1.) * (bn[EE]
            - .5 * ( sqr(bx  ) + sqr(by  ) + sqr(bz  )) // /mu0
            - .5 * ( sqr(rvx ) + sqr(rvy ) + sqr(rvz )) / rr);
        mrc_fld_data_t pp_ = sub->const_pp? pp : ppbn;

        M3(fld, EE, ix,iy,iz, p) = pp_ / (gamm - 1.)
          - .5 * ( sqr(bx_ ) + sqr(by_ ) + sqr(bz_ )) // /mu0
          + .5 * ( sqr(rvx_) + sqr(rvy_) + sqr(rvz_)) / rr_;
      } else {
        assert(0);
      }
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_fill_ghosts

static void
ggcm_mhd_bnd_sphere_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld_base,
				float bntim)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  if (map->cc_n_map == 0) {
    return;
  }

  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);
  assert(mhd_type == MT);

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);
  if (mhd_type == MT_GKEYLL)
    assert(fld->_aos && fld->_c_order);
  if (sub->test == 0) {
    sphere_fill_ghosts(bnd, fld);
  } else if (sub->test == 1) {
    sphere_fill_ghosts_test_1(bnd, fld);
  } else if (sub->test == 2) {
    sphere_fill_ghosts_test_2(bnd, fld);
  } else if (sub->test == 3) {
    sphere_fill_ghosts_test_3(bnd, fld);
  } else if (sub->test == 4) {
    sphere_fill_ghosts_test_4(bnd, fld);
  } else {
    assert(0);
  }
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts_E

static void
sphere_fill_ghosts_E(struct ggcm_mhd_bnd *bnd, struct mrc_fld *E)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->ec_n_map[d]; i++) {
      int ix = MRC_I2(map->ec_imap[d], 0, i);
      int iy = MRC_I2(map->ec_imap[d], 1, i);
      int iz = MRC_I2(map->ec_imap[d], 2, i);
      int p  = MRC_I2(map->ec_imap[d], 3, i);

      M3(E, d, ix,iy,iz, p) = 0.f;
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_fill_ghosts_E

static void _mrc_unused
ggcm_mhd_bnd_sphere_fill_ghosts_E(struct ggcm_mhd_bnd *bnd, struct mrc_fld *E_base)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  if (map->ec_n_map[0] + map->ec_n_map[1] + map->ec_n_map[2] == 0) {
    return;
  }

  struct mrc_fld *E = mrc_fld_get_as(E_base, FLD_TYPE);
  sphere_fill_ghosts_E(bnd, E);
  mrc_fld_put_as(E, E_base);
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts_reconstr

static void
sphere_fill_ghosts_reconstr(struct ggcm_mhd_bnd *bnd, struct mrc_fld *U_l[],
			    struct mrc_fld *U_r[], int p)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;
  struct mrc_crds *crds = mrc_domain_get_crds(bnd->mhd->domain);

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      // FIXME, we really should have maps per patch, or some other way
      // to not go through the entire list here
      if (p != MRC_I2(map->fc_imap[d], 3, i)) {
	continue;
      }

      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      float crd_fc[3];
      mrc_crds_at_fc(crds, ix,iy,iz, p, d, crd_fc);

      int n_comps = mrc_fld_nr_comps(U_l[0]);
      // find true (inside) face values
      mrc_fld_data_t U_ghost[n_comps], U_true[n_comps];
      for (int m = 0; m < n_comps; m++) {
	if (bndp) {
	  U_true[m] = M3(U_l[d], m, ix,iy,iz, p);
	} else {
	  U_true[m] = M3(U_r[d], m, ix,iy,iz, p);
	}
      }

      // from BATSRUS / ganymede b.c.
      // for inflow float everything 
      for (int m = 0; m < n_comps; m++) {
	U_ghost[m] = U_true[m];
      }

      if (sub->radial_velocity != 0) {
	// calculate r^2
	mrc_fld_data_t r2 = 0.;
	for (int m = 0; m < 3; m++) {
	  r2 += sqr(crd_fc[m]);
	}
	
	// for outflow reflect radial velocity: uG = u - 2*(u.r)*r/r^2
	mrc_fld_data_t UdotR = 0.;
	for (int m = 0; m < 3; m++) {
	  UdotR += U_true[RVX + m] * crd_fc[m];
	}
	if (sub->radial_velocity == 1 || UdotR > 0.) {
	  for (int m = 0; m < 3; m++) {
	    U_ghost[RVX + m] = U_true[RVX + m] - 2. * UdotR / r2 * crd_fc[m];
	  }
	}
      }

      // store ghost values back
      for (int m = 0; m < n_comps; m++) {
	if (bndp) {
	  M3(U_r[d], m, ix, iy,iz, p) = U_ghost[m];
	} else {
	  M3(U_l[d], m, ix, iy,iz, p) = U_ghost[m];
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_fill_ghosts_reconstr

static void
ggcm_mhd_bnd_sphere_fill_ghosts_reconstr(struct ggcm_mhd_bnd *bnd, struct mrc_fld *U_l_base[],
					 struct mrc_fld *U_r_base[], int p)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  if (map->fc_n_map[0] + map->fc_n_map[1] + map->fc_n_map[2] == 0) {
    return;
  }

  struct mrc_fld *U_l[3], *U_r[3];
  for (int d = 0; d < 3; d++) {
    U_l[d] = mrc_fld_get_as(U_l_base[d], FLD_TYPE);
    U_r[d] = mrc_fld_get_as(U_r_base[d], FLD_TYPE);
  }

  for (int p = 0; p < mrc_fld_nr_patches(U_l[0]); p++) {
    sphere_fill_ghosts_reconstr(bnd, U_l, U_r, p);
  }

  for (int d = 0; d < 3; d++) {
    mrc_fld_put_as(U_l[d], U_l_base[d]);
    mrc_fld_put_as(U_r[d], U_r_base[d]);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd "sphere" subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bnd_sphere, x)
static struct param ggcm_mhd_bnd_sphere_descr[] = {
  { "radius"          , VAR(radius)          , PARAM_DOUBLE(1.)          },
  { "rr"              , VAR(bnvals[RR])      , PARAM_DOUBLE(1.)          },
  { "pp"              , VAR(bnvals[PP])      , PARAM_DOUBLE(1.)          },
  { "vx"              , VAR(bnvals[VX])      , PARAM_DOUBLE(0.)          },
  { "vy"              , VAR(bnvals[VY])      , PARAM_DOUBLE(0.)          },
  { "vz"              , VAR(bnvals[VZ])      , PARAM_DOUBLE(0.)          },
  { "bx"              , VAR(bnvals[BX])      , PARAM_DOUBLE(0.)          },
  { "by"              , VAR(bnvals[BY])      , PARAM_DOUBLE(0.)          },
  { "bz"              , VAR(bnvals[BZ])      , PARAM_DOUBLE(0.)          },
  { "test"            , VAR(test)            , PARAM_INT(0)              },
  { "radial_velocity" , VAR(radial_velocity) , PARAM_INT(0)              },
  { "dr"              , VAR(dr)              , PARAM_DOUBLE(.01)         },
  { "extra_dr"        , VAR(extra_dr)        , PARAM_DOUBLE(0.)          },

  { "map_dr"          , VAR(map.dr)          , MRC_VAR_DOUBLE            },
  { "map_extra_dr"    , VAR(map.extra_dr)    , MRC_VAR_DOUBLE            },
  { "min_dr"          , VAR(map.min_dr)      , MRC_VAR_DOUBLE            },
  { "map_radius"      , VAR(map.radius)      , MRC_VAR_DOUBLE            },
  { "r1"              , VAR(map.r1)          , MRC_VAR_DOUBLE            },
  { "cc_n_map"        , VAR(map.cc_n_map)    , MRC_VAR_INT               },
  { "cc_imap"         , VAR(map.cc_imap)     , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_n_map[0]"     , VAR(map.ec_n_map[0]) , MRC_VAR_INT               },
  { "ec_n_map[1]"     , VAR(map.ec_n_map[1]) , MRC_VAR_INT               },
  { "ec_n_map[2]"     , VAR(map.ec_n_map[2]) , MRC_VAR_INT               },
  { "ec_imap[0]"      , VAR(map.ec_imap[0])  , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_imap[1]"      , VAR(map.ec_imap[1])  , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_imap[2]"      , VAR(map.ec_imap[2])  , MRC_VAR_OBJ(mrc_fld)      },
  { "fc_n_map[0]"     , VAR(map.fc_n_map[0]) , MRC_VAR_INT               },
  { "fc_n_map[1]"     , VAR(map.fc_n_map[1]) , MRC_VAR_INT               },
  { "fc_n_map[2]"     , VAR(map.fc_n_map[2]) , MRC_VAR_INT               },
  { "fc_imap[0]"      , VAR(map.fc_imap[0])  , MRC_VAR_OBJ(mrc_fld)      },
  { "fc_imap[1]"      , VAR(map.fc_imap[1])  , MRC_VAR_OBJ(mrc_fld)      },
  { "fc_imap[2]"      , VAR(map.fc_imap[2])  , MRC_VAR_OBJ(mrc_fld)      },
  { "const_rr"        , VAR(const_rr)        , PARAM_BOOL(true)          },
  { "const_pp"        , VAR(const_pp)        , PARAM_BOOL(false)         },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_bnd subclass "sphere"

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere = {
  .name             = ggcm_mhd_bnd_sub_name,
  .size             = sizeof(struct ggcm_mhd_bnd_sphere),
  .param_descr      = ggcm_mhd_bnd_sphere_descr,
  .setup            = ggcm_mhd_bnd_sphere_setup,
  .destroy          = ggcm_mhd_bnd_sphere_destroy,
  .fill_ghosts      = ggcm_mhd_bnd_sphere_fill_ghosts,
  //  .fill_ghosts_E    = ggcm_mhd_bnd_sphere_fill_ghosts_E,
  .fill_ghosts_reconstr = ggcm_mhd_bnd_sphere_fill_ghosts_reconstr,
};

