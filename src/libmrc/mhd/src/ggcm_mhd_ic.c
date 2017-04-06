
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_gkeyll.h"

#include <mrc_io.h>
#include <mrc_ddc.h>
#include <mrc_fld_as_double.h>

#include "pde/pde_defs.h"
#include "pde/pde_mhd_convert.c"

#include <assert.h>

// ======================================================================
// ggcm_mhd_ic class

typedef double (*vector_potential_f)(struct ggcm_mhd_ic *ic, int m, double crd[3]);
typedef double (*primitive_f)(struct ggcm_mhd_ic *ic, int m, double crd[3]);

// ----------------------------------------------------------------------
// ggcm_mhd_ic_B_from_vector_potential_fc
//
// initialize face-centered B from edge-centered vector potential

static void
ggcm_mhd_ic_B_from_vector_potential_fc(struct ggcm_mhd_ic *ic, struct mrc_fld *b,
				       vector_potential_f vector_potential, int S)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double xx[3];

  /* initialize face-centered B from edge-centered vector potential */
  for (int p = 0; p < mrc_fld_nr_patches(b); p++) {
    mrc_fld_foreach(b, ix,iy,iz, 2, 2) {
      mrc_dcrds_at_ec(crds, ix  , iy  , iz  , p, 0, xx);
      mrc_fld_data_t Ax    = vector_potential(ic, 0, xx);
      double x_y0 = xx[1], x_z0 = xx[2];

      mrc_dcrds_at_ec(crds, ix  , iy+1, iz  , p, 0, xx);
      mrc_fld_data_t Ax_yp = vector_potential(ic, 0, xx);
      double x_yp = xx[1];
 
      mrc_dcrds_at_ec(crds, ix  , iy  , iz+1, p, 0, xx);
      mrc_fld_data_t Ax_zp = vector_potential(ic, 0, xx);
      double x_zp = xx[2];


      mrc_dcrds_at_ec(crds, ix  , iy  , iz  , p, 1, xx);
      mrc_fld_data_t Ay    = vector_potential(ic, 1, xx);
      double y_x0 = xx[0], y_z0 = xx[2];

      mrc_dcrds_at_ec(crds, ix  , iy  , iz+1, p, 1, xx);
      mrc_fld_data_t Ay_zp = vector_potential(ic, 1, xx);
      double y_zp = xx[2];

      mrc_dcrds_at_ec(crds, ix+1, iy  , iz  , p, 1, xx);
      mrc_fld_data_t Ay_xp = vector_potential(ic, 1, xx);
      double y_xp = xx[0];


      mrc_dcrds_at_ec(crds, ix  , iy  , iz  , p, 2, xx);
      mrc_fld_data_t Az    = vector_potential(ic, 2, xx);
      double z_x0 = xx[0], z_y0 = xx[1];

      mrc_dcrds_at_ec(crds, ix+1, iy  , iz  , p, 2, xx);
      mrc_fld_data_t Az_xp = vector_potential(ic, 2, xx);
      double z_xp = xx[0];

      mrc_dcrds_at_ec(crds, ix  , iy+1, iz  , p, 2, xx);
      mrc_fld_data_t Az_yp = vector_potential(ic, 2, xx);
      double z_yp = xx[1];

      if (ix > -2) {
	M3(b, 0, ix-S,iy,iz, p) += (Az_yp - Az) / (z_yp - z_y0) - (Ay_zp - Ay) / (y_zp - y_z0);
      }
      if (iy > -2) {
	M3(b, 1, ix,iy-S,iz, p) += (Ax_zp - Ax) / (x_zp - x_z0) - (Az_xp - Az) / (z_xp - z_x0);
      }
      if (iz > -2) {
	M3(b, 2, ix,iy,iz-S, p) += (Ay_xp - Ay) / (y_xp - y_x0) - (Ax_yp - Ax) / (x_yp - x_y0);
      }
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_B_from_vector_potential_cc
//
// initialize face-centered B from edge-centered vector potential

static void
ggcm_mhd_ic_B_from_vector_potential_cc(struct ggcm_mhd_ic *ic, struct mrc_fld *b,
				       vector_potential_f vector_potential)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double crd_xp[3], crd_xm[3], crd_yp[3], crd_ym[3], crd_zp[3], crd_zm[3];

  for (int p = 0; p < mrc_fld_nr_patches(b); p++) {
    /* initialize face-centered fields */
    mrc_fld_foreach(b, ix,iy,iz, 1, 1) {
      mrc_dcrds_at_cc(crds, ix+1, iy  , iz  , p, crd_xp);
      mrc_dcrds_at_cc(crds, ix-1, iy  , iz  , p, crd_xm);
      mrc_dcrds_at_cc(crds, ix  , iy+1, iz  , p, crd_yp);
      mrc_dcrds_at_cc(crds, ix  , iy-1, iz  , p, crd_ym);
      mrc_dcrds_at_cc(crds, ix  , iy  , iz+1, p, crd_zp);
      mrc_dcrds_at_cc(crds, ix  , iy  , iz-1, p, crd_zm);

      mrc_fld_data_t Ax_yp = vector_potential(ic, 0, crd_yp);
      mrc_fld_data_t Ax_ym = vector_potential(ic, 0, crd_ym);
      mrc_fld_data_t Ax_zp = vector_potential(ic, 0, crd_zp);
      mrc_fld_data_t Ax_zm = vector_potential(ic, 0, crd_zm);

      mrc_fld_data_t Ay_xp = vector_potential(ic, 1, crd_xp);
      mrc_fld_data_t Ay_xm = vector_potential(ic, 1, crd_xm);
      mrc_fld_data_t Ay_zp = vector_potential(ic, 1, crd_zp);
      mrc_fld_data_t Ay_zm = vector_potential(ic, 1, crd_zm);

      mrc_fld_data_t Az_xp = vector_potential(ic, 2, crd_xp);
      mrc_fld_data_t Az_xm = vector_potential(ic, 2, crd_xm);
      mrc_fld_data_t Az_yp = vector_potential(ic, 2, crd_yp);
      mrc_fld_data_t Az_ym = vector_potential(ic, 2, crd_ym);

      M3(b, 0, ix,iy,iz, p) += (Az_yp - Az_ym) / (crd_yp[1] - crd_ym[1]) - (Ay_zp - Ay_zm) / (crd_zp[2] - crd_zm[2]);
      M3(b, 1, ix,iy,iz, p) += (Ax_zp - Ax_zm) / (crd_zp[2] - crd_zm[2]) - (Az_xp - Az_xm) / (crd_xp[0] - crd_xm[0]);
      M3(b, 2, ix,iy,iz, p) += (Ay_xp - Ay_xm) / (crd_xp[0] - crd_xm[0]) - (Ax_yp - Ax_ym) / (crd_yp[1] - crd_ym[1]);
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_B_from_vector_potential

static void
ggcm_mhd_ic_B_from_vector_potential(struct ggcm_mhd_ic *ic, struct mrc_fld *b,
				    vector_potential_f vector_potential)
{
  struct ggcm_mhd *mhd = ic->mhd;
  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);

  if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
    ggcm_mhd_ic_B_from_vector_potential_fc(ic, b, vector_potential, 0);
  } else if (MT_BGRID(mhd_type) == MT_BGRID_FC_GGCM) {
    ggcm_mhd_ic_B_from_vector_potential_fc(ic, b, vector_potential, 1);
  } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
    ggcm_mhd_ic_B_from_vector_potential_cc(ic, b, vector_potential);
  } else {
    mprintf("mhd_type %d unhandled\n", mhd_type);
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_B_from_primitive_fc
//
// initialize face-centered B directly

static void
ggcm_mhd_ic_B_from_primitive_fc(struct ggcm_mhd_ic *ic, struct mrc_fld *b, primitive_f primitive)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  for (int p = 0; p < mrc_fld_nr_patches(b); p++) {
    mrc_fld_foreach(b, ix,iy,iz, 0, 0) {
      for (int m = 0; m < 3; m++) {
	double dcrd_fc[3];
	mrc_dcrds_at_fc(crds, ix,iy,iz, p, m, dcrd_fc);
	
	M3(b, m, ix,iy,iz, p) += primitive(ic, BX + m, dcrd_fc);
      }
    } mrc_fld_foreach_end;    
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_B_from_primitive_cc
//
// initialize cell-centered B directly

static void
ggcm_mhd_ic_B_from_primitive_cc(struct ggcm_mhd_ic *ic, struct mrc_fld *b, primitive_f primitive)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  for (int p = 0; p < mrc_fld_nr_patches(b); p++) {
    mrc_fld_foreach(b, ix,iy,iz, 0, 0) {
      double dcrd_cc[3];
      mrc_dcrds_at_cc(crds, ix,iy,iz, p, dcrd_cc);
	
      for (int m = 0; m < 3; m++) {
	M3(b, m, ix,iy,iz, p) += primitive(ic, BX + m, dcrd_cc);
      }
    } mrc_fld_foreach_end;    
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_B_from_primitive

static void
ggcm_mhd_ic_B_from_primitive(struct ggcm_mhd_ic *ic, struct mrc_fld *b, primitive_f primitive)
{
  struct ggcm_mhd *mhd = ic->mhd;
  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);
  
  if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
    ggcm_mhd_ic_B_from_primitive_fc(ic, b, primitive);
  } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
    ggcm_mhd_ic_B_from_primitive_cc(ic, b, primitive);
  } else {
    mprintf("mhd_type %d unhandled\n", mhd_type);
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydro_from_primitive_semi
// 
// init semi conservative hydro state

static void
ggcm_mhd_ic_hydro_from_primitive_semi(struct ggcm_mhd_ic *ic, struct mrc_fld *fld)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      double dcrd_cc[3];
      mrc_dcrds_at_cc(crds, ix,iy,iz, p, dcrd_cc);
      
      mrc_fld_data_t prim[8], state[8];
      for (int m = 0; m < 5; m++) {
	prim[m] = ops->primitive(ic, m, dcrd_cc);
      }
      convert_scons_from_prim(state, prim);
      for (int m = 0; m < 5; m++) {
	M3(fld, m, ix,iy,iz, p) = state[m];
      }
    } mrc_fld_foreach_end;    
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydro_from_primitive_fully
//
// init fully conservative hydro state assuming that B is face centered

static void
ggcm_mhd_ic_hydro_from_primitive_fully(struct ggcm_mhd_ic *ic, struct mrc_fld *fld)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);

  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      double dcrd_cc[3];
      mrc_dcrds_at_cc(crds, ix,iy,iz, p, dcrd_cc);
      
      mrc_fld_data_t prim[8], state[8];
      for (int m = 0; m < 5; m++) {
	prim[m] = ops->primitive(ic, m, dcrd_cc);
      }
      prim[BX] = .5f * (BX_(fld, ix,iy,iz, p) + BX_(fld, ix+dx,iy,iz, p));
      prim[BY] = .5f * (BY_(fld, ix,iy,iz, p) + BY_(fld, ix,iy+dy,iz, p));
      prim[BZ] = .5f * (BZ_(fld, ix,iy,iz, p) + BZ_(fld, ix,iy,iz+dz, p));

      convert_fcons_from_prim(state, prim);

      for (int m = 0; m < 5; m++) {
	M3(fld, m, ix,iy,iz, p) = state[m];
      }
    } mrc_fld_foreach_end;    
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydro_from_primitive_fully_cc
//
// init fully conservative hydro state assuming that B is cell centered

static void
ggcm_mhd_ic_hydro_from_primitive_fully_cc(struct ggcm_mhd_ic *ic, struct mrc_fld *fld)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      double dcrd_cc[3];
      mrc_dcrds_at_cc(crds, ix,iy,iz, p, dcrd_cc);
      
      mrc_fld_data_t prim[8], state[8];
      for (int m = 0; m < 5; m++) {
	prim[m] = ops->primitive(ic, m, dcrd_cc);
      }
      prim[BX] = BX_(fld, ix,iy,iz, p);
      prim[BY] = BY_(fld, ix,iy,iz, p);
      prim[BZ] = BZ_(fld, ix,iy,iz, p);
      
      convert_fcons_from_prim(state, prim);

      for (int m = 0; m < 5; m++) {
	M3(fld, m, ix,iy,iz, p) = state[m];
      }
    } mrc_fld_foreach_end;    
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydro_from_primitive_gkeyll
//
// init gkeyll fluid state (assuming that B is cell centered)

static void
ggcm_mhd_ic_hydro_from_primitive_gkeyll(struct ggcm_mhd_ic *ic, struct mrc_fld *fld)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);

  // FIXME, generalize / use macro / whatever
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  mrc_fld_data_t gamma_m1 = mhd->par.gamm - 1.;

  float *qq = mhd->par.gk_charge.vals;
  float *mm = mhd->par.gk_mass.vals;
  float *pressure_ratios = ggcm_mhd_gkeyll_pressure_ratios(mhd);

  int nr_moments = ggcm_mhd_gkeyll_nr_moments(mhd);
  int nr_fluids = ggcm_mhd_gkeyll_nr_fluids(mhd);

  int idx[nr_fluids];
  ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);
  int idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);

  assert(nr_moments == 5);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      double dcrd_cc[3];
      mrc_dcrds_at_cc(crds, ix,iy,iz, p, dcrd_cc);

      double crd_xp[3], crd_xm[3], crd_yp[3], crd_ym[3], crd_zp[3], crd_zm[3];
      mrc_dcrds_at_cc(crds, ix+1, iy  , iz  , p, crd_xp);
      mrc_dcrds_at_cc(crds, ix-1, iy  , iz  , p, crd_xm);
      mrc_dcrds_at_cc(crds, ix  , iy+1, iz  , p, crd_yp);
      mrc_dcrds_at_cc(crds, ix  , iy-1, iz  , p, crd_ym);
      mrc_dcrds_at_cc(crds, ix  , iy  , iz+1, p, crd_zp);
      mrc_dcrds_at_cc(crds, ix  , iy  , iz-1, p, crd_zm);
      
      mrc_fld_data_t prim[5];
      for (int m = 0; m < 5; m++) {
	prim[m] = ops->primitive(ic, m, dcrd_cc);
      }

      mrc_fld_data_t Bz_yp = M3(fld, idx_em + GK_BZ, ix,iy+dy,iz, p);
      mrc_fld_data_t Bz_ym = M3(fld, idx_em + GK_BZ, ix,iy-dy,iz, p);
      mrc_fld_data_t By_zp = M3(fld, idx_em + GK_BY, ix,iy,iz+dz, p);
      mrc_fld_data_t By_zm = M3(fld, idx_em + GK_BY, ix,iy,iz-dz, p);

      mrc_fld_data_t Bx_zp = M3(fld, idx_em + GK_BX, ix,iy,iz+dz, p);
      mrc_fld_data_t Bx_zm = M3(fld, idx_em + GK_BX, ix,iy,iz-dz, p);
      mrc_fld_data_t Bz_xp = M3(fld, idx_em + GK_BZ, ix+dx,iy,iz, p);
      mrc_fld_data_t Bz_xm = M3(fld, idx_em + GK_BZ, ix-dx,iy,iz, p);
      
      mrc_fld_data_t By_xp = M3(fld, idx_em + GK_BY, ix+dx,iy,iz, p);
      mrc_fld_data_t By_xm = M3(fld, idx_em + GK_BY, ix-dx,iy,iz, p);
      mrc_fld_data_t Bx_yp = M3(fld, idx_em + GK_BX, ix,iy+dy,iz, p);
      mrc_fld_data_t Bx_ym = M3(fld, idx_em + GK_BX, ix,iy-dy,iz, p);

      mrc_fld_data_t jx = (Bz_yp - Bz_ym) / (crd_yp[1] - crd_ym[1]) - (By_zp - By_zm) / (crd_zp[2] - crd_zm[2]);
      mrc_fld_data_t jy = (Bx_zp - Bx_zm) / (crd_zp[2] - crd_zm[2]) - (Bz_xp - Bz_xm) / (crd_xp[0] - crd_xm[0]);
      mrc_fld_data_t jz = (By_xp - By_xm) / (crd_xp[0] - crd_xm[0]) - (Bx_yp - Bx_ym) / (crd_yp[1] - crd_ym[1]);

      // fluid quantities for each species
      assert(nr_fluids == 2);
      // FIXME, for more than two species, needs more information, and more thinking
      // (e.g., distribute momentum to match diamagnetic drift, but
      // also need to make sure that net (one-fluid) flow matches what MHD says.
      for (int s = 0; s < nr_fluids; s++) {
	M3(fld, idx[s] + G5M_RRS,  ix,iy,iz, p) = prim[RR] * mm[s] / (mm[0] + mm[1]);
	M3(fld, idx[s] + G5M_RVXS, ix,iy,iz, p) = (jx - prim[RR] * prim[VX] * qq[1-s]/mm[1-s]) / (qq[s] / mm[s] - qq[1-s] / mm[1-s]);
	M3(fld, idx[s] + G5M_RVYS, ix,iy,iz, p) = (jy - prim[RR] * prim[VY] * qq[1-s]/mm[1-s]) / (qq[s] / mm[s] - qq[1-s] / mm[1-s]);
	M3(fld, idx[s] + G5M_RVZS, ix,iy,iz, p) = (jz - prim[RR] * prim[VZ] * qq[1-s]/mm[1-s]) / (qq[s] / mm[s] - qq[1-s] / mm[1-s]);
	M3(fld, idx[s] + G5M_UUS,  ix,iy,iz, p) = 
	  prim[PP] * pressure_ratios[s] / gamma_m1
	  + .5 * (sqr(M3(fld, idx[s] + G5M_RVXS, ix,iy,iz, p)) +
		  sqr(M3(fld, idx[s] + G5M_RVYS, ix,iy,iz, p)) +
		  sqr(M3(fld, idx[s] + G5M_RVZS, ix,iy,iz, p)))
	  / M3(fld, idx[s] + G5M_RRS, ix,iy,iz,p);
      }

      // E=-vxB, i.e., only convection E field
      // FIXME, should really also contain diamagnetic drift
      // incl electron pressure gradient
      M3(fld, idx_em + GK_EX, ix,iy,iz, p) =
	- prim[VY] * M3(fld, idx_em + GK_BZ, ix,iy,iz, p)
	+ prim[VZ] * M3(fld, idx_em + GK_BY, ix,iy,iz, p);
      M3(fld, idx_em + GK_EY, ix,iy,iz, p) =
	- prim[VZ] * M3(fld, idx_em + GK_BX, ix,iy,iz, p)
	+ prim[VX] * M3(fld, idx_em + GK_BZ, ix,iy,iz, p);
      M3(fld, idx_em + GK_EZ, ix,iy,iz, p) =
	- prim[VX] * M3(fld, idx_em + GK_BY, ix,iy,iz, p)
	+ prim[VY] * M3(fld, idx_em + GK_BX, ix,iy,iz, p);
      
      // FIXME sensible correction potentials?
      // e.g., copy ghost vals for inoutflow bnd
      M3(fld, idx_em + GK_PHI, ix,iy,iz, p) = .0;
      M3(fld, idx_em + GK_PSI, ix,iy,iz, p) = .0;
    } mrc_fld_foreach_end;    
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydro_from_primitive

static void
ggcm_mhd_ic_hydro_from_primitive(struct ggcm_mhd_ic *ic, struct mrc_fld *fld)
{
  struct ggcm_mhd *mhd = ic->mhd;
  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS) {
    ggcm_mhd_ic_hydro_from_primitive_semi(ic, fld);
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
      ggcm_mhd_ic_hydro_from_primitive_fully(ic, fld);
    } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
      ggcm_mhd_ic_hydro_from_primitive_fully_cc(ic, fld);
    } else {
      assert(0);
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    ggcm_mhd_ic_hydro_from_primitive_gkeyll(ic, fld);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_run

void
ggcm_mhd_ic_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  assert(mhd);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);

  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);

  // FIXME, this should probably go somewhere where it's reusable
  int idx_BX;
  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS ||
      MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    idx_BX = BX;
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    int idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);
    idx_BX = idx_em + GK_BX;
  } else {
    assert(0);
  }

  struct mrc_fld *b = mrc_fld_make_view(fld, idx_BX, idx_BX + 3);

  /* initialize background magnetic field */
  if (mhd->b0) {
    struct mrc_fld *b0 = mrc_fld_get_as(mhd->b0, FLD_TYPE);
    if (ops->vector_potential_bg) {
      ggcm_mhd_ic_B_from_vector_potential(ic, b0, ops->vector_potential_bg);
    } else if (ops->primitive_bg) {
      ggcm_mhd_ic_B_from_primitive(ic, b0, ops->primitive_bg);
    }
    mrc_fld_put_as(b0, mhd->b0);
  } else {
    if (ops->vector_potential_bg) {
      ggcm_mhd_ic_B_from_vector_potential(ic, b, ops->vector_potential_bg);
    } else if (ops->primitive_bg) {
      ggcm_mhd_ic_B_from_primitive(ic, b, ops->primitive_bg);
    }
  }

  /* initialize magnetic field */
  if (ops->vector_potential) {
    ggcm_mhd_ic_B_from_vector_potential(ic, b, ops->vector_potential);
  } else if (ops->primitive) {
    ggcm_mhd_ic_B_from_primitive(ic, b, ops->primitive);
  }

  mrc_fld_destroy(b);

  /* initialize density, velocity, pressure, or corresponding
     conservative quantities */
  if (ops->primitive) {
    ggcm_mhd_ic_hydro_from_primitive(ic, fld);
  }

  mrc_fld_put_as(fld, mhd->fld);

  if (ops->run) {
    ops->run(ic);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_setup

static void
_ggcm_mhd_ic_setup(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_destroy

static void
_ggcm_mhd_ic_destroy(struct ggcm_mhd_ic *ic)
{
  pde_free();
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_init

static void
ggcm_mhd_ic_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mirdip_float_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mirdip_double_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_obstacle_double_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic, x)
static struct param ggcm_mhd_ic_descr[] = {
  { "mhd"             , VAR(mhd)             , PARAM_OBJ(ggcm_mhd)      },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic class description

struct mrc_class_ggcm_mhd_ic mrc_class_ggcm_mhd_ic = {
  .name             = "ggcm_mhd_ic",
  .size             = sizeof(struct ggcm_mhd_ic),
  .param_descr      = ggcm_mhd_ic_descr,
  .init             = ggcm_mhd_ic_init,
  .setup            = _ggcm_mhd_ic_setup,
  .destroy          = _ggcm_mhd_ic_destroy,
};

