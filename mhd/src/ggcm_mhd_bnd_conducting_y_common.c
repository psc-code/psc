#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <assert.h>

#include "ggcm_mhd_gkeyll.h"

static void
conducting_y_lo_mhd(struct mrc_fld *x, struct ggcm_mhd_bnd *bnd,
    int p, int m, int nx, int ny, int nz, int sw)
{
  // x lo
  // either the ig loop is the outer loop, or the iz/ix loops have to go backward
  for (int ig = 0; ig < sw; ig++) {      
#if 0
    int iyy = -ig;  // index of edge centered B_y's interior neighbors
#endif
    for (int iz = 0; iz < nz; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {

        M3(x, m+RR , ix, -1 - ig, iz, p) =   M3(x, m+RR,  ix, ig, iz, p);
        M3(x, m+RVX, ix, -1 - ig, iz, p) =   M3(x, m+RVX, ix, ig, iz, p);
        M3(x, m+RVY, ix, -1 - ig, iz, p) = - M3(x, m+RVY, ix, ig, iz, p);
        M3(x, m+RVZ, ix, -1 - ig, iz, p) =   M3(x, m+RVZ, ix, ig, iz, p);
        M3(x, m+UU , ix, -1 - ig, iz, p) =   M3(x, m+UU,  ix, ig, iz, p);
        M3(x, m+BX , ix, -1 - ig, iz, p) =   M3(x, m+BX,  ix, ig, iz, p);
        M3(x, m+BZ , ix, -1 - ig, iz, p) =   M3(x, m+BZ,  ix, ig, iz, p);

#if 0
        // to make div B = 0
        // for double / triple ghost points, we don't have the Bx/Bz we need
        // to make divB=0, but these By points are used by Ohm's Law, so let's
        // extrapolate
        if (ix + 1 == nx + sw || iz + 1 == nz + sw) {
          if (ix + 1 == nx + sw) {
            // extrapolate Bx value
            Bxp = 2.0f * M3(x, m+BX, ix, iyy, iz, p) - M3(x, m+BX, ix-1, iyy, iz, p);
          } else {
            // use the Bx value we have
            Bxp = M3(x, m+BX, ix+1, iyy, iz, p);
          }

          if (iz + 1 == nz + sw) {
            // extrapolate Bz value
            Bzp = 2 * M3(x, m+BZ, ix, iyy, iz, p) - M3(x, m+BZ, ix, iyy, iz-1);
          } else {
            // use the Bz value we have
            Bzp = M3(x, m+BZ, ix, iyy, iz+1);
          }

          M3(x, m+BY, ix, iyy, iz, p) = M3(x, m+BY, ix, iyy + 1, iz, p) + dx[1] * 
              ((Bxp - M3(x, m+BX, ix, iyy, iz, p)) / dx[0] +
               (Bzp - M3(x, m+BZ, ix, iyy, iz, p)) / dx[2]);
        } else {
          M3(x, m+BY, ix, iyy, iz, p) = M3(x, m+BY, ix, iyy + 1, iz, p) + dx[1] * 
              ((M3(x, m+BX, ix+1, iyy, iz  ) - M3(x, m+BX, ix, iyy, iz, p)) / dx[0] +
               (M3(x, m+BZ, ix  , iyy, iz+1) - M3(x, m+BZ, ix, iyy, iz, p)) / dx[2]);
        }
#endif
      }
    }
  }
}

static void
conducting_y_hi_mhd(struct mrc_fld *x, struct ggcm_mhd_bnd *bnd,
    int p, int m, int nx, int ny, int nz, int sw)
{
  // x hi
  for (int ig = 0; ig < sw; ig++) {      
#if 0
    int iyy = ny + ig - 1;  // index of edge centered B_y's interior neighbors
#endif
    for (int iz = 0; iz < nz; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {

        M3(x, m+RR , ix, ny + ig, iz, p) =   M3(x, m+RR , ix, ny - 1 - ig, iz, p);
        M3(x, m+RVX, ix, ny + ig, iz, p) =   M3(x, m+RVX, ix, ny - 1 - ig, iz, p);
        M3(x, m+RVY, ix, ny + ig, iz, p) = - M3(x, m+RVY, ix, ny - 1 - ig, iz, p);
        M3(x, m+RVZ, ix, ny + ig, iz, p) =   M3(x, m+RVZ, ix, ny - 1 - ig, iz, p);
        M3(x, m+UU , ix, ny + ig, iz, p) =   M3(x, m+UU , ix, ny - 1 - ig, iz, p);
        M3(x, m+BX , ix, ny + ig, iz, p) =   M3(x, m+BX , ix, ny - 1 - ig, iz, p);
        M3(x, m+BZ , ix, ny + ig, iz, p) =   M3(x, m+BZ , ix, ny - 1 - ig, iz, p);

        // to make div B = 0
        // for double / triple ghost points, we don't have the Bx/Bz we need
        // to make divB=0, but these By points are used by Ohm's Law, so let's
        // extrapolate
#if 0
        if (ix + 1 == nx + sw || iz + 1 == nz + sw) {
          if (ix + 1 == nx + sw) {
            // extrapolate Bx value
            Bxp = 2 * M3(x, m+BX, ix, iyy, iz, p) - M3(x, m+BX, ix-1, iyy, iz, p);
          } else {
            // use the Bx value we have
            Bxp = M3(x, m+BX, ix+1, iyy, iz, p);
          }

          if (iz + 1 == nz + sw) {
            // extrapolate Bz value
            Bzp = 2 * M3(x, m+BZ, ix, iyy, iz, p) - M3(x, m+BZ, ix, iyy, iz-1);
          } else {
            // use the Bz value we have
            Bzp = M3(x, m+BZ, ix, iyy, iz+1);
          }
          M3(x, m+BY, ix, iyy + 1, iz, p) = M3(x, m+BY, ix, iyy, iz, p) - dx[1] * 
              ((Bxp - M3(x, m+BX, ix, iyy, iz, p)) / dx[0] +
               (Bzp - M3(x, m+BZ, ix, iyy, iz, p)) / dx[2]);
        } else{
          M3(x, m+BY, ix, iyy + 1, iz, p) = M3(x, m+BY, ix, iyy, iz, p) - dx[1] * 
              ((M3(x, m+BX, ix+1, iyy, iz    ) - M3(x, m+BX, ix, iyy, iz, p)) / dx[0] +
               (M3(x, m+BZ, ix  , iyy, iz + 1) - M3(x, m+BZ, ix, iyy, iz, p)) / dx[2]);
        }
#endif
      }
    }
  }
}

static void
conducting_y_lo_gkeyll(struct mrc_fld *x, struct ggcm_mhd_bnd *bnd,
    int p, int m, int nx, int ny, int nz, int sw)
{
  struct ggcm_mhd *mhd = bnd->mhd;
  int nr_fluids = mhd->par.gk_nr_fluids;
  int nr_moments = mhd->par.gk_nr_moments;

  assert(nr_moments == 5);
  int idx[nr_fluids];
  ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);
  int idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);

  // x lo
  // either the ig loop is the outer loop, or the iz/ix loops have to go backward
  for (int ig = 0; ig < sw; ig++) {      
#if 0
    int iyy = -ig;  // index of edge centered B_y's interior neighbors
#endif
    for (int iz = 0; iz < nz; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {

        for (int s = 0; s < nr_fluids; s++) {
          M3(x, idx[s]+G5M_RRS , ix, -1 - ig, iz, p) =   M3(x, idx[s]+G5M_RRS,  ix, ig, iz, p);
          M3(x, idx[s]+G5M_RVXS, ix, -1 - ig, iz, p) =   M3(x, idx[s]+G5M_RVXS, ix, ig, iz, p);
          M3(x, idx[s]+G5M_RVYS, ix, -1 - ig, iz, p) = - M3(x, idx[s]+G5M_RVYS, ix, ig, iz, p);
          M3(x, idx[s]+G5M_RVZS, ix, -1 - ig, iz, p) =   M3(x, idx[s]+G5M_RVZS, ix, ig, iz, p);
          M3(x, idx[s]+G5M_UUS , ix, -1 - ig, iz, p) =   M3(x, idx[s]+G5M_UUS,  ix, ig, iz, p);
        }

        M3(x, idx_em+GK_EX , ix, -1 - ig, iz, p) = - M3(x, idx_em+GK_EX,  ix, ig, iz, p);
        M3(x, idx_em+GK_EY , ix, -1 - ig, iz, p) =   M3(x, idx_em+GK_EY,  ix, ig, iz, p);
        M3(x, idx_em+GK_EZ , ix, -1 - ig, iz, p) = - M3(x, idx_em+GK_EZ,  ix, ig, iz, p);

        M3(x, idx_em+GK_BX , ix, -1 - ig, iz, p) =   M3(x, idx_em+GK_BX,  ix, ig, iz, p);
        M3(x, idx_em+GK_BY , ix, -1 - ig, iz, p) = - M3(x, idx_em+GK_BY,  ix, ig, iz, p);
        M3(x, idx_em+GK_BZ , ix, -1 - ig, iz, p) =   M3(x, idx_em+GK_BZ,  ix, ig, iz, p);

        M3(x, idx_em+GK_PHI, ix, -1 - ig, iz, p) =   M3(x, idx_em+GK_PHI, ix, ig, iz, p);
        M3(x, idx_em+GK_PSI, ix, -1 - ig, iz, p) =   M3(x, idx_em+GK_PSI, ix, ig, iz, p);

      }
    }
  }
}

static void
conducting_y_hi_gkeyll(struct mrc_fld *x, struct ggcm_mhd_bnd *bnd,
    int p, int m, int nx, int ny, int nz, int sw)
{
  struct ggcm_mhd *mhd = bnd->mhd;
  int nr_fluids = mhd->par.gk_nr_fluids;
  int nr_moments = mhd->par.gk_nr_moments;

  assert(nr_moments == 5);
  int idx[nr_fluids];
  ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);
  int idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);

  // x hi
  for (int ig = 0; ig < sw; ig++) {      
#if 0
    int iyy = ny + ig - 1;  // index of edge centered B_y's interior neighbors
#endif
    for (int iz = 0; iz < nz; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {

        for (int s = 0; s < nr_fluids; s++) {
          M3(x, idx[s]+G5M_RRS , ix, ny + ig, iz, p) =   M3(x, idx[s]+G5M_RRS , ix, ny - 1 - ig, iz, p);
          M3(x, idx[s]+G5M_RVXS, ix, ny + ig, iz, p) =   M3(x, idx[s]+G5M_RVXS, ix, ny - 1 - ig, iz, p);
          M3(x, idx[s]+G5M_RVYS, ix, ny + ig, iz, p) = - M3(x, idx[s]+G5M_RVYS, ix, ny - 1 - ig, iz, p);
          M3(x, idx[s]+G5M_RVZS, ix, ny + ig, iz, p) =   M3(x, idx[s]+G5M_RVZS, ix, ny - 1 - ig, iz, p);
          M3(x, idx[s]+G5M_UUS , ix, ny + ig, iz, p) =   M3(x, idx[s]+G5M_UUS , ix, ny - 1 - ig, iz, p);
        }

        M3(x, idx_em+GK_EX , ix, ny + ig, iz, p) = - M3(x, idx_em+GK_EX , ix, ny - 1 - ig, iz, p);
        M3(x, idx_em+GK_EY , ix, ny + ig, iz, p) =   M3(x, idx_em+GK_EY , ix, ny - 1 - ig, iz, p);
        M3(x, idx_em+GK_EZ , ix, ny + ig, iz, p) = - M3(x, idx_em+GK_EZ , ix, ny - 1 - ig, iz, p);

        M3(x, idx_em+GK_BX , ix, ny + ig, iz, p) =   M3(x, idx_em+GK_BX , ix, ny - 1 - ig, iz, p);
        M3(x, idx_em+GK_BY , ix, ny + ig, iz, p) = - M3(x, idx_em+GK_BY , ix, ny - 1 - ig, iz, p);
        M3(x, idx_em+GK_BZ , ix, ny + ig, iz, p) =   M3(x, idx_em+GK_BZ , ix, ny - 1 - ig, iz, p);

        M3(x, idx_em+GK_PHI, ix, ny + ig, iz, p) =   M3(x, idx_em+GK_PHI, ix, ny - 1 - ig, iz, p);
        M3(x, idx_em+GK_PSI, ix, ny + ig, iz, p) =   M3(x, idx_em+GK_PSI, ix, ny - 1 - ig, iz, p);

      }
    }
  }
}

void
GGCM_MHD_BND_CONDUCTING_Y_FILL_GHOST(struct ggcm_mhd_bnd *bnd,
                          struct mrc_fld *fld_base,
                          int m, float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  struct mrc_fld *x = mrc_fld_get_as(fld_base, FLD_TYPE);
  int mhd_type;
  mrc_fld_get_param_int(x, "mhd_type", &mhd_type);
  assert(MT_BGRID(mhd_type) != MT_BGRID_FC_GGCM);  // written for C staggering only

  const int *dims = mrc_fld_spatial_dims(x);
  int nx = dims[0], ny = dims[1], nz = dims[2];
  int sw = x->_nr_ghosts;
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

  int bc[3];
  mrc_domain_get_param_int(mhd->domain, "bcx", &bc[0]); // FIXME in libmrc
  mrc_domain_get_param_int(mhd->domain, "bcy", &bc[1]);
  mrc_domain_get_param_int(mhd->domain, "bcz", &bc[2]);

  // struct mrc_patch_info info;
  // mrc_domain_get_local_patch_info(fld->_domain, 0, &info);

  /* float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2); */
  /* float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2); */
  /* float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2); */

  double dx[3];
  mrc_crds_get_dx_base(mrc_domain_get_crds(mhd->domain), dx);

#if 0
  // used for extrapolation for double / triple ghosts
  mrc_fld_data_t Bxp, Bzp;
#endif

  assert(mrc_fld_nr_patches(x) == 1);
  int p = 0;
  // lower boundary
  if (bc[1] != BC_PERIODIC && info.off[1] == 0) {
    if (mhd_type == MT_GKEYLL)
      conducting_y_lo_gkeyll(x, bnd, p, m, nx, ny, nz, sw);
    else
      conducting_y_lo_mhd(x, bnd, p, m, nx, ny, nz, sw);
  }

  // upper boundary
  if (bc[1] != BC_PERIODIC && info.off[1] + info.ldims[1] == gdims[1]) {
    if (mhd_type == MT_GKEYLL)
      conducting_y_hi_gkeyll(x, bnd, p, m, nx, ny, nz, sw);
    else
      conducting_y_hi_mhd(x, bnd, p, m, nx, ny, nz, sw);
  }

  mrc_fld_put_as(x, fld_base);
}
