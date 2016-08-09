#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <assert.h>

void
GGCM_MHD_BND_CONDUCTING_Y_FILL_GHOST(struct ggcm_mhd_bnd *bnd,
                          struct mrc_fld *fld_base,
                          int m, float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  struct mrc_fld *x = mrc_fld_get_as(fld_base, FLD_TYPE);
  int mhd_type;
  mrc_fld_get_param_int(x, "mhd_type", &mhd_type);
  assert(mhd_type != MT_SEMI_CONSERVATIVE_GGCM);  // written for C staggering only

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

  // used for extrapolation for double / triple ghosts
  mrc_fld_data_t Bxp, Bzp;

  // lower boundary
  if (bc[1] != BC_PERIODIC && info.off[1] == 0) { // x lo
    // either the ig loop is the outer loop, or the iz/ix loops have to go backward
    for (int ig = 0; ig < sw; ig++) {      
      int iyy = -ig;  // index of edge centered B_y's interior neighbors
      for (int iz = 0; iz < nz; iz++) {
        for (int ix = -sw; ix < nx + sw; ix++) {

          F3(x, m+RR , ix, -1 - ig, iz) =   F3(x, m+RR,  ix, ig, iz);
          F3(x, m+RVX, ix, -1 - ig, iz) =   F3(x, m+RVX, ix, ig, iz);
          F3(x, m+RVY, ix, -1 - ig, iz) = - F3(x, m+RVY, ix, ig, iz);
          F3(x, m+RVZ, ix, -1 - ig, iz) =   F3(x, m+RVZ, ix, ig, iz);
          F3(x, m+UU , ix, -1 - ig, iz) =   F3(x, m+UU,  ix, ig, iz);
          F3(x, m+BX , ix, -1 - ig, iz) =   F3(x, m+BX,  ix, ig, iz);
          F3(x, m+BZ , ix, -1 - ig, iz) =   F3(x, m+BZ,  ix, ig, iz);

#if 0
          // to make div B = 0
          // for double / triple ghost points, we don't have the Bx/Bz we need
          // to make divB=0, but these By points are used by Ohm's Law, so let's
          // extrapolate
          if (ix + 1 == nx + sw || iz + 1 == nz + sw) {
            if (ix + 1 == nx + sw) {
              // extrapolate Bx value
              Bxp = 2.0f * F3(x, m+BX, ix, iyy, iz) - F3(x, m+BX, ix-1, iyy, iz);
            } else {
              // use the Bx value we have
              Bxp = F3(x, m+BX, ix+1, iyy, iz);
            }
  
            if (iz + 1 == nz + sw) {
              // extrapolate Bz value
              Bzp = 2 * F3(x, m+BZ, ix, iyy, iz) - F3(x, m+BZ, ix, iyy, iz-1);
            } else {
              // use the Bz value we have
              Bzp = F3(x, m+BZ, ix, iyy, iz+1);
            }

            F3(x, m+BY, ix, iyy, iz) = F3(x, m+BY, ix, iyy + 1, iz) + dx[1] * 
                ((Bxp - F3(x, m+BX, ix, iyy, iz)) / dx[0] +
                 (Bzp - F3(x, m+BZ, ix, iyy, iz)) / dx[2]);
          } else {
            F3(x, m+BY, ix, iyy, iz) = F3(x, m+BY, ix, iyy + 1, iz) + dx[1] * 
                ((F3(x, m+BX, ix+1, iyy, iz  ) - F3(x, m+BX, ix, iyy, iz)) / dx[0] +
                 (F3(x, m+BZ, ix  , iyy, iz+1) - F3(x, m+BZ, ix, iyy, iz)) / dx[2]);
          }
#endif
        }
      }
    }
  }

  // upper boundary
  if (bc[1] != BC_PERIODIC && info.off[1] + info.ldims[1] == gdims[1]) { // x hi
    for (int ig = 0; ig < sw; ig++) {      
      int iyy = ny + ig - 1;  // index of edge centered B_y's interior neighbors
      for (int iz = 0; iz < nz; iz++) {
        for (int ix = -sw; ix < nx + sw; ix++) {

          F3(x, m+RR , ix, ny + ig, iz) =   F3(x, m+RR , ix, ny - 1 - ig, iz);
          F3(x, m+RVX, ix, ny + ig, iz) =   F3(x, m+RVX, ix, ny - 1 - ig, iz);
          F3(x, m+RVY, ix, ny + ig, iz) = - F3(x, m+RVY, ix, ny - 1 - ig, iz);
          F3(x, m+RVZ, ix, ny + ig, iz) =   F3(x, m+RVZ, ix, ny - 1 - ig, iz);
          F3(x, m+UU , ix, ny + ig, iz) =   F3(x, m+UU , ix, ny - 1 - ig, iz);
          F3(x, m+BX , ix, ny + ig, iz) =   F3(x, m+BX , ix, ny - 1 - ig, iz);
          F3(x, m+BZ , ix, ny + ig, iz) =   F3(x, m+BZ , ix, ny - 1 - ig, iz);

          // to make div B = 0
          // for double / triple ghost points, we don't have the Bx/Bz we need
          // to make divB=0, but these By points are used by Ohm's Law, so let's
          // extrapolate
#if 0
          if (ix + 1 == nx + sw || iz + 1 == nz + sw) {
            if (ix + 1 == nx + sw) {
              // extrapolate Bx value
              Bxp = 2 * F3(x, m+BX, ix, iyy, iz) - F3(x, m+BX, ix-1, iyy, iz);
            } else {
              // use the Bx value we have
              Bxp = F3(x, m+BX, ix+1, iyy, iz);
            }
  
            if (iz + 1 == nz + sw) {
              // extrapolate Bz value
              Bzp = 2 * F3(x, m+BZ, ix, iyy, iz) - F3(x, m+BZ, ix, iyy, iz-1);
            } else {
              // use the Bz value we have
              Bzp = F3(x, m+BZ, ix, iyy, iz+1);
            }
            F3(x, m+BY, ix, iyy + 1, iz) = F3(x, m+BY, ix, iyy, iz) - dx[1] * 
                ((Bxp - F3(x, m+BX, ix, iyy, iz)) / dx[0] +
                 (Bzp - F3(x, m+BZ, ix, iyy, iz)) / dx[2]);
          } else{
            F3(x, m+BY, ix, iyy + 1, iz) = F3(x, m+BY, ix, iyy, iz) - dx[1] * 
                ((F3(x, m+BX, ix+1, iyy, iz    ) - F3(x, m+BX, ix, iyy, iz)) / dx[0] +
                 (F3(x, m+BZ, ix  , iyy, iz + 1) - F3(x, m+BZ, ix, iyy, iz)) / dx[2]);
          }
#endif
        }
      }
    }
  }

  mrc_fld_put_as(x, fld_base);
}
