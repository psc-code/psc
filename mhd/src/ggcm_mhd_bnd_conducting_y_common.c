
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
  mrc_crds_get_dx(mrc_domain_get_crds(mhd->domain), dx);

  // lower boundary
  if (bc[1] != BC_PERIODIC && info.off[1] == 0) { // x lo
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
        for (int ig = 0; ig < sw; ig++) {      
          F3(x, m+_RR1 , ix, -1 - ig, iz) =   F3(x, m+_RR1,  ix, ig, iz);
          F3(x, m+_RV1X, ix, -1 - ig, iz) =   F3(x, m+_RV1X, ix, ig, iz);
          F3(x, m+_RV1Y, ix, -1 - ig, iz) = - F3(x, m+_RV1Y, ix, ig, iz);
          F3(x, m+_RV1Z, ix, -1 - ig, iz) =   F3(x, m+_RV1Z, ix, ig, iz);
          F3(x, m+UU   , ix, -1 - ig, iz) =   F3(x, m+UU  ,  ix, ig, iz);
          F3(x, m+BX   , ix, -1 - ig, iz) =   F3(x, m+BX  ,  ix, ig, iz);
          F3(x, m+BZ   , ix, -1 - ig, iz) =   F3(x, m+BZ  ,  ix, ig, iz);

          // // nothing special, B != 0, but Bz stays 0
          // F3(x, m+BY , ix, -1 - ig, iz) =   F3(x, m+BY,  ix, ig, iz);

          // to make div B = 0, but bz != 0 and it creeps in over time
          if (ix + 1 == nx + sw || iz + 1 == nz + sw) {
            // F3(x, m+BY, ix, iy + 1, iz) = 1.0 / bd2x[ix] + 1.0 / bd2z[iz];
            continue;
          }
          int iy = - ig;
          F3(x, m+BY, ix, iy, iz) = F3(x, m+BY, ix, iy + 1, iz) +
              ((F3(x, m+BX, ix+1, iy, iz    ) - F3(x, m+BX, ix, iy, iz)) / dx[0] +
               (F3(x, m+BZ, ix  , iy, iz + 1) - F3(x, m+BZ, ix, iy, iz)) / dx[2]) * dx[1];
        }
      }
    }
  }

  // upper boundary
  if (bc[1] != BC_PERIODIC && info.off[1] + info.ldims[1] == gdims[1]) { // x hi
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
        for (int ig = 0; ig < sw; ig++) {      
          F3(x, m+_RR1 , ix, ny + ig, iz) =   F3(x, m+_RR1 , ix, ny - 1 - ig, iz);
          F3(x, m+_RV1X, ix, ny + ig, iz) =   F3(x, m+_RV1X, ix, ny - 1 - ig, iz);
          F3(x, m+_RV1Y, ix, ny + ig, iz) = - F3(x, m+_RV1Y, ix, ny - 1 - ig, iz);
          F3(x, m+_RV1Z, ix, ny + ig, iz) =   F3(x, m+_RV1Z, ix, ny - 1 - ig, iz);
          F3(x, m+UU   , ix, ny + ig, iz) =   F3(x, m+UU   , ix, ny - 1 - ig, iz);
          F3(x, m+BX   , ix, ny + ig, iz) =   F3(x, m+BX   , ix, ny - 1 - ig, iz);
          F3(x, m+BZ   , ix, ny + ig, iz) =   F3(x, m+BZ   , ix, ny - 1 - ig, iz);

          // nothing special, B != 0, but Bz stays 0
          // F3(x, m+BY , ix, ny + ig, iz) =   F3(x, m+BY , ix, ny - 1 - ig, iz);

          // to make div B = 0, but bz != 0 and it creeps in over time
          if (ix + 1 == nx + sw || iz + 1 == nz + sw) {
            // F3(x, m+BY, ix, iy + 1, iz) = 1.0 / bd2x[ix] + 1.0 / bd2z[iz];
            continue;
          }
          int iy = ny - 1 + ig;
          /* F3(x, m+BY, ix, iy + 1, iz) = F3(x, m+BY, ix, iy, iz) - */
          /*     ((F3(x, m+BX, ix+1, iy, iz    ) - F3(x, m+BX, ix, iy, iz)) / bd2x[ix] + */
          /*      (F3(x, m+BZ, ix  , iy, iz + 1) - F3(x, m+BZ, ix, iy, iz)) / bd2z[iz]) * bd2y[iy]; */
          F3(x, m+BY, ix, iy + 1, iz) = F3(x, m+BY, ix, iy, iz) -
              ((F3(x, m+BX, ix+1, iy, iz    ) - F3(x, m+BX, ix, iy, iz)) / dx[0] +
               (F3(x, m+BZ, ix  , iy, iz + 1) - F3(x, m+BZ, ix, iy, iz)) / dx[2]) * dx[1];
        }
      }
    }
  }

  mrc_fld_put_as(x, fld_base);
}
