#include <mrc_profile.h>
#include <mrc_bits.h>

// ----------------------------------------------------------------------
// zmaskn

static inline void
zmaskn_inl(struct ggcm_mhd *mhd, struct mrc_fld *zmask, int m_zmask,
	   struct mrc_fld *ymask, int m_ymask, struct mrc_fld *x, struct mrc_fld *b0)
{
  float va02i = 1.f / sqr(mhd->par.speedlimit / mhd->vvnorm);
  float eps   = 1e-15f;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  for (int p = 0; p < mrc_fld_nr_patches(zmask); p++) {
    mrc_fld_foreach(zmask, ix,iy,iz, 1, 1) {
      mrc_fld_data_t bb = (sqr(.5f * (BT(x, 0, ix,iy,iz, p) + BT(x, 0, ix+dx,iy,iz, p))) +
			   sqr(.5f * (BT(x, 1, ix,iy,iz, p) + BT(x, 1, ix,iy+dy,iz, p))) +
			   sqr(.5f * (BT(x, 2, ix,iy,iz, p) + BT(x, 2, ix,iy,iz+dz, p))));
      float rrm = fmaxf(eps, bb * va02i);
      M3(zmask, m_zmask, ix,iy,iz, p) = M3(ymask, m_ymask, ix,iy,iz, p) *
	fminf(1.f, RR_(x, ix,iy,iz, p) / rrm);
    } mrc_fld_foreach_end;
  }
}

static void _mrc_unused
zmaskn(struct ggcm_mhd *mhd, struct mrc_fld *zmask, int m_zmask,
       struct mrc_fld *ymask, int m_ymask, struct mrc_fld *x)
{
  if (mhd->b0) {
    zmaskn_inl(mhd, zmask, m_zmask, ymask, m_ymask, x, mhd->b0);
  } else {
    zmaskn_inl(mhd, zmask, m_zmask, ymask, m_ymask, x, NULL);
  }
}

// ----------------------------------------------------------------------
// badval_checks_sc
// 
// Check for bad values and call wrongful death if we see any,
// abort with values:
//  Err Code       Reason
//     3           NaNs in first 8 components of x
//     4           density smaller than 0.0
//     5           pressure smaller than 0.0

static void _mrc_unused
badval_checks_sc(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *prim)
{
  static int pr = 0;
  float *crdx, *crdy, *crdz;
  int local_has_badval = 0;
  int global_has_badval = 0;
  
  int max_comp = MIN(mrc_fld_nr_comps(x), 8);

  mrc_fld_data_t ppmin = 0.0;
  mrc_fld_data_t rrmin = 0.0;  // mhd->par.rrmin / mhd->rrnorm

  if (mhd->do_badval_checks) {
    if (!pr) {
      pr = prof_register("badval_checks", 0.0, 0, 0);
    }
    
    prof_start(pr);

    crdx = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX1);
    crdy = ggcm_mhd_crds_get_crd(mhd->crds, 1, FX1);
    crdz = ggcm_mhd_crds_get_crd(mhd->crds, 2, FX1);

    // Check for negative pressure if we have a valid prim
    if (prim) {
      for (int p = 0; p < mrc_fld_nr_patches(prim); p++) {
        mrc_fld_foreach(prim, i,j,k, 0, 0) {
          if (PP_(prim, i,j,k, p) < ppmin) {
            local_has_badval = 5;
            mprintf("pressure @ (x=%g y=%g z=%g) = %lg < %lg\n",
                    crdx[i], crdy[j], crdz[k], PP_(prim, i,j,k, p), ppmin);
          }
        } mrc_fld_foreach_end;
        if (local_has_badval) {
          break;
        }
      }
    }
  
    // check for NaNs and negative density
    if (!local_has_badval && x) {
      for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
        mrc_fld_foreach(x, i,j,k, 0, 0) {
          // Check for negative density
          if (RR_(x, i,j,k, p) < rrmin) {
            local_has_badval = 4;
            mprintf("density @ (x=%g y=%g z=%g) = %lg < %lg\n",
                    crdx[i], crdy[j], crdz[k], RR_(prim, i,j,k, p), rrmin);
          }
          
          // Check for NaN
          for (int comp=0; comp < max_comp; comp++) {
            if isnan(M3(x, comp, i,j,k, p)) {
              local_has_badval = 3;
              mprintf("NaN in field %d @ (x=%g y=%g z=%g)\n",
                      comp, crdx[i], crdy[j], crdz[k]);
            }
          }
        } mrc_fld_foreach_end;
        if (local_has_badval) {
          break;
        }
      }
    }
    MPI_Allreduce(&local_has_badval, &global_has_badval, 1, MPI_INT, MPI_MAX,
                  ggcm_mhd_comm(mhd));
    if (global_has_badval) {
      ggcm_mhd_wrongful_death(mhd, x, global_has_badval);
    }
    
    prof_stop(pr);
  }
}

// ----------------------------------------------------------------------
// enforce_rrmin_sc
//
// nudge rr and uu such that rr >= rrmin if needed

static void _mrc_unused
enforce_rrmin_sc(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  static int pr = 0;
  mrc_fld_data_t rrmin, s;
  mrc_fld_data_t rr, rrvv, pp, uu, new_rr, new_uu;
  float *crdx, *crdy, *crdz;
  
  if (!pr) {
    pr = prof_register("enforce_rrmin", 0, 0, 0);
  }
  
  prof_start(pr);
  crdx = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX1);
  crdy = ggcm_mhd_crds_get_crd(mhd->crds, 1, FX1);
  crdz = ggcm_mhd_crds_get_crd(mhd->crds, 2, FX1);
  
  rrmin = mhd->par.rrmin / mhd->rrnorm;
  s = mhd->par.gamm - 1.0;
  
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    mrc_fld_foreach(x, ix,iy,iz, 0, 0) {
      rr = RR_(x, ix, iy, iz, p);
      if (rr < rrmin) {
        // get pressure
        rrvv = (sqr(RVX_(x, ix, iy, iz, p)) +
                sqr(RVY_(x, ix, iy, iz, p)) +
                sqr(RVZ_(x, ix, iy, iz, p)));
        uu = UU_(x, ix, iy, iz, p);
        pp = s * (uu - 0.5 * rrvv / rr);
        
        // set new values using rrmin
        new_rr = rrmin;
        new_uu = (pp / s) + (0.5 * rrvv / rrmin);
        RR_(x, ix, iy, iz, p) = new_rr;
        UU_(x, ix, iy, iz, p) = new_uu;
        
        mprintf("!! Note: enforcing min density at (x=%g y=%g z=%g): "
                "rr %lg -> %lg, uu %lg -> %lg\n",
                crdx[ix], crdy[iy], crdz[iz], rr, new_rr, uu, new_uu);
      }
    } mrc_fld_foreach_end;
  }
  prof_stop(pr);
}
