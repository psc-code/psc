
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"

#include <mrc_io.h>
#include <mrc_ddc.h>
#include <mrc_fld_as_double.h>

#include <assert.h>

// ======================================================================
// ggcm_mhd_ic class

// ----------------------------------------------------------------------
// ggcm_mhd_ic_run

void
ggcm_mhd_ic_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  assert(mhd);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);

  if (ops->init_b0) {
    mhd->b0 = ggcm_mhd_get_3d_fld(mhd, 3);
    ops->init_b0(ic, mhd->b0);
    // FIXME, this doesn't set B values in exterior ghost points
    mrc_ddc_fill_ghosts_fld(mrc_domain_get_ddc(mhd->domain), 0, 3, mhd->b0);
  }

  if (ops->vector_potential) {
    /* initialize magnetic field from vector potential */
    struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);

    struct mrc_fld *A = mrc_domain_fld_create(mhd->domain, 2, "Ax:Ay:Az");
    mrc_fld_set_type(A, FLD_TYPE);
    mrc_fld_setup(A);
    mrc_fld_view(A);
    
    int gdims[3], p1x, p1y;
    mrc_domain_get_global_dims(mhd->domain, gdims);
    p1x = (gdims[0] > 1);
    p1y = (gdims[1] > 1);

    for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
      struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
      double dx[3];
      mrc_crds_get_dx(crds, p, dx);
      
      /* initialize vector potential */
      mrc_fld_foreach(fld, ix,iy,iz, 1, 2) {
	for (int m = 0; m < 3; m++) {
	  // FIXME, want double precision crds natively here
	  float crd_ec[3];
	  ggcm_mhd_get_crds_ec(mhd, ix,iy,iz, p, m, crd_ec);
	  double dcrd_ec[3] = { crd_ec[0], crd_ec[1], crd_ec[2] };
	  
	  M3(A, m, ix,iy,iz, p) = ops->vector_potential(ic, m, dcrd_ec);
	}
      } mrc_fld_foreach_end;
      
      /* initialize face-centered fields */
      mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
	BX_(fld, ix,iy,iz, p) =  (M3(A, 2, ix    , iy+p1y, iz, p) - M3(A, 2, ix,iy,iz, p)) / dx[1];
	BY_(fld, ix,iy,iz, p) = -(M3(A, 2, ix+p1x, iy    , iz, p) - M3(A, 2, ix,iy,iz, p)) / dx[0];
      } mrc_fld_foreach_end;
    }
    
    mrc_fld_destroy(A);
    mrc_fld_put_as(fld, mhd->fld);
  }

  if (ops->primitive) {
    /* initialize density, velocity, pressure */
    struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);

    for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
      mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
	float crd_cc[3];
	ggcm_mhd_get_crds_cc(mhd, ix,iy,iz, p, crd_cc);
	double dcrd_cc[3] = { crd_cc[0], crd_cc[1], crd_cc[2] };
	
	mrc_fld_data_t prim[5];
	for (int m = 0; m < 5; m++) {
	  prim[m] = ops->primitive(ic, m, dcrd_cc);
	}
	
	RR_(fld, ix,iy,iz, p) = prim[RR];
	VX_(fld, ix,iy,iz, p) = prim[VX];
	VY_(fld, ix,iy,iz, p) = prim[VY];
	PP_(fld, ix,iy,iz, p) = prim[PP];
      } mrc_fld_foreach_end;    
    }

    mrc_fld_put_as(fld, mhd->fld);

    ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
  }

  if (ops->run) {
    ops->run(ic);
  }

  if (ops->init_b0) {
    if (!ggcm_mhd_step_supports_b0(mhd->step)) {
      // if the stepper doesn't support a separate b0, 
      // add b0 into b, destroy b0 again.
      struct mrc_fld *b0 = mrc_fld_get_as(mhd->b0, FLD_TYPE);
      struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
      
      // FIXME, could use some axpy kinda thing
      for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
	mrc_fld_foreach(fld, ix,iy,iz, 2, 2) {
	  for (int d = 0; d < 3; d++) {
	    M3(fld, BX+d, ix,iy,iz, p) += M3(b0, d, ix,iy,iz, p);
	  }
	} mrc_fld_foreach_end;
      }
      
      mrc_fld_put_as(b0, mhd->b0);
      mrc_fld_put_as(fld, mhd->fld);
      
      ggcm_mhd_put_3d_fld(mhd, mhd->b0);
      mhd->b0 = NULL;
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_init

static void
ggcm_mhd_ic_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mirdip_float_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mirdip_double_ops);
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
};

