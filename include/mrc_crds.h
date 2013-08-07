
#ifndef MRC_CRDS_H
#define MRC_CRDS_H

#include <mrc_obj.h>
#include <mrc_fld.h>

struct mrc_crds_params {
  float xl[3];
  float xh[3];
  int sw;
};

struct mrc_crds {
  struct mrc_obj obj;
  struct mrc_crds_params par;
  struct mrc_domain *domain;
  struct mrc_fld *crd[3];
  struct mrc_m1 *mcrd[3];
  struct mrc_m1 *mcrd_nc[3];

  // temporally available between patch_get() and patch_put()
  struct mrc_fld_patch *mcrd_p[3];
};

#define MRC_CRD(crds, d, ix) MRC_F1((crds)->crd[d],0, ix)
#define MRC_CRDX(crds, ix) MRC_CRD(crds, 0, ix)
#define MRC_CRDY(crds, iy) MRC_CRD(crds, 1, iy)
#define MRC_CRDZ(crds, iz) MRC_CRD(crds, 2, iz)

#define MRC_MCRD(crds, d, ix) MRC_M1((crds)->mcrd_p[d],0, ix)
#define MRC_MCRDX(crds, ix) MRC_MCRD(crds, 0, ix)
#define MRC_MCRDY(crds, iy) MRC_MCRD(crds, 1, iy)
#define MRC_MCRDZ(crds, iz) MRC_MCRD(crds, 2, iz)

MRC_CLASS_DECLARE(mrc_crds, struct mrc_crds);

void mrc_crds_set_domain(struct mrc_crds *crds, struct mrc_domain *domain);
void mrc_crds_set_values(struct mrc_crds *crds, float *crdx, int mx,
			 float *crdy, int my, float *crdz, int mz);
void mrc_crds_get_xl_xh(struct mrc_crds *crds, float xl[3], float xh[3]);
void mrc_crds_get_dx(struct mrc_crds *crds, float dx[3]);
void mrc_crds_patch_get(struct mrc_crds *crds, int p);
void mrc_crds_patch_put(struct mrc_crds *crds);

struct mrc_crds_ops {
  MRC_SUBCLASS_OPS(struct mrc_crds);
  void (*set_values)(struct mrc_crds *crds, float *crdx, int mx,
		     float *crdy, int my, float *crdz, int mz);
};

#endif

