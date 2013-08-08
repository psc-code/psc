
#ifndef MRC_CRDS_H
#define MRC_CRDS_H

#include <mrc_obj.h>
#include <mrc_fld.h>

struct mrc_crds {
  struct mrc_obj obj;
  // parameters
  float xl[3];
  float xh[3];
  int sw;

  // state
  struct mrc_domain *domain;
  struct mrc_fld *crd[3];
  struct mrc_fld *crd_nc[3];

  struct mrc_crds_gen *crds_gen[3];
};

#define MRC_CRD(crds, d, ix) MRC_F1((crds)->crd[d],0, ix)
#define MRC_CRDX(crds, ix) MRC_CRD(crds, 0, ix)
#define MRC_CRDY(crds, iy) MRC_CRD(crds, 1, iy)
#define MRC_CRDZ(crds, iz) MRC_CRD(crds, 2, iz)

#define MRC_MCRD(crds, d, ix, p) MRC_M1((crds)->crd[d],0, ix, p)
#define MRC_MCRDX(crds, ix, p) MRC_MCRD(crds, 0, ix, p)
#define MRC_MCRDY(crds, iy, p) MRC_MCRD(crds, 1, iy, p)
#define MRC_MCRDZ(crds, iz, p) MRC_MCRD(crds, 2, iz, p)

MRC_CLASS_DECLARE(mrc_crds, struct mrc_crds);

void mrc_crds_get_xl_xh(struct mrc_crds *crds, float xl[3], float xh[3]);
void mrc_crds_get_dx(struct mrc_crds *crds, float dx[3]);

struct mrc_crds_ops {
  MRC_SUBCLASS_OPS(struct mrc_crds);
};

#endif

