
#ifndef MRC_CRDS_H
#define MRC_CRDS_H

#include <mrc_obj.h>
#include <mrc_fld.h>

struct mrc_crds {
  struct mrc_obj obj;
  // parameters
  double xl[3];
  double xh[3];
  int sw;

  // state
  struct mrc_domain *domain;
  struct mrc_fld *crd[3];
  struct mrc_fld *dcrd[3]; // Double version of the coordinates
                           // not fully supported in io yet
  struct mrc_fld *crd_nc[3];
  struct mrc_fld *global_crd[3];

  struct mrc_crds_gen *crds_gen[3];
};

#define MRC_CRD(crds, d, ix) MRC_S2((crds)->crd[d], ix, 0)
#define MRC_CRDX(crds, ix) MRC_CRD(crds, 0, ix)
#define MRC_CRDY(crds, iy) MRC_CRD(crds, 1, iy)
#define MRC_CRDZ(crds, iz) MRC_CRD(crds, 2, iz)

#define MRC_DCRD(crds, d, ix) MRC_D2((crds)->dcrd[d], ix, 0)
#define MRC_DCRDX(crds, ix) MRC_DCRD(crds, 0, ix)
#define MRC_DCRDY(crds, iy) MRC_DCRD(crds, 1, iy)
#define MRC_DCRDZ(crds, iz) MRC_DCRD(crds, 2, iz)

#define MRC_MCRD(crds, d, ix, p) MRC_M1((crds)->crd[d],0, ix, p)
#define MRC_MCRDX(crds, ix, p) MRC_MCRD(crds, 0, ix, p)
#define MRC_MCRDY(crds, iy, p) MRC_MCRD(crds, 1, iy, p)
#define MRC_MCRDZ(crds, iz, p) MRC_MCRD(crds, 2, iz, p)

// Awkward macro for double precision multi-patch coords
#define MRC_DMCRD(crds, d, ix, p) MRC_D3((crds)->dcrd[d],ix, 0, p)
#define MRC_DMCRDX(crds, ix, p) MRC_DMCRD(crds, 0, ix, p)
#define MRC_DMCRDY(crds, iy, p) MRC_DMCRD(crds, 1, iy, p)
#define MRC_DMCRDZ(crds, iz, p) MRC_DMCRD(crds, 2, iz, p)

#define MRC_MCRD_NC(crds, d, ix, p) MRC_M1((crds)->crd_nc[d],0, ix, p)
#define MRC_MCRDX_NC(crds, ix, p) MRC_MCRD(crds, 0, ix, p)
#define MRC_MCRDY_NC(crds, iy, p) MRC_MCRD(crds, 1, iy, p)
#define MRC_MCRDZ_NC(crds, iz, p) MRC_MCRD(crds, 2, iz, p)


MRC_CLASS_DECLARE(mrc_crds, struct mrc_crds);

void mrc_crds_get_dx_base(struct mrc_crds *crds, double dx[3]);
void mrc_crds_get_dx(struct mrc_crds *crds, int p, double dx[3]);

struct mrc_crds_ops {
  MRC_SUBCLASS_OPS(struct mrc_crds);
};

#endif

