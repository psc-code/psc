
#ifndef MRC_CRDS_H
#define MRC_CRDS_H

#include <mrc_obj.h>
#include <mrc_fld.h>

struct mrc_crds {
  struct mrc_obj obj;
  // parameters
  // domain limits in I/O units
  double l[3];
  double h[3];
  int sw;
  // normalization factors
  // typically, this might work like:
  // length input is in units of "cm".
  // then norm_length_scale transforms that back to the base SI units,
  // and norm_length transforms it into normalized units.
  // output will be scaled back into terms of "cm".
  double norm_length_scale;
  double norm_length;
  struct mrc_domain *domain;

  // state
  // domain limits in code units
  double xnorm;
  double lo_code[3];
  double hi_code[3];
  struct mrc_fld *crd[3];
  struct mrc_fld *dcrd[3]; // Double version of the coordinates
                           // not fully supported in io yet
  struct mrc_fld *crd_nc[3];
  struct mrc_fld *dcrd_nc[3];

  struct mrc_ndarray *global_crd[3];

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

#define MRC_MCRD_NC(crds, d, ix, p) MRC_M1((crds)->crd_nc[d], 0, ix, p)
#define MRC_MCRDX_NC(crds, ix, p) MRC_MCRD_NC(crds, 0, ix, p)
#define MRC_MCRDY_NC(crds, iy, p) MRC_MCRD_NC(crds, 1, iy, p)
#define MRC_MCRDZ_NC(crds, iz, p) MRC_MCRD_NC(crds, 2, iz, p)

#define MRC_DMCRD_NC(crds, d, ix, p) MRC_D3((crds)->dcrd_nc[d], ix, 0, p)
#define MRC_DMCRDX_NC(crds, ix, p) MRC_DMCRD_NC(crds, 0, ix, p)
#define MRC_DMCRDY_NC(crds, iy, p) MRC_DMCRD_NC(crds, 1, iy, p)
#define MRC_DMCRDZ_NC(crds, iz, p) MRC_DMCRD_NC(crds, 2, iz, p)


MRC_CLASS_DECLARE(mrc_crds, struct mrc_crds);

// get coordinate limits lo, hi (per dimension) in code units
const double *mrc_crds_lo(struct mrc_crds *crds);
const double *mrc_crds_hi(struct mrc_crds *crds);
void mrc_crds_get_dx_base(struct mrc_crds *crds, double dx[3]);
void mrc_crds_get_dx(struct mrc_crds *crds, int p, double dx[3]);

struct mrc_crds_ops {
  MRC_SUBCLASS_OPS(struct mrc_crds);
};

// ----------------------------------------------------------------------
// mrc_crds_at_cc

static inline void
mrc_crds_at_cc(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       float crd_cc[3])
{
  crd_cc[0] = MRC_MCRDX(crds, ix, p);
  crd_cc[1] = MRC_MCRDY(crds, iy, p);
  crd_cc[2] = MRC_MCRDZ(crds, iz, p);
}

// ----------------------------------------------------------------------
// mrc_crds_at_nc

static inline void
mrc_crds_at_nc(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       float crd_nc[3])
{
  crd_nc[0] = MRC_MCRDX_NC(crds, ix, p);
  crd_nc[1] = MRC_MCRDY_NC(crds, iy, p);
  crd_nc[2] = MRC_MCRDZ_NC(crds, iz, p);
}

// ----------------------------------------------------------------------
// mrc_crds_at_fc

static inline void
mrc_crds_at_fc(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       int d, float crd_fc[3])
{
  if (d == 0) {
    // Bx located at i, j+.5, k+.5
    crd_fc[0] = MRC_MCRDX_NC(crds, ix, p);
    crd_fc[1] = MRC_MCRDY   (crds, iy, p);
    crd_fc[2] = MRC_MCRDZ   (crds, iz, p);
  } else if (d == 1) {
    // By located at i+.5, j, k+.5
    crd_fc[0] = MRC_MCRDX   (crds, ix, p);
    crd_fc[1] = MRC_MCRDY_NC(crds, iy, p);
    crd_fc[2] = MRC_MCRDZ   (crds, iz, p);
  } else if (d == 2) {
    // Bz located at i+.5, j+.5, k
    crd_fc[0] = MRC_MCRDX   (crds, ix, p);
    crd_fc[1] = MRC_MCRDY   (crds, iy, p);
    crd_fc[2] = MRC_MCRDZ_NC(crds, iz, p);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_crds_at_ec

static inline void
mrc_crds_at_ec(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       int d, float crd_ec[3])
{
  if (d == 0) {
    // Ex located at i+.5, j, k
    crd_ec[0] = MRC_MCRDX   (crds, ix, p);
    crd_ec[1] = MRC_MCRDY_NC(crds, iy, p);
    crd_ec[2] = MRC_MCRDZ_NC(crds, iz, p);
  } else if (d == 1) {
    // Ey located at i, j+.5, k
    crd_ec[0] = MRC_MCRDX_NC(crds, ix, p);
    crd_ec[1] = MRC_MCRDY   (crds, iy, p);
    crd_ec[2] = MRC_MCRDZ_NC(crds, iz, p);
  } else if (d == 2) {
    // Ez located at i, j, k+.5
    crd_ec[0] = MRC_MCRDX_NC(crds, ix, p);
    crd_ec[1] = MRC_MCRDY_NC(crds, iy, p);
    crd_ec[2] = MRC_MCRDZ   (crds, iz, p);
  } else {
    assert(0);
  }
}

// ======================================================================
// and now everything repeated in double precision
//
// FIXME, this can probably be done more nicely using the "mrc_fld_as_double.h"
// mechanism...


// ----------------------------------------------------------------------
// mrc_dcrds_at_cc

static inline void
mrc_dcrds_at_cc(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       double crd_cc[3])
{
  crd_cc[0] = MRC_DMCRDX(crds, ix, p);
  crd_cc[1] = MRC_DMCRDY(crds, iy, p);
  crd_cc[2] = MRC_DMCRDZ(crds, iz, p);
}

// ----------------------------------------------------------------------
// mrc_dcrds_at_nc

static inline void
mrc_dcrds_at_nc(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       double crd_nc[3])
{
  crd_nc[0] = MRC_DMCRDX_NC(crds, ix, p);
  crd_nc[1] = MRC_DMCRDY_NC(crds, iy, p);
  crd_nc[2] = MRC_DMCRDZ_NC(crds, iz, p);
}

// ----------------------------------------------------------------------
// mrc_dcrds_at_fc

static inline void
mrc_dcrds_at_fc(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       int d, double crd_fc[3])
{
  if (d == 0) {
    // Bx located at i, j+.5, k+.5
    crd_fc[0] = MRC_DMCRDX_NC(crds, ix, p);
    crd_fc[1] = MRC_DMCRDY   (crds, iy, p);
    crd_fc[2] = MRC_DMCRDZ   (crds, iz, p);
  } else if (d == 1) {
    // By located at i+.5, j, k+.5
    crd_fc[0] = MRC_DMCRDX   (crds, ix, p);
    crd_fc[1] = MRC_DMCRDY_NC(crds, iy, p);
    crd_fc[2] = MRC_DMCRDZ   (crds, iz, p);
  } else if (d == 2) {
    // Bz located at i+.5, j+.5, k
    crd_fc[0] = MRC_DMCRDX   (crds, ix, p);
    crd_fc[1] = MRC_DMCRDY   (crds, iy, p);
    crd_fc[2] = MRC_DMCRDZ_NC(crds, iz, p);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_dcrds_at_ec

static inline void
mrc_dcrds_at_ec(struct mrc_crds *crds, int ix, int iy, int iz, int p,
	       int d, double crd_ec[3])
{
  if (d == 0) {
    // Ex located at i+.5, j, k
    crd_ec[0] = MRC_DMCRDX   (crds, ix, p);
    crd_ec[1] = MRC_DMCRDY_NC(crds, iy, p);
    crd_ec[2] = MRC_DMCRDZ_NC(crds, iz, p);
  } else if (d == 1) {
    // Ey located at i, j+.5, k
    crd_ec[0] = MRC_DMCRDX_NC(crds, ix, p);
    crd_ec[1] = MRC_DMCRDY   (crds, iy, p);
    crd_ec[2] = MRC_DMCRDZ_NC(crds, iz, p);
  } else if (d == 2) {
    // Ez located at i, j, k+.5
    crd_ec[0] = MRC_DMCRDX_NC(crds, ix, p);
    crd_ec[1] = MRC_DMCRDY_NC(crds, iy, p);
    crd_ec[2] = MRC_DMCRDZ   (crds, iz, p);
  } else {
    assert(0);
  }
}

#endif

