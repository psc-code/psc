
#ifndef MRC_TRAFO_H
#define MRC_TRAFO_H

#include <mrc_domain.h>
#include <mrc_obj.h>

#define MRC_CARTESIAN   (0)
#define MRC_SINUSOIDAL  (1)
#define MRC_BUTTERFLY   (2)
#define MRC_CYLINDRICAL (3)

// ----------------------------------------------------------------------
// mrc_trafo
//

#define TRAFO_CRD0(trafo, jx,jy,jz,patch) MRC_D5(trafo->_cc[0], 0, jx,jy,jz, patch)
#define TRAFO_CRD1(trafo, jx,jy,jz,patch) MRC_D5(trafo->_cc[1], 0, jx,jy,jz, patch)
#define TRAFO_CRD2(trafo, jx,jy,jz,patch) MRC_D5(trafo->_cc[2], 0, jx,jy,jz, patch)

#define TRAFO_JAC(trafo, jx,jy,jz,patch)  MRC_D5(trafo->_jac, 0, jx,jy,jz, patch)

#define TRAFO_M3(x,i,j, jx,jy,jz,patch) MRC_D5(x, (i)*3+j, jx,jy,jz, patch)
#define TRAFO_M33(x, i,j,k, jx,jy,jz,patch) MRC_D5(x, ((i)*3+(j))*3+(k), jx,jy,jz, patch)

// i: unit vector nr, j: component (x,y,z)
#define TRAFO_EU(trafo, i,j, jx,jy,jz,patch) TRAFO_M3(trafo->_eu,i,j, jx,jy,jz,patch)
#define TRAFO_EL(trafo, i,j, jx,jy,jz,patch) TRAFO_M3(trafo->_el,i,j, jx,jy,jz,patch)

#define TRAFO_GUU(trafo, i,j, jx,jy,jz,patch) TRAFO_M3(trafo->_guu,i,j, jx,jy,jz,patch)
#define TRAFO_GLL(trafo, i,j, jx,jy,jz,patch) TRAFO_M3(trafo->_gll,i,j, jx,jy,jz,patch)

#define TRAFO_GAM(trafo, i,j,k, jx,jy,jz,patch) TRAFO_M33(trafo->_gam, i,j,k, jx,jy,jz,patch)

// ----------------------------------------------------------------------


struct mrc_trafo {
  struct mrc_obj obj;
  struct mrc_domain *_domain;

  // specific to mrc_trafo_curvilin, actually
  struct mrc_fld *_cc[3];
  struct mrc_fld *_nc;
  struct mrc_fld *_jac;
  struct mrc_fld *_eu;
  struct mrc_fld *_el;
  struct mrc_fld *_guu;
  struct mrc_fld *_gll;
  struct mrc_fld *_gam;
};

MRC_CLASS_DECLARE(mrc_trafo, struct mrc_trafo);

struct mrc_trafo_ops {
  MRC_SUBCLASS_OPS(struct mrc_trafo);
  void (*calcCRD)(struct mrc_trafo *trafo, int block, const double xi[3], double xx[3]);
  void (*calcJAC)(struct mrc_trafo *trafo, int block, const double xi[3], double *pJ);
  void (*calcEU) (struct mrc_trafo *trafo, int block, const double xi[3], int d, double eu[3]);
  void (*calcEL) (struct mrc_trafo *trafo, int block, const double xi[3], int d, double el[3]);
  void (*calcGAM)(struct mrc_trafo *trafo, int block, const double xi[3], int i, int k, int l, double *pGam);

  const char *_block_factory; //< the default block factory to be used with this trafo
  // I think these must be defined if we want to use the finite difference
  // approximations in the crds, but I haven't the foggiest idea how to set
  // up and use them at the moment, and they never seem to be used so we'll
  // keep them around for now.
  struct MB_bc *bc_jac; // FIXME, this is not the right place.
  struct MB_bc *bc_eu;
  struct MB_bc *bc_el;

};

/* int mrc_trafo_ReadH5(struct mrc_domain_mb *mb, hid_t loc, const char *path, */
/* 		    struct mrc_trafo **ptrafo); */
/* int mrc_trafo_WriteH5(struct mrc_trafo *trafo, hid_t loc, const char *name); */

const char *mrc_trafo_block_factory_type(struct mrc_trafo *trafo);
// helpers for actual implementations
//int __mrc_trafo_WriteH5(struct mrc_trafo *trafo, hid_t loc, const char *path);

#endif
