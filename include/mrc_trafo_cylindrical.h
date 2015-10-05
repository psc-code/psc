
#ifndef MRC_TRAFO_CYLINDRICAL_H
#define MRC_TRAFO_CYLINDRICAL_H

#include <mrc_trafo.h>

// ----------------------------------------------------------------------
// mrc_trafo_cylindrical
//
// cylindrical trafo -- physical coords

#define CY_CRD0(jx,jy,jz) F3(mrc_domain_mb(__mb)->_trafo->_cc[0], 0, jx,jy,jz)
#define CY_CRD1(jx,jy,jz) F3(mrc_domain_mb(__mb)->_trafo->_cc[1], 0, jx,jy,jz)
#define CY_CRD2(jx,jy,jz) F3(mrc_domain_mb(__mb)->_trafo->_cc[2], 0, jx,jy,jz)

#define CY_JAC(jx,jy,jz)  F3(mrc_domain_mb(__mb)->_trafo->_jac,0, jx,jy,jz)

#define CY_M3(x,i,j, jx,jy,jz) F3(x,(i)*3+j, jx,jy,jz)
#define CY_M33(x, i,j,k, jx,jy,jz) F3(x, ((i)*3+(j))*3+(k), jx,jy,jz)

// i: unit vector nr, j: component (x,y,z)
#define CY_EU(i,j, jx,jy,jz) CY_M3(mrc_domain_mb(__mb)->_trafo->_eu,i,j, jx,jy,jz)
#define CY_EL(i,j, jx,jy,jz) CY_M3(mrc_domain_mb(__mb)->_trafo->_el,i,j, jx,jy,jz)

#define CY_GUU(i,j, jx,jy,jz) CY_M3(mrc_domain_mb(__mb)->_trafo->_guu,i,j, jx,jy,jz)
#define CY_GLL(i,j, jx,jy,jz) CY_M3(mrc_domain_mb(__mb)->_trafo->_gll,i,j, jx,jy,jz)

#define CY_GAM(i,j,k, jx,jy,jz) CY_M33(mrc_domain_mb(__mb)->_trafo->_gam, i,j,k, jx,jy,jz)


#endif
