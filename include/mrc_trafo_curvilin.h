
#ifndef MRC_TRAFO_CURVILIN_H
#define MRC_TRAFO_CURVILIN_H

#include <mrc_trafo.h>

// ----------------------------------------------------------------------
// mrc_trafo_curvilin
//
// curvilinear trafo -- physical coords

#define CL_CRD0(jx,jy,jz) F3(mrc_domain_mb(__mb)->_trafo->_cc[0], 0, jx,jy,jz)
#define CL_CRD1(jx,jy,jz) F3(mrc_domain_mb(__mb)->_trafo->_cc[1], 0, jx,jy,jz)
#define CL_CRD2(jx,jy,jz) F3(mrc_domain_mb(__mb)->_trafo->_cc[2], 0, jx,jy,jz)

#define CL_JAC(jx,jy,jz)  F3(mrc_domain_mb(__mb)->_trafo->_jac,0, jx,jy,jz)

#define CL_M3(x,i,j, jx,jy,jz) F3(x,(i)*3+j, jx,jy,jz)
#define CL_M33(x, i,j,k, jx,jy,jz) F3(x, ((i)*3+(j))*3+(k), jx,jy,jz)

// i: unit vector nr, j: component (x,y,z)
#define CL_EU(i,j, jx,jy,jz) CL_M3(mrc_domain_mb(__mb)->_trafo->_eu,i,j, jx,jy,jz)
#define CL_EL(i,j, jx,jy,jz) CL_M3(mrc_domain_mb(__mb)->_trafo->_el,i,j, jx,jy,jz)

#define CL_GUU(i,j, jx,jy,jz) CL_M3(mrc_domain_mb(__mb)->_trafo->_guu,i,j, jx,jy,jz)
#define CL_GLL(i,j, jx,jy,jz) CL_M3(mrc_domain_mb(__mb)->_trafo->_gll,i,j, jx,jy,jz)

#define CL_GAM(i,j,k, jx,jy,jz) CL_M33(mrc_domain_mb(__mb)->_trafo->_gam, i,j,k, jx,jy,jz)


#endif
