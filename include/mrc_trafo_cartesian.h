
#ifndef MRC_TRAFO_CARTESIAN_H
#define MRC_TRAFO_CARTESIAN_H

#include <mrc_trafo.h>

// ----------------------------------------------------------------------
// mrc_trafo_cartesian
//
// cartesian (identity) trafo -- physical coords
// Umm.. I'm pretty sure this isn't right...
// Maybe these just never get used? Code gen
// uses XI0 instead...
#define CA_CRD0(trafo,jx,jy,jz,patch) XI0(mrc_domain_get_crds(trafo->_domain),jx,patch)
#define CA_CRD1(trafo,jx,jy,jz,patch) XI1(mrc_domain_get_crds(trafo->_domain),jx,patch)
#define CA_CRD2(trafo,jx,jy,jz,patch) XI2(mrc_domain_get_crds(trafo->_domain),jx,patch)

#define CA_JAC(trafo,jx,jy,jz,patch)  (1.)

// i: unit vector nr, j: component (x,y,z)
#define CA_EU(trafo,i,j, jx,jy,jz,patch) (i == j ? 1. : 0.)
#define CA_EL(trafo,i,j, jx,jy,jz,patch) (i == j ? 1. : 0.)
#define CA_GUU(trafo,i,j, jx,jy,jz,patch) (i == j ? 1. : 0.)
#define CA_GLL(trafo,i,j, jx,jy,jz,patch) (i == j ? 1. : 0.)
#define CA_GAM(trafo,i,j,k, jx,jy,jz,patch) (0.)

#endif
