
#ifndef MRC_TRAFO_FAST_H
#define MRC_TRAFO_FAST_H

#include <mrc_trafo.h>

#if MRC_TRAFO == MRC_CARTESIAN
// If we're cartesian, we can use some quick logic instead
// of actually reading the trafo information.

#include <mrc_trafo_cartesian.h>
#define CRD0(trafo,jx,jy,jz,patch) CA_CRD0(trafo,jx,jy,jz,patch)
#define CRD1(trafo,jx,jy,jz,patch) CA_CRD1(trafo,jx,jy,jz,patch)
#define CRD2(trafo,jx,jy,jz,patch) CA_CRD2(trafo,jx,jy,jz,patch)

#define JAC(trafo,jx,jy,jz,patch) CA_JAC(trafo,jx,jy,jz,patch)
#define EU(trafo,i,j, jx,jy,jz,patch) CA_EU(trafo,i,j, jx,jy,jz,patch)
#define EL(trafo,i,j, jx,jy,jz,patch) CA_EL(trafo,i,j, jx,jy,jz,patch)
#define GUU(trafo,i,j, jx,jy,jz,patch) CA_GUU(trafo,i,j, jx,jy,jz,patch)
#define GLL(trafo,i,j, jx,jy,jz,patch) CA_GLL(trafo,i,j, jx,jy,jz,patch)
#define GAM(trafo,i,j,k, jx,jy,jz,patch) CA_GAM(trafo,i,j,k, jx,jy,jz,patch)

#elif MRC_TRAFO == MRC_CYLINDRICAL || MRC_TRAFO == MRC_BUTTERFLY

// otherwise we need to read the full trafo.

#define CRD0(trafo,jx,jy,jz,patch) TRAFO_CRD0(trafo,jx,jy,jz,patch)
#define CRD1(trafo,jx,jy,jz,patch) TRAFO_CRD1(trafo,jx,jy,jz,patch)
#define CRD2(trafo,jx,jy,jz,patch) TRAFO_CRD2(trafo,jx,jy,jz,patch)

#define JAC(trafo,jx,jy,jz,patch) TRAFO_JAC(trafo,jx,jy,jz,patch)
#define EU(trafo,i,j, jx,jy,jz,patch) TRAFO_EU(trafo,i,j, jx,jy,jz,patch)
#define EL(trafo,i,j, jx,jy,jz,patch) TRAFO_EL(trafo,i,j, jx,jy,jz,patch)
#define GUU(trafo,i,j, jx,jy,jz,patch) TRAFO_GUU(trafo,i,j, jx,jy,jz,patch)
#define GLL(trafo,i,j, jx,jy,jz,patch) TRAFO_GLL(trafo,i,j, jx,jy,jz,patch)
#define GAM(trafo,i,j,k, jx,jy,jz,patch) TRAFO_GAM(trafo,i,j,k, jx,jy,jz,patch)

#else
#error unknown MRC_TRAFO
#endif

#endif
