
#include <mrc_trafo_cartesian.h>
#include <mrc_crds.h>

// ======================================================================
// Cartesian trafo (analytic)

static void
_trafo_cart_coord(struct mrc_trafo *trafo, int block, const double xi[3], double xx[3])
{
  xx[0] = xi[0];
  xx[1] = xi[1];
  xx[2] = xi[2];
}

static void
_trafo_cart_jac(struct mrc_trafo *trafo, int block, const double xi[3], double *pJ)
{
  *pJ = 1.;
}

static void
_trafo_cart_elu(struct mrc_trafo *trafo, int block, const double xi[3], int d, double el[3])
{
  switch (d) {
  case 0: el[0] = 1.; el[1] = 0.; el[2] = 0.; break;
  case 1: el[0] = 0.; el[1] = 1.; el[2] = 0.; break;
  case 2: el[0] = 0.; el[1] = 0.; el[2] = 1.; break;
  }
}

static void
_trafo_cart_gam(struct mrc_trafo *trafo, int block, const double xi[3], int i, int k, int l, double *pGam)
{
  *pGam = 0.;
}


// ======================================================================
// mrc_trafo_cartesian_get*
// Deprecated, but mabe still useful at some point.

/* static double */
/* mrc_trafo_cartesian_getCRD(struct mrc_trafo *trafo, int d, int jx, int jy, int jz) */
/* { */
/*   struct mrc_domain *mb = trafo->_domain; */

/*   switch (d) { */
/*   case 0: return F3L(mb->mb_coord->mbc_xi[0], 0, jx,0 ,0 ); */
/*   case 1: return F3L(mb->mb_coord->mbc_xi[1], 0, 0 ,jy,0 ); */
/*   case 2: return F3L(mb->mb_coord->mbc_xi[2], 0, 0 ,0 ,jz); */
/*   } */
/*   ASSERT(0); */
/* } */

/* static double */
/* mrc_trafo_cartesian_getJAC(struct mrc_trafo *trafo, int jx, int jy, int jz) */
/* { */
/*   return 1.; */
/* } */

/* static double */
/* mrc_trafo_cartesian_getEU(struct mrc_trafo *trafo, int i, int j, */
/* 			int jx, int jy, int jz) */
/* { */
/*   return i == j ? 1. : 0.; */
/* } */

/* static double */
/* mrc_trafo_cartesian_getEL(struct mrc_trafo *trafo, int i, int j, */
/* 			int jx, int jy, int jz) */
/* { */
/*   return i == j ? 1. : 0.; */
/* } */

/* static double */
/* mrc_trafo_cartesian_getGUU(struct mrc_trafo *trafo, int i, int j, */
/* 			 int jx, int jy, int jz) */
/* { */
/*   return i == j ? 1. : 0.; */
/* } */

/* static double */
/* mrc_trafo_cartesian_getGLL(struct mrc_trafo *trafo, int i, int j, */
/* 			 int jx, int jy, int jz) */
/* { */
/*   return i == j ? 1. : 0.; */
/* } */

/* static double */
/* mrc_trafo_cartesian_getGAM(struct mrc_trafo *trafo, int i, int j, int k, */
/* 			 int jx, int jy, int jz) */
/* { */
/*   return 0.; */
/* } */

// ======================================================================
// mrc_trafo_cartesian_ops

struct mrc_trafo_ops mrc_trafo_cartesian_ops = {
  .name             = "cartesian",
  .size             = sizeof(struct mrc_trafo_ops),
  ._block_factory   = "simple2d",
  .calcCRD  = _trafo_cart_coord,
  .calcJAC  = _trafo_cart_jac,
  .calcEU   = _trafo_cart_elu,
  .calcEL   = _trafo_cart_elu,
  .calcGAM  = _trafo_cart_gam,
};

