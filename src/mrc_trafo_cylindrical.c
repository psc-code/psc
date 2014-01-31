
#include <mrc_trafo_cylindrical.h>
#include <mrc_crds.h>
#include <math.h>
  
#define mrc_trafo_cyl(trafo) mrc_to_subobj(trafo, struct mrc_trafo_cylindrical)

static const double dr_hack = 1e-15;  // FIXME ugly sp workaround

// ======================================================================
// cylindrical trafo (analytic)

static void
_trafo_cyl_coord(struct mrc_trafo *trafo, int block, const double xi[3], double xx[3])
{
  
  double r = xi[0], phi = xi[1];// + M_PI * (block == 1);
  if (r < 0) {
    r += dr_hack;
  }
  xx[0] = r * cos(phi);
  xx[1] = r * sin(phi);
  xx[2] = xi[2];
  
}

static void
_trafo_cyl_jac(struct mrc_trafo *trafo, int block, const double xi[3], double *pJ)
{
  
  double r = xi[0];
  if (r < 0) {
    r += dr_hack;
  }
  *pJ = r;
  
}

// FIXME, this is really a broken interface as far as order
// (component vs vector #) is concerned

static void
_trafo_cyl_el(struct mrc_trafo *trafo, int block, const double xi[3], int d, double el[3])
{
  
  double r = xi[0], phi = xi[1];// + M_PI * (block == 1);
  if (r < 0) {
    r += dr_hack;
  }
  switch (d) {
  case 0: 
    el[0] = 1./r * cos(phi); 
    el[1] =      - sin(phi);
    el[2] = 0.;
    break;
  case 1: 
    el[0] = 1./r * sin(phi);
    el[1] =        cos(phi);
    el[2] = 0.;
    break;
  case 2: 
    el[0] = 0.;
    el[1] = 0.;
    el[2] = 1./r;
    break;
  }
  
}

static void
_trafo_cyl_eu(struct mrc_trafo *trafo, int block, const double xi[3], int d, double eu[3])
{
  double r = xi[0], phi = xi[1];// + M_PI * (block == 1);
  if (r < 0) {
    r += dr_hack;
  }
  switch (d) {
  case 0: 
    eu[0] =         cos(phi); 
    eu[1] = -1./r * sin(phi);
    eu[2] =  0.;
    break;
  case 1: 
    eu[0] =         sin(phi);
    eu[1] =  1./r * cos(phi);
    eu[2] =  0.;
    break;
  case 2: 
    eu[0] =  0.;
    eu[1] =  0.;
    eu[2] =  1.;
    break;
  }
}

static void
_trafo_cyl_gam(struct mrc_trafo *trafo, int block, const double xi[3], int i, int k, int l, double *pGam)
{
  double r = xi[0];
  if (r < 0) {
    r += dr_hack;
  }
  switch (i*100+10*k+l) {
  case ( 11): *pGam = -r; break;
  case (101): *pGam = 1./r; break;
  case (110): *pGam = 1./r; break;
  default:    *pGam = 0.;
  }
}

// ======================================================================
// mrc_trafo_cylindrical_ops

struct mrc_trafo_ops mrc_trafo_cylindrical_ops = {
  .name             = "cylindrical",
  ._block_factory   = "cylindrical",
  .calcCRD  = _trafo_cyl_coord,
  .calcJAC  = _trafo_cyl_jac,
  .calcEL   = _trafo_cyl_el,
  .calcEU   = _trafo_cyl_eu,
  .calcGAM  = _trafo_cyl_gam,
};

