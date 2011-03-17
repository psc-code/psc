
#include "psc.h"
#include "psc_case_private.h"
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// WFox: Plasma simulation parameters
//       needed because they determine length scales for initial conditions
// BB:  peak Harris magnetic field  (magnitude gives ratio w_ce/w_pe)
// nnb: number of background particles (density max == 1)
// TTe,TTi:  bulk temperature of electrons and ions (units m_e c^2)
// MMi: ion mass / electron mass
// LLL = reversal length scale (units of c/wpe)
// LLz, LLx = simulation box size (units of c/wpe)
// AA = perturbation (units of B * de)

// FIXME (description), below parameters don't include scaling factors

struct psc_case_harris {
  double BB;
  double nnb;
  double Te, Ti;
  double MMi;
  double lambda;
  double lx, lz;
  double pert;
};

#define VAR(x) (void *)offsetof(struct psc_case_harris, x)

static struct param psc_case_harris_descr[] = {
  { "BB"            , VAR(BB)              , PARAM_DOUBLE(1.)     },
  { "MMi"           , VAR(MMi)             , PARAM_DOUBLE(25.)    },
  { "nnb"           , VAR(nnb)             , PARAM_DOUBLE(.2)     },
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(1./12.) },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(5./12.) },
  { "lambda"        , VAR(lambda)          , PARAM_DOUBLE(.5)     },
  { "lx"            , VAR(lx)              , PARAM_DOUBLE(25.6)   },
  { "lz"            , VAR(lz)              , PARAM_DOUBLE(12.8)   },
  { "pert"          , VAR(pert)            , PARAM_DOUBLE(.1)     },
  {},
};

#undef VAR

static void
psc_case_harris_set_from_options(struct psc_case *_case)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);

  psc.prm.qq = 1.;
  psc.prm.mm = 1.;
  psc.prm.tt = 1.;
  psc.prm.cc = 1.;
  psc.prm.eps0 = 1.;

  psc.prm.nmax = 16000;
  psc.prm.cpum = 5*24.0*60*60;
  psc.prm.lw = 2.*M_PI;
  psc.prm.i0 = 0.;
  psc.prm.n0 = 1.;
  psc.prm.e0 = 1.;

  psc.prm.nicell = 50;

  real d_i = sqrt(harris->MMi); // in units of d_e
  psc.domain.length[0] = harris->lx * d_i;
  psc.domain.length[1] = 1.; // no y dependence 
  psc.domain.length[2] = 2. * harris->lz * d_i; // double tearing

  psc.domain.gdims[0] = 640;
  psc.domain.gdims[1] = 1;
  psc.domain.gdims[2] = 640;

  psc.domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc.domain.bnd_part[0] = BND_PART_PERIODIC;
  psc.domain.bnd_part[1] = BND_PART_PERIODIC;
  psc.domain.bnd_part[2] = BND_PART_PERIODIC;
}

static void
psc_case_harris_init_field(struct psc_case *_case, mfields_base_t *flds)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);
  struct psc *psc = _case->psc;

  real d_i = sqrt(harris->MMi); // in units of d_e
  double BB = harris->BB;
  double LLx = harris->lx * d_i, LLz = harris->lz * d_i;
  double LLL = harris->lambda * d_i;
  double AA = harris->pert * BB * d_i;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dx = psc->dx[0], dz = psc->dx[2];
      double xx = CRDX(p, jx), zz = CRDZ(p, jz);
    
      F3_BASE(pf, HX, jx,jy,jz) = 
	BB * (-1. 
	      + tanh((zz + .5*dz - 0.5*LLz)/LLL)
	      - tanh((zz + .5*dz - 1.5*LLz)/LLL))
	+ AA*M_PI/LLz * sin(2.*M_PI*xx/LLx) * cos(M_PI*(zz+.5*dz)/LLz);
      
      F3_BASE(pf, HZ, jx,jy,jz) =
	- AA*2.*M_PI/LLx * cos(2.*M_PI*(xx+.5*dx)/LLx) * sin(M_PI*zz/LLz);
      
      F3_BASE(pf, JYI, jx,jy,jz) = BB/LLL *
	(1./sqr(cosh((zz - 0.5*LLz)/LLL)) -1./sqr(cosh((zz - 1.5*LLz)/LLL)))
	- (AA*sqr(M_PI) * (1./sqr(LLz) + 4./sqr(LLx)) 
	   * sin(2.*M_PI*xx/LLx) * sin(M_PI*zz/LLz));
    } foreach_3d_g_end;
  }
}

static void
psc_case_harris_init_npt(struct psc_case *_case, int kind, double x[3],
			 struct psc_particle_npt *npt)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);

  real d_i = sqrt(harris->MMi); // in units of d_e
  double BB = harris->BB;
  double LLz = harris->lz * d_i;
  double LLL = harris->lambda * d_i;
  double nnb = harris->nnb;
  double TTi = harris->Ti * sqr(BB);
  double TTe = harris->Te * sqr(BB);

  double jy0 = 1./sqr(cosh((x[2]-0.5*LLz)/LLL)) - 1./sqr(cosh((x[2]-1.5*LLz)/LLL));

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->n = nnb + 1./sqr(cosh((x[2]-0.5*LLz)/LLL)) + 1./sqr(cosh((x[2]-1.5*LLz)/LLL));
    npt->p[1] = - 2. * TTe / BB / LLL * jy0 / npt->n;
    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = harris->MMi;
    npt->n = nnb + 1./sqr(cosh((x[2]-0.5*LLz)/LLL)) + 1./sqr(cosh((x[2]-1.5*LLz)/LLL));
    npt->p[1] = 2. * TTi / BB / LLL * jy0 / npt->n;
    npt->T[0] = TTi;
    npt->T[1] = TTi;
    npt->T[2] = TTi;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_harris_ops = {
  .name             = "harris",
  .size             = sizeof(struct psc_case_harris),
  .param_descr      = psc_case_harris_descr,
  .set_from_options = psc_case_harris_set_from_options,
  .init_field       = psc_case_harris_init_field,
  .init_npt         = psc_case_harris_init_npt,
};

// ----------------------------------------------------------------------
// case test_xz:
//
// basically the same as harris

static void
psc_case_test_xz_set_from_options(struct psc_case *_case)
{
  psc_case_harris_set_from_options(_case);
  
  psc.prm.nicell = 100;

  psc.domain.gdims[0] = 64;
  psc.domain.gdims[1] = 1;
  psc.domain.gdims[2] = 64;
  
}

struct psc_case_ops psc_case_test_xz_ops = {
  .name             = "test_xz",
  .size             = sizeof(struct psc_case_harris),
  .param_descr      = psc_case_harris_descr,
  .set_from_options = psc_case_test_xz_set_from_options,
  .init_field       = psc_case_harris_init_field,
  .init_npt         = psc_case_harris_init_npt,
};

// ----------------------------------------------------------------------
// case test_yz:
//
// basically the same as harris, but with different coordinates

static void
psc_case_test_yz_set_from_options(struct psc_case *_case)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);

  psc_case_harris_set_from_options(_case);
  
  psc.prm.nicell = 100;

  real d_i = sqrt(harris->MMi); // in units of d_e
  psc.domain.length[0] = 10000000;
  psc.domain.length[1] = harris->lx * d_i;
  psc.domain.length[2] = 2. * harris->lz * d_i; // double tearing

  psc.domain.gdims[0] = 1;
  psc.domain.gdims[1] = 64;
  psc.domain.gdims[2] = 64;
}

static void
psc_case_test_yz_init_field(struct psc_case *_case, mfields_base_t *flds)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);
  struct psc *psc = _case->psc;

  real d_i = sqrt(harris->MMi); // in units of d_e
  double BB = harris->BB;
  double LLy = harris->lx * d_i, LLz = harris->lz * d_i;
  double LLL = harris->lambda * d_i;
  double AA = harris->pert * BB * d_i;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dy = psc->dx[1], dz = psc->dx[2];
      double yy = CRDY(p, jy), zz = CRDZ(p, jz);
      
      F3_BASE(pf, HY, jx,jy,jz) = 
	BB * (-1. 
	      + tanh((zz + .5*dz - 0.5*LLz)/LLL)
	      - tanh((zz + .5*dz - 1.5*LLz)/LLL))
	+ AA*M_PI/LLz * sin(2.*M_PI*yy/LLy) * cos(M_PI*(zz+.5*dz)/LLz);
      
      F3_BASE(pf, HZ, jx,jy,jz) =
	- AA*2.*M_PI/LLy * cos(2.*M_PI*(yy+.5*dy)/LLy) * sin(M_PI*zz/LLz);
      
      F3_BASE(pf, JXI, jx,jy,jz) = - BB/LLL *
	(1./sqr(cosh((zz - 0.5*LLz)/LLL)) -1./sqr(cosh((zz - 1.5*LLz)/LLL)))
	- (AA*sqr(M_PI) * (1./sqr(LLz) + 4./sqr(LLy)) 
	   * sin(2.*M_PI*yy/LLy) * sin(M_PI*zz/LLz));
    } foreach_3d_g_end;
  }
}

static void
psc_case_test_yz_init_npt(struct psc_case *_case, int kind, double x[3],
			   struct psc_particle_npt *npt)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);

  real d_i = sqrt(harris->MMi); // in units of d_e
  double BB = harris->BB;
  double LLz = harris->lz * d_i;
  double LLL = harris->lambda * d_i;
  double nnb = harris->nnb;
  double TTi = harris->Ti * sqr(BB);
  double TTe = harris->Te * sqr(BB);

  double jx0 = -1./sqr(cosh((x[2]-0.5*LLz)/LLL)) - 1./sqr(cosh((x[2]-1.5*LLz)/LLL));

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->n = nnb + 1./sqr(cosh((x[2]-0.5*LLz)/LLL)) + 1./sqr(cosh((x[2]-1.5*LLz)/LLL));
    npt->p[0] = - 2. * TTe / BB / LLL * jx0 / npt->n;
    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = harris->MMi;
    npt->n = nnb + 1./sqr(cosh((x[2]-0.5*LLz)/LLL)) + 1./sqr(cosh((x[2]-1.5*LLz)/LLL));
    npt->p[0] = 2. * TTi / BB / LLL * jx0 / npt->n;
    npt->T[0] = TTi;
    npt->T[1] = TTi;
    npt->T[2] = TTi;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_test_yz_ops = {
  .name                  = "test_yz",
  .size                  = sizeof(struct psc_case_harris),
  .param_descr           = psc_case_harris_descr,
  .set_from_options      = psc_case_test_yz_set_from_options,
  .init_field            = psc_case_test_yz_init_field,
  .init_npt              = psc_case_test_yz_init_npt,
};


// ----------------------------------------------------------------------
// case test_xy:
//
// basically the same as harris, but with different coordinates

static void
psc_case_test_xy_set_from_options(struct psc_case *_case)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);

  psc_case_harris_set_from_options(_case);

  psc.prm.nicell = 200;

  real d_i = sqrt(harris->MMi); // in units of d_e
  psc.domain.length[0] = 2. * harris->lz * d_i; // double tearing
  psc.domain.length[1] = harris->lx * d_i;
  psc.domain.length[2] = 10000000;

  // I hacked this in to test the cell pusher, which is why 
  // the domain is sort of funky shaped. 
  psc.domain.gdims[0] = 26;
  psc.domain.gdims[1] = 26;
  psc.domain.gdims[2] = 1;
  
}

psc_case_test_xy_init_field(struct psc_case *_case, mfields_base_t *flds)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);
  struct psc *psc = _case->psc;

  real d_i = sqrt(harris->MMi); // in units of d_e
  double BB = harris->BB;
  double LLy = harris->lx * d_i, LLx = harris->lz * d_i;
  double LLL = harris->lambda * d_i;
  double AA = harris->pert * BB * d_i;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    psc_foreach_3d_g(psc,p, jx, jy, jz) {
      double dy = psc.dx[1], dx = psc.dx[0];
      double yy = CRDY(p, jy), xx = CRDX(p, jx);
      
      F3_BASE(pf, HY, jx,jy,jz) = 
	BB * (-1. 
	      + tanh((xx + .5*dx - 0.5*LLx)/LLL)
	      - tanh((xx + .5*dx - 1.5*LLx)/LLL))
	+ AA*M_PI/LLx * sin(2.*M_PI*yy/LLy) * cos(M_PI*(xx+.5*dx)/LLx);
      
      F3_BASE(pf, HZ, jx,jy,jz) =
	- AA*2.*M_PI/LLy * cos(2.*M_PI*(yy+.5*dy)/LLy) * sin(M_PI*xx/LLx);
      
      F3_BASE(pf, JXI, jx,jy,jz) = - BB/LLL *
	(1./sqr(cosh((xx - 0.5*LLx)/LLL)) -1./sqr(cosh((xx - 1.5*LLx)/LLL)))
	- (AA*sqr(M_PI) * (1./sqr(LLx) + 4./sqr(LLy)) 
	   * sin(2.*M_PI*yy/LLy) * sin(M_PI*xx/LLx));
    } psc_foreach_3d_g_end;
  }
}

static void
psc_case_test_xy_init_npt(struct psc_case *_case, int kind, double x[3],
			   struct psc_particle_npt *npt)
{
  struct psc_case_harris *harris = mrc_to_subobj(_case, struct psc_case_harris);


  real d_i = sqrt(harris->MMi); // in units of d_e
  double BB = harris->BB;
  double LLx = harris->lz * d_i;
  double LLL = harris->lambda * d_i;
  double nnb = harris->nnb;
  double TTi = harris->Ti * sqr(BB);
  double TTe = harris->Te * sqr(BB);

  double jz0 = -1./sqr(cosh((x[0]-0.5*LLx)/LLL)) - 1./sqr(cosh((x[0]-1.5*LLx)/LLL));

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->n = nnb + 1./sqr(cosh((x[0]-0.5*LLx)/LLL)) + 1./sqr(cosh((x[0]-1.5*LLx)/LLL));
    npt->p[2] = - 2. * TTe / BB / LLL * jz0 / npt->n;
    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = harris->MMi;
    npt->n = nnb + 1./sqr(cosh((x[0]-0.5*LLx)/LLL)) + 1./sqr(cosh((x[0]-1.5*LLx)/LLL));
    npt->p[2] = 2. * TTi / BB / LLL * jz0 / npt->n;
    npt->T[0] = TTi;
    npt->T[1] = TTi;
    npt->T[2] = TTi;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_test_xy_ops = {
  .name                  = "test_xy",
  .size                  = sizeof(struct psc_case_harris),
  .param_descr           = psc_case_harris_descr,
  .set_from_options      = psc_case_test_xy_set_from_options,
  .init_field            = psc_case_test_xy_init_field,
  .init_npt              = psc_case_test_xy_init_npt,
};


