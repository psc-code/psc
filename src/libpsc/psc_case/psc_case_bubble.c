
#include "psc.h"
#include "psc_case_private.h"
#include "psc_push_fields.h"
#include "psc_bnd_fields.h"
#include "psc_sort.h"
#include "psc_fields_as_c.h"
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// FIXME (description), below parameters don't include scaling factors
// BB:  peak bubble magnetic field  (magnitude gives ratio w_ce/w_pe)
// nnb: number of background particles (density max == 1)
// TTe,TTi:  bulk temperature of electrons and ions (units m_e c^2)
// MMi: ion mass / electron mass
// MMach: Mach number of initial bubble flow
// LLn, LLB: length scales for initial bubble (units of de)

struct psc_case_bubble {
  double BB;
  double nnb;
  double TTe;
  double TTi;
  double MMi;
  double MMach;
  double LLn;
  double LLB;
};

#define VAR(x) (void *)offsetof(struct psc_case_bubble, x)

static struct param psc_case_bubble_descr[] = {
  { "BB"            , VAR(BB)              , PARAM_DOUBLE(.07)    },
  { "nnb"           , VAR(nnb)             , PARAM_DOUBLE(.1)     },
  { "MMi"           , VAR(MMi)             , PARAM_DOUBLE(100.)   },
  { "MMach"         , VAR(MMach)           , PARAM_DOUBLE(3.)     },
  { "LLn"           , VAR(LLn)             , PARAM_DOUBLE(200.)   },
  { "LLB"           , VAR(LLB)             , PARAM_DOUBLE(200./6.)},
  { "TTe"           , VAR(TTe)             , PARAM_DOUBLE(.02)    },
  { "TTi"           , VAR(TTi)             , PARAM_DOUBLE(.02)    },
  {},
};

#undef VAR

static void
psc_case_bubble_set_from_options(struct psc_case *_case)
{
  struct psc_case_bubble *bubble = mrc_to_subobj(_case, struct psc_case_bubble);

  ppsc->prm.qq = 1.;
  ppsc->prm.mm = 1.;
  ppsc->prm.tt = 1.;
  ppsc->prm.cc = 1.;
  ppsc->prm.eps0 = 1.;

  ppsc->prm.nmax = 32000;
  ppsc->prm.cpum = 5*24.0*60*60;
  ppsc->prm.lw = 2.*M_PI;
  ppsc->prm.i0 = 0.;
  ppsc->prm.n0 = 1.;
  ppsc->prm.e0 = 1.;

  ppsc->prm.nicell = 50;

  ppsc->domain.length[0] = 3. * bubble->LLn;
  ppsc->domain.length[1] = 10. * bubble->LLn; // no y dependence 
  ppsc->domain.length[2] = 2. * bubble->LLn;

  ppsc->domain.corner[0] = -1.5 * bubble->LLn;
  ppsc->domain.corner[1] = -5. * bubble->LLn;
  ppsc->domain.corner[2] = -1. * bubble->LLn;

  ppsc->domain.gdims[0] = 2400;
  ppsc->domain.gdims[1] = 1;
  ppsc->domain.gdims[2] = 1600;

  ppsc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_part[0] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part[1] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part[2] = BND_PART_PERIODIC;

  struct psc_bnd_fields *bnd_fields = 
    psc_push_fields_get_bnd_fields(ppsc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "none");

  psc_sort_set_type(ppsc->sort, "countsort2");
}

static void
psc_case_bubble_init_field(struct psc_case *_case, mfields_t *flds)
{
  struct psc_case_bubble *bubble = mrc_to_subobj(_case, struct psc_case_bubble);
  struct psc *psc = _case->psc;

  double BB = bubble->BB;
  double LLn = bubble->LLn;
  double LLB = bubble->LLB;
  double MMi = bubble->MMi;
  double MMach = bubble->MMach;
  double TTe = bubble->TTe;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_t *pf = psc_mfields_get_patch(flds, p);
    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dx = psc->dx[0], dz = psc->dx[2];
      double xx = CRDX(p, jx), zz = CRDZ(p, jz);

      double x1 = xx;
      double z1 = zz + 0.5*dz + LLn;
      double r1 = sqrt(sqr(x1) + sqr(z1));
      double x2 = xx;
      double z2 = zz + 0.5*dz - LLn;
      double r2 = sqrt(sqr(x2) + sqr(z2));

      if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	F3(pf, HX, jx,jy,jz) += - BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * z1 / r1;
      }
      if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	F3(pf, HX, jx,jy,jz) += - BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * z2 / r2;
      }

      x1 = xx + 0.5*dx;
      z1 = zz + LLn;
      r1 = sqrt(sqr(x1) + sqr(z1));
      x2 = xx + 0.5*dx;
      z2 = zz - LLn;
      r2 = sqrt(sqr(x2) + sqr(z2));

      if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	F3(pf, HZ, jx,jy,jz) += BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * x1 / r1;
      }
      if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	F3(pf, HZ, jx,jy,jz) += BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * x2 / r2;
      }

      x1 = xx;
      z1 = zz + LLn;
      r1 = sqrt(sqr(x1) + sqr(z1));
      x2 = xx;
      z2 = zz - LLn;
      r2 = sqrt(sqr(x2) + sqr(z2));

      if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	F3(pf, EY, jx,jy,jz) += MMach * sqrt(TTe/MMi) * BB *
	  sin(M_PI * (LLn - r1)/(2.*LLB)) *
	  sin(M_PI * r1 / LLn);
	F3(pf, JYI, jx,jy,jz) += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r1)/(2.*LLB));
      }
      if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	F3(pf, EY, jx,jy,jz) += MMach * sqrt(TTe/MMi) * BB *
	  sin(M_PI * (LLn - r2)/(2.*LLB)) *
	  sin(M_PI * r2 / LLn);
	F3(pf, JYI, jx,jy,jz) += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r2)/(2.*LLB));
      }
    } foreach_3d_g_end;
  }
}

static void
psc_case_bubble_init_npt(struct psc_case *_case, int kind, double x[3],
			 struct psc_particle_npt *npt)
{
  struct psc_case_bubble *bubble = mrc_to_subobj(_case, struct psc_case_bubble);

  double BB = bubble->BB;
  double LLn = bubble->LLn;
  double LLB = bubble->LLB;
  double V0 = bubble->MMach * sqrt(bubble->TTe / bubble->MMi);

  double nnb = bubble->nnb;
  double TTi = bubble->TTi;
  double TTe = bubble->TTe;

  double r1 = sqrt(sqr(x[0]) + sqr(x[2] + LLn));
  double r2 = sqrt(sqr(x[0]) + sqr(x[2] - LLn));

  npt->n = nnb;
  if (r1 < LLn) {
    npt->n += (1. - nnb) * sqr(cos(M_PI / 2. * r1 / LLn));
    if (r1 > 0.0) {
      npt->p[0] += V0 * sin(M_PI * r1 / LLn) * x[0] / r1;
      npt->p[2] += V0 * sin(M_PI * r1 / LLn) * (x[2] + 1.*LLn) / r1;
    }
  }
  if (r2 < LLn) {
    npt->n += (1. - nnb) * sqr(cos(M_PI / 2. * r2 / LLn));
    if (r2 > 0.0) {
      npt->p[0] += V0 * sin(M_PI * r2 / LLn) * x[0] / r2;
      npt->p[2] += V0 * sin(M_PI * r2 / LLn) * (x[2] - 1.*LLn) / r2;
    }
  }

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;

    // electron drift consistent with initial current
    if ((r1 <= LLn) && (r1 >= LLn - 2.*LLB)) {
      npt->p[1] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r1)/(2.*LLB)) / npt->n;
    }
    if ((r2 <= LLn) && (r2 >= LLn - 2.*LLB)) {
      npt->p[1] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r2)/(2.*LLB)) / npt->n;
    }

    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = bubble->MMi;

    npt->T[0] = TTi;
    npt->T[1] = TTi;
    npt->T[2] = TTi;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_bubble_ops = {
  .name             = "bubble",
  .size             = sizeof(struct psc_case_bubble),
  .param_descr      = psc_case_bubble_descr,
  .set_from_options = psc_case_bubble_set_from_options,
  .init_field       = psc_case_bubble_init_field,
  .init_npt         = psc_case_bubble_init_npt,
};

