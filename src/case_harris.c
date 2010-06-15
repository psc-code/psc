
#include "psc.h"
#include "util/params.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

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

struct psc_harris {
  double BB;
  double nnb;
  double Te, Ti;
  double MMi;
  double lambda;
  double lx, lz;
  double pert;
};

#define VAR(x) (void *)offsetof(struct psc_harris, x)

static struct param psc_harris_descr[] = {
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
harris_create()
{
  struct psc_harris *harris = malloc(sizeof(*harris));
  memset(harris, 0, sizeof(*harris));

  params_parse_cmdline(harris, psc_harris_descr, "PSC Harris", MPI_COMM_WORLD);
  params_print(harris, psc_harris_descr, "PSC Harris", MPI_COMM_WORLD);

  psc.case_data = harris;
}

static void
harris_destroy()
{
  free(psc.case_data);
  psc.case_data = NULL;
}

static void
harris_init_param()
{
  struct psc_harris *harris = psc.case_data;

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

  psc.domain.length[0] = harris->lx * sqrt(harris->MMi);
  psc.domain.length[1] = 10000000.; // no y dependence 
  psc.domain.length[2] = 2. * harris->lz * sqrt(harris->MMi); // double tearing

  psc.domain.itot[0] = 640;
  psc.domain.itot[1] = 1;
  psc.domain.itot[2] = 640;
  psc.domain.ilo[0] = 0;
  psc.domain.ilo[1] = 0;
  psc.domain.ilo[2] = 0;
  psc.domain.ihi[0] = 640;
  psc.domain.ihi[1] = 1;
  psc.domain.ihi[2] = 640;

  psc.domain.bnd_fld[0] = 1;
  psc.domain.bnd_fld[1] = 1;
  psc.domain.bnd_fld[2] = 1;
  psc.domain.bnd_part[0] = 1;
  psc.domain.bnd_part[1] = 1;
  psc.domain.bnd_part[2] = 1;
}

struct psc_case_ops psc_case_ops_harris = {
  .name       = "harris",
  .create     = harris_create,
  .destroy    = harris_destroy,
  .init_param = harris_init_param,
};
