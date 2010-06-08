
#include "psc.h"

#include <stdlib.h>

struct psc_harris {
};

static void
harris_create()
{
  psc.case_data = malloc(sizeof(struct psc_harris));
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
  psc.prm.qq = 1.;
  psc.prm.mm = 1.;
  psc.prm.tt = 1.;
  psc.prm.cc = 1.;
  psc.prm.eps0 = 1.;

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

  psc.prm.nmax = 16000;
  psc.prm.cpum = 5*24.0*60*60;
  psc.prm.lw = 2.*M_PI;
  psc.prm.i0 = 0.;
  psc.prm.n0 = 1.;
  psc.prm.e0 = 1.;
}

struct psc_case_ops psc_case_ops_harris = {
  .name       = "harris",
  .create     = harris_create,
  .destroy    = harris_destroy,
  .init_param = harris_init_param,
};
