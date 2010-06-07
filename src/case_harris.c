
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
}

struct psc_case_ops psc_case_ops_harris = {
  .name       = "harris",
  .create     = harris_create,
  .destroy    = harris_destroy,
  .init_param = harris_init_param,
};
