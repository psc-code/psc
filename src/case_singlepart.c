
#include "psc.h"
#include "util/params.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// FIXME description

struct psc_singlepart {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct psc_singlepart, x)

static struct param psc_singlepart_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.)             },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.)             },
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(10.0 * 1e-6)    },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(20.0 * 1e-6)    },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(10.0 * 1e-6)    },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(1836.)          },
  {},
};

#undef VAR

static void
singlepart_create()
{
  struct psc_singlepart *singlepart = malloc(sizeof(*singlepart));
  memset(singlepart, 0, sizeof(*singlepart));

  params_parse_cmdline(singlepart, psc_singlepart_descr, "PSC Singlepart", MPI_COMM_WORLD);
  params_print(singlepart, psc_singlepart_descr, "PSC Singlepart", MPI_COMM_WORLD);

  psc.case_data = singlepart;
}

static void
singlepart_destroy()
{
  free(psc.case_data);
  psc.case_data = NULL;
}

static void
singlepart_init_param()
{
  psc.prm.nmax = 1000;
  psc.prm.cpum = 20000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 2.0e22;
  psc.prm.n0 = 1.0e25;

  psc.prm.nicell = 1;

  psc.domain.length[0] = 20. * 1e-6;
  psc.domain.length[1] = 40. * 1e-6;
  psc.domain.length[2] = 60. * 1e-6;

  psc.domain.itot[0] = 200;
  psc.domain.itot[1] = 400;
  psc.domain.itot[2] = 600;
  psc.domain.ilo[0] = 99;
  psc.domain.ilo[1] = 0;
  psc.domain.ilo[2] = 0;
  psc.domain.ihi[0] = 100;
  psc.domain.ihi[1] = 400;
  psc.domain.ihi[2] = 600;

  psc.domain.bnd_fld_lo[0] = BND_FLD_UPML;
  psc.domain.bnd_fld_hi[0] = BND_FLD_UPML;
  psc.domain.bnd_fld_lo[1] = BND_FLD_UPML;
  psc.domain.bnd_fld_hi[1] = BND_FLD_UPML;
  psc.domain.bnd_fld_lo[2] = BND_FLD_UPML;
  psc.domain.bnd_fld_hi[2] = BND_FLD_UPML;
  psc.domain.bnd_part[0] = BND_PART_REFLECTING;
  psc.domain.bnd_part[1] = BND_PART_REFLECTING;
  psc.domain.bnd_part[2] = BND_PART_REFLECTING;
}

static void
singlepart_init_field(void)
{
  // FIXME, do we need the ghost points?
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {

	// FIXME, why this time?
	FF3(EY, jx,jy,jz) = 0.0;

	// FIXME, this pulse needs a - to propagate in the right direction (+z)
	FF3(BX, jx,jy,jz) = 0.0;
      }
    }
  }
}

static void
singlepart_init_nvt(int kind, double x[3], double *q, double *m, double *n,
		    double v[3], double T[3])
{
  struct psc_singlepart *singlepart = psc.case_data;

  real Te = singlepart->Te, Ti = singlepart->Ti;

  real ld = psc.coeff.ld;

  real x0[3] = { singlepart->x0 / ld,
		 singlepart->y0 / ld,
		 singlepart->z0 / ld };

  real dens = 0.0;
  if ((psc.domain.ihi[0] - psc.domain.ilo[0] == 1 ||
       (int) (x[0] / psc.dx[0] + .5) == (int) (x0[0] / psc.dx[0] + .5)) && 
      (psc.domain.ihi[1] - psc.domain.ilo[1] == 1 ||
       (int) (x[1] / psc.dx[1] + .5) == (int) (x0[1] / psc.dx[1] + .5)) && 
      (psc.domain.ihi[2] - psc.domain.ilo[2] == 1 ||
       (int) (x[2] / psc.dx[2] + .5) == (int) (x0[2] / psc.dx[2] + .5))) {
    dens = 1.0;
  }

  switch (kind) {
  case 0: // electrons
    *q = -1.;
    *m = 1.;
    *n = dens;
    v[0] = 0.;
    v[1] = 0.;
    v[2] = 1.0e3;
    T[0] = Te;
    T[1] = Te;
    T[2] = Te;
    break;
  case 1: // ions
    *q = 1.;
    *m = singlepart->mass_ratio;
    *n = dens;
    v[0] = 0.;
    v[1] = 0.;
    v[2] = 0.;
    T[0] = Ti;
    T[1] = Ti;
    T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_ops_singlepart = {
  .name       = "singlepart",
  .create     = singlepart_create,
  .destroy    = singlepart_destroy,
  .init_param = singlepart_init_param,
  .init_field = singlepart_init_field,
  .init_nvt   = singlepart_init_nvt,
};
