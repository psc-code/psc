
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>

#include <mrc_params.h>

#include <math.h>

// originally by WFox

// Plasma simulation parameters
//   needed because they determine length scales for initial conditions
// B0:  peak Harris magnetic field  (magnitude gives ratio w_ce/w_pe)
// nb:  background particle density
// Te, Ti:  bulk temperature of electrons and ions (units m_e c^2)
// mi_over_me: ion mass / electron mass
// lambda = shear length scale (units of d_i)
// ly, lz = simulation box size (units of d_i)
// pert = perturbation (units of B * d_i)

// FIXME (description), below parameters don't include scaling factors

struct psc_harris {
  double B0;
  double nb;
  double Te, Ti;
  double mi_over_me;
  double lambda;
  double lz, ly;
  double pert;
};

#define to_psc_harris(psc) mrc_to_subobj(psc, struct psc_harris)

#define VAR(x) (void *)offsetof(struct psc_harris, x)
static struct param psc_harris_descr[] = {
  { "B0"            , VAR(B0)              , PARAM_DOUBLE(1.)     },
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(25.)    },
  { "nb"            , VAR(nb)              , PARAM_DOUBLE(.2)     },
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(1./12.) },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(5./12.) },
  { "lambda"        , VAR(lambda)          , PARAM_DOUBLE(.5)     },
  { "lz"            , VAR(lz)              , PARAM_DOUBLE(25.6)   },
  { "ly"            , VAR(ly)              , PARAM_DOUBLE(12.8)   },
  { "pert"          , VAR(pert)            , PARAM_DOUBLE(.1)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_harris_create

static void
psc_harris_create(struct psc *psc)
{
  struct psc_harris *harris = to_psc_harris(psc);

  // new defaults (dimensionless) for this case
  ppsc->prm.qq = 1.;
  ppsc->prm.mm = 1.;
  ppsc->prm.tt = 1.;
  ppsc->prm.cc = 1.;
  ppsc->prm.eps0 = 1.;

  ppsc->prm.nmax = 16000;
  ppsc->prm.cpum = 5*24.0*60*60;
  ppsc->prm.lw = 2.*M_PI;
  ppsc->prm.i0 = 0.;
  ppsc->prm.n0 = 1.;
  ppsc->prm.e0 = 1.;

  ppsc->prm.nicell = 50;

  ppsc->domain.length[0] = 1.; // no x dependence 
  ppsc->domain.length[1] = 2. * harris->ly; // double tearing
  ppsc->domain.length[2] = harris->lz;

  ppsc->domain.gdims[0] = 1;
  ppsc->domain.gdims[1] = 640;
  ppsc->domain.gdims[2] = 640;

  ppsc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_part[0] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part[1] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part[2] = BND_PART_PERIODIC;
}

// ----------------------------------------------------------------------
// psc_harris_init_field

static double
psc_harris_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_harris *harris = to_psc_harris(psc);

  double B0 = harris->B0;
  double lz = harris->lz, ly = harris->ly;
  double lambda = harris->lambda;
  double AA = harris->pert * B0;

  switch (m) {
  case HZ:
    return
      B0 * (-1. + tanh((x[1] - 0.5*ly) / lambda)- tanh((x[1] - 1.5*ly) / lambda))
      + AA * M_PI/ly * sin(2.*M_PI * x[2] / lz) * cos(M_PI * x[1] / ly);

  case HY:
    return - AA * 2.*M_PI / lz * cos(2.*M_PI * x[2] / lz) * sin(M_PI * x[1] / ly);

  case JXI:
    return B0 / lambda *
      (1./sqr(cosh((x[1] - 0.5*ly) / lambda)) - 1./sqr(cosh((x[1] - 1.5*ly) / lambda)))
      - (AA*sqr(M_PI) * (1./sqr(ly) + 4./sqr(lz)) 
	 * sin(2.*M_PI * x[2] / lz) * sin(M_PI * x[1] / ly));

  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_harris_init_npt

static void
psc_harris_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt)
{
  struct psc_harris *harris = to_psc_harris(psc);

  double B0 = harris->B0;
  double ly = harris->ly;
  double lambda = harris->lambda;
  double nb = harris->nb;
  double TTi = harris->Ti * sqr(B0) * harris->mi_over_me; // FIXME, why???
  double TTe = harris->Te * sqr(B0) * harris->mi_over_me;

  double jx0 = B0 / lambda *
    (1./sqr(cosh((x[1] - 0.5*ly) / lambda)) - 1./sqr(cosh((x[1] - 1.5*ly) / lambda)));

  npt->n = nb + 1./sqr(cosh((x[1] - 0.5*ly) / lambda)) + 1./sqr(cosh((x[1] - 1.5*ly) / lambda));
  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1. / harris->mi_over_me;
    npt->p[0] = - 2. * TTe / B0 / lambda * (jx0 / (B0 / lambda)) / npt->n;
    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 1.;
    npt->p[0] = 2. * TTi / B0 / lambda * (jx0 / B0 / lambda) / npt->n;
    npt->T[0] = TTi;
    npt->T[1] = TTi;
    npt->T[2] = TTi;
    break;
  default:
    assert(0);
  }
}

// ======================================================================
// psc_harris_ops

struct psc_ops psc_harris_ops = {
  .name             = "harris",
  .size             = sizeof(struct psc_harris),
  .param_descr      = psc_harris_descr,
  .create           = psc_harris_create,
  .init_field       = psc_harris_init_field,
  .init_npt         = psc_harris_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_harris_ops);
}
