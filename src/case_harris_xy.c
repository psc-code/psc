
#include "psc.h"
#include <mrc_common.h>
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// Plasma simulation parameters like Che et al '10
//
// c_over_vA : ratio of speed of light / Alfven velocity
// mass_ratio: ion mass / electron mass
// Te,Ti     : bulk temperature of electrons and ions (units m_i v_A^2)
// lambda    : B reversal length scale (units of d_i)
// xl[3]     : simulation box size (units of d_i)
// pert      : perturbation

// ----------------------------------------------------------------------
// case harris_xy
//
// same as harris, but with Bx = tanh(y) type of equilibrium

struct harris_xy {
  // parameters
  double c_over_vA;
  double Bguide;
  double Te, Ti;
  double mass_ratio;
  double lambda;
  double xl[3];
  double pert;

  // renormalized, such that c = 1, me = 1
  double _me, _mi;
  double _B0;
  double _Bguide;
  double _lambda;
  double _Te, _Ti;
  double _pert;
};

#define VAR(x) (void *)offsetof(struct harris_xy, x)

static struct param harris_xy_descr[] = {
  { "c_over_vA"     , VAR(c_over_vA)       , PARAM_DOUBLE(20.)    },
  { "Bguide"        , VAR(Bguide)          , PARAM_DOUBLE(5.)     },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(25.)    },
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(.04)    },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(.04)    },
  { "lambda"        , VAR(lambda)          , PARAM_DOUBLE(.5)     },
  { "lx"            , VAR(xl[0])           , PARAM_DOUBLE(4.)     },
  { "ly"            , VAR(xl[1])           , PARAM_DOUBLE(2.)     },
  { "lz"            , VAR(xl[2])           , PARAM_DOUBLE(8.)     },
  { "pert"          , VAR(pert)            , PARAM_DOUBLE(.1)     },
  {},
};

#undef VAR

static void
harris_xy_init_param(struct psc_case *Case)
{
  struct harris_xy *harris = Case->ctx;

  real c = 1.;
  real eps0 = 1.;
  real e = 1.;
  real mu0 = 1.;
  real n0 = 1.;
  real m_e = 1.;
  real m_i = harris->mass_ratio;

  psc.prm.qq = e;
  psc.prm.mm = m_e;
  psc.prm.tt = 1.;
  psc.prm.cc = c;
  psc.prm.eps0 = eps0;

  psc.prm.nmax = 16000;
  psc.prm.cpum = 5*24.0*60*60;
  psc.prm.lw = 2.*M_PI;
  psc.prm.i0 = 0.;
  psc.prm.n0 = 1.;
  psc.prm.e0 = 1.;

  psc.prm.nicell = 50;

  harris->_me = m_e;
  harris->_mi = m_i;
  // vA = B / \sqrt{\mu_0 m_i n0}
  harris->_B0 = c / harris->c_over_vA * sqrt(mu0 * m_i * n0);
  harris->_Bguide = harris->Bguide * harris->_B0;

  real d_e = sqrt(m_e);
  real d_i = sqrt(m_i);
  for (int d = 0; d < 3; d++) {
    psc.domain.length[d] = harris->xl[d] * d_i;
  }
  harris->_lambda = harris->lambda * d_i;

  // normalized temperatures are in terms of m_s c^2,
  // parameters in terms of m_i v_A ^ 2
  real v_A = harris->_B0 / sqrt(mu0 * m_i * n0);
  harris->_Ti = harris->Ti * m_i * sqr(v_A) / (m_i * sqr(c));
  harris->_Te = harris->Te * m_i * sqr(v_A) / (m_e * sqr(c));

  harris->_pert = harris->pert * harris->_B0 * d_e; // correct ?

  mpi_printf(MPI_COMM_WORLD, "::: B0     = %g\n", harris->_B0);
  mpi_printf(MPI_COMM_WORLD, "::: Bguide = %g\n", harris->_Bguide);
  mpi_printf(MPI_COMM_WORLD, "::: domain = %g d_e x %g d_e x %g d_e\n",
	     psc.domain.length[0], psc.domain.length[1], psc.domain.length[2]);
  mpi_printf(MPI_COMM_WORLD, "::: v_A    = %g c\n", v_A);
  mpi_printf(MPI_COMM_WORLD, "::: lambda = %g d_e\n", harris->_lambda);
  mpi_printf(MPI_COMM_WORLD, "::: om_pe  = %g\n", sqrt(n0 * sqr(e) / (m_e * eps0)));
  mpi_printf(MPI_COMM_WORLD, "::: om_pi  = %g\n", sqrt(n0 * sqr(e) / (m_i * eps0)));
  mpi_printf(MPI_COMM_WORLD, "::: om_ce  = %g\n", harris->_Bguide * e / m_e);
  mpi_printf(MPI_COMM_WORLD, "::: om_ci  = %g\n", harris->_Bguide * e / m_i);
  mpi_printf(MPI_COMM_WORLD, "::: om_cep  = %g\n", harris->_B0 * e / m_e);
  mpi_printf(MPI_COMM_WORLD, "::: om_cip  = %g\n", harris->_B0 * e / m_i);
  mpi_printf(MPI_COMM_WORLD, "::: d_e    = %g d_e\n", c / sqrt(n0 * sqr(e) / (m_e * eps0)));
  mpi_printf(MPI_COMM_WORLD, "::: d_i    = %g d_e\n", c / sqrt(n0 * sqr(e) / (m_i * eps0)));
  mpi_printf(MPI_COMM_WORLD, "::: l_Debye= %g d_e\n", sqrt(eps0 * harris->_Te / (n0 * sqr(e))));
  mpi_printf(MPI_COMM_WORLD, "\n");

  psc.domain.gdims[0] = 128;
  psc.domain.gdims[1] = 128;
  psc.domain.gdims[2] = 1;

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

static real
Bx0(struct harris_xy *harris, real y)
{
  double ly = psc.domain.length[1];
  return harris->_B0 * (-1. 
			+ tanh((y - 0.25*ly) / harris->_lambda)
			- tanh((y - 0.75*ly) / harris->_lambda)
			+ tanh((y - 1.25*ly) / harris->_lambda)
			- tanh((y - 1.75*ly) / harris->_lambda)
			+ tanh((y + 0.75*ly) / harris->_lambda)
			- tanh((y + 0.25*ly) / harris->_lambda));
}

static real
Az1(struct harris_xy *harris, real x, real y)
{
  double lx = psc.domain.length[0], ly = psc.domain.length[1];
  return harris->_pert / (2.*M_PI / lx) * sin(2.*M_PI * x / lx) * .5 * (1 - cos(2*2.*M_PI * y / ly));
}

static real
Bx(struct harris_xy *harris, real x, real y)
{
  real dy = 1e-3;
  return Bx0(harris, y) + 
    (Az1(harris, x, y + dy) - Az1(harris, x, y - dy)) / (2.*dy);
}

static real
By(struct harris_xy *harris, real x, real y)
{
  real dx = 1e-3;
  return - (Az1(harris, x + dx, y) - Az1(harris, x - dx, y)) / (2.*dx);
}

static real
Bz(struct harris_xy *harris, real x, real y)
{
  real B0 = harris->_B0;
  return sqrt(sqr(harris->_Bguide) + sqr(B0) - sqr(Bx0(harris, y)));
}

static real
Jx(struct harris_xy *harris, real x, real y)
{
  real dy = 1e-3;
  return (Bz(harris, x, y + dy) - Bz(harris, x, y - dy)) / (2.*dy);
}

static real
Jy(struct harris_xy *harris, real x, real y)
{
  real dx = 1e-3;
  return - (Bz(harris, x + dx, y) - Bz(harris, x - dx, y)) / (2.*dx);
}

static real
Jz(struct harris_xy *harris, real x, real y)
{
  real dx = 1e-3, dy = 1e-3;
  return
    (By(harris, x + dx, y) - By(harris, x - dx, y)) / (2.*dx) -
    (Bx(harris, x, y + dy) - Bx(harris, x, y - dy)) / (2.*dy);
  /* return - B0 / lambda * */
  /*   (1./sqr(cosh((y - 0.25*ly)/lambda)) - 1./sqr(cosh((y - .75*ly)/lambda))); */
}

static void
harris_xy_init_field(struct psc_case *Case)
{
  struct harris_xy *harris = Case->ctx;

  // FIXME, do we need the ghost points?
  foreach_patch(p) {
    fields_base_t *pf = &psc.flds.f[p];
    foreach_3d_g(p, jx, jy, jz) {
      double dx = psc.dx[0], dy = psc.dx[1];
      double xx = CRDX(p, jx), yy = CRDY(p, jy);
      
      F3_BASE(pf, HX, jx,jy,jz) = Bx(harris, xx, yy + .5*dy        );
      F3_BASE(pf, HY, jx,jy,jz) = By(harris, xx + .5*dx, yy        );
      F3_BASE(pf, HZ, jx,jy,jz) = Bz(harris, xx + .5*dx, yy + .5*dy);
      
      // FIXME centering
      F3_BASE(pf, JXI, jx,jy,jz) = Jx(harris, xx, yy);
      F3_BASE(pf, JYI, jx,jy,jz) = Jy(harris, xx, yy);
      F3_BASE(pf, JZI, jx,jy,jz) = Jz(harris, xx, yy);
    } foreach_3d_g_end;
  }
}

static void
harris_xy_init_npt(struct psc_case *Case, int kind, double x[3],
		   struct psc_particle_npt *npt)
{
  struct harris_xy *harris = Case->ctx;

  real Ti = harris->Ti, Te = harris->Te;

  double jx = Jx(harris, x[0], x[1]);
  double jy = Jy(harris, x[0], x[1]);
  double jz = Jz(harris, x[0], x[1]);
  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = harris->_me;
    npt->n = 1.;
    npt->p[0] = - Te / (Te + Ti) * jx / npt->n;
    npt->p[1] = - Te / (Te + Ti) * jy / npt->n;
    npt->p[2] = - Te / (Te + Ti) * jz / npt->n;
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = harris->_mi;
    npt->n = 1.;
    npt->p[0] = Ti / (Te + Ti) * jx / npt->n;
    npt->p[1] = Ti / (Te + Ti) * jy / npt->n;
    npt->p[2] = Ti / (Te + Ti) * jz / npt->n;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_ops_harris_xy = {
  .name       = "harris_xy",
  .ctx_size   = sizeof(struct harris_xy),
  .ctx_descr  = harris_xy_descr,
  .init_param = harris_xy_init_param,
  .init_field = harris_xy_init_field,
  .init_npt   = harris_xy_init_npt,
};

