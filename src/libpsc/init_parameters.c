
#include "psc.h"
#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <assert.h>

struct psc_cmdline {
  const char *mod_particle;
};

// ======================================================================

#define VAR(x) (void *)offsetof(struct psc_param, x)

static struct param psc_param_descr[] = {
  { "qq"            , VAR(qq)              , PARAM_DOUBLE(1.6021e-19)   },
  { "mm"            , VAR(mm)              , PARAM_DOUBLE(9.1091e-31)   },
  { "tt"            , VAR(tt)              , PARAM_DOUBLE(1.6021e-16)   },
  { "cc"            , VAR(cc)              , PARAM_DOUBLE(3.0e8)        },
  { "eps0"          , VAR(eps0)            , PARAM_DOUBLE(8.8542e-12)   },
  { "nmax"          , VAR(nmax)            , PARAM_INT(0)               },
  { "cpum"          , VAR(cpum)            , PARAM_DOUBLE(100000.)      },
  { "lw"            , VAR(lw)              , PARAM_DOUBLE(3.2e-6)       },
  { "i0"            , VAR(i0)              , PARAM_DOUBLE(1e21)         },
  { "n0"            , VAR(n0)              , PARAM_DOUBLE(1e26)         },
  { "e0"            , VAR(e0)              , PARAM_DOUBLE(0.)           },
  { "nicell"        , VAR(nicell)          , PARAM_INT(200)             },
  // by default, we put the # of particles per cell according to the
  // density, using the weights (~ 1) only to fine-tune to the
  // right density.
  // if this parameter is set, we always use nicell particles / cell,
  // and adjust to the right density via the weights.
  { "const_num_particles_per_cell"
                    , VAR(const_num_particles_per_cell), PARAM_BOOL(0)  },
  // a hack which allows to set the particle weight equal the density,
  // even though the # of particles was already set according to it.
  // this is for compatibility testing between Fortran and C initial
  // conditions, but I believe it's incorrect and should go away, eventually.
  { "fortran_particle_weight_hack"
                    , VAR(fortran_particle_weight_hack), PARAM_BOOL(0)  },
  // yet another hackish thing for compatibility
  // adjust dt so that the laser period is an integer multiple of dt
  // only useful when actually doing lasers.
  { "adjust_dt_to_cycles"
                    , VAR(adjust_dt_to_cycles), PARAM_BOOL(0)  },
  { "wallclock_limit"
                    , VAR(wallclock_limit)    , PARAM_DOUBLE(0.) },
  { "from_checkpoint"
                    , VAR(from_checkpoint)    , PARAM_BOOL(false) },
  {},
};

#undef VAR

void
psc_set_default_psc(struct psc *psc)
{
  mrc_params_set_default(&psc->prm, psc_param_descr);
}

void
psc_set_from_options_psc(struct psc *psc)
{
  mrc_params_parse_nodefault(&psc->prm, psc_param_descr, "PSC parameters",
			     MPI_COMM_WORLD);
}

void
psc_view_psc(struct psc *psc)
{
  mrc_params_print(&psc->prm, psc_param_descr, "PSC parameters", MPI_COMM_WORLD);
}

// ======================================================================

void
psc_setup_coeff(struct psc *psc)
{
  assert(psc->prm.nicell > 0);
  psc->coeff.cori = 1. / psc->prm.nicell;
  psc->coeff.wl = 2. * M_PI * psc->prm.cc / psc->prm.lw;
  psc->coeff.ld = psc->prm.cc / psc->coeff.wl;
  if (psc->prm.e0 == 0.) {
    psc->prm.e0 = sqrt(2.0 * psc->prm.i0 / psc->prm.eps0 / psc->prm.cc) /
      psc->prm.lw / 1.0e6;
  }
  psc->prm.b0 = psc->prm.e0 / psc->prm.cc;
  psc->prm.rho0 = psc->prm.eps0 * psc->coeff.wl * psc->prm.b0;
  psc->prm.phi0 = psc->coeff.ld * psc->prm.e0;
  psc->prm.a0 = psc->prm.e0 / psc->coeff.wl;
  psc->coeff.vos = psc->prm.qq * psc->prm.e0 / (psc->prm.mm * psc->coeff.wl);
  psc->coeff.vt = sqrt(psc->prm.tt / psc->prm.mm);
  psc->coeff.wp = sqrt(sqr(psc->prm.qq) * psc->prm.n0 / psc->prm.eps0 / psc->prm.mm);
  psc->coeff.alpha = psc->coeff.wp / psc->coeff.wl;
  psc->coeff.beta = psc->coeff.vt / psc->prm.cc;
  psc->coeff.eta = psc->coeff.vos / psc->prm.cc;

  for (int d = 0; d < 3; d++) {
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_PERIODIC){
      psc->dx[d] = psc->domain.length[d] / psc->coeff.ld / psc->domain.gdims[d];
    } else {
      psc->dx[d] = psc->domain.length[d] / psc->coeff.ld / (psc->domain.gdims[d] - 1);
    }
  }
  psc->dt = .75 * sqrt(1./(1./sqr(psc->dx[0]) + 1./sqr(psc->dx[1]) + 1./sqr(psc->dx[2])));
  mpi_printf(MPI_COMM_WORLD, "::: dt      = %g\n", psc->dt);
  mpi_printf(MPI_COMM_WORLD, "::: dx      = %g %g %g\n", psc->dx[0], psc->dx[1], psc->dx[2]);

  // adjust to match laser cycles FIXME, this isn't a good place,
  // and hardcoded params (2, 30.)
  psc->coeff.nnp = ceil(2.*M_PI / psc->dt);
  if (psc->prm.adjust_dt_to_cycles) {
    psc->dt = 2.*M_PI / psc->coeff.nnp;
  }
  psc->coeff.np = 2 * psc->coeff.nnp;
  if (psc->prm.nmax == 0) {
    psc->prm.nmax = 30. * psc->coeff.nnp;
  }
}

