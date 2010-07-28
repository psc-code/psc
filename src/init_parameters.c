
#include "psc.h"
#include "util/params.h"

#include <math.h>
#include <assert.h>

struct psc_cmdline {
  const char *mod_particle;
};

#define VAR(x) (void *)offsetof(struct psc_domain, x)

static struct param_select bnd_fld_descr[] = {
  { .val = BND_FLD_OPEN        , .str = "open"        },
  { .val = BND_FLD_PERIODIC    , .str = "periodic"    },
  { .val = BND_FLD_UPML        , .str = "upml"        },
  { .val = BND_FLD_TIME        , .str = "time"        },
  {},
};

static struct param_select bnd_part_descr[] = {
  { .val = BND_PART_REFLECTING , .str = "reflecting"  },
  { .val = BND_PART_PERIODIC   , .str = "periodic"    },
  {},
};

static struct param psc_domain_descr[] = {
  { "length_x"      , VAR(length[0])       , PARAM_DOUBLE(1e-6)   },
  { "length_y"      , VAR(length[1])       , PARAM_DOUBLE(1e-6)   },
  { "length_z"      , VAR(length[2])       , PARAM_DOUBLE(20e-6)  },
  { "itot_x"        , VAR(itot[0])         , PARAM_INT(10)        },
  { "itot_y"        , VAR(itot[1])         , PARAM_INT(10)        },
  { "itot_z"        , VAR(itot[2])         , PARAM_INT(400)       },
  { "ilo_x"         , VAR(ilo[0])          , PARAM_INT(8)         },
  { "ilo_y"         , VAR(ilo[1])          , PARAM_INT(8)         },
  { "ilo_z"         , VAR(ilo[2])          , PARAM_INT(0)         },
  { "ihi_x"         , VAR(ihi[0])          , PARAM_INT(9)         },
  { "ihi_y"         , VAR(ihi[1])          , PARAM_INT(9)         },
  { "ihi_z"         , VAR(ihi[2])          , PARAM_INT(400)       },

  { "bnd_field_lo_x", VAR(bnd_fld_lo[0])   , PARAM_SELECT(BND_FLD_PERIODIC,
							  bnd_fld_descr) },
  { "bnd_field_lo_y", VAR(bnd_fld_lo[1])   , PARAM_SELECT(BND_FLD_PERIODIC,
							  bnd_fld_descr) },
  { "bnd_field_lo_z", VAR(bnd_fld_lo[2])   , PARAM_SELECT(BND_FLD_PERIODIC,
							  bnd_fld_descr) },
  { "bnd_field_hi_x", VAR(bnd_fld_hi[0])   , PARAM_SELECT(BND_FLD_PERIODIC,
							  bnd_fld_descr) },
  { "bnd_field_hi_y", VAR(bnd_fld_hi[1])   , PARAM_SELECT(BND_FLD_PERIODIC,
							  bnd_fld_descr) },
  { "bnd_field_hi_z", VAR(bnd_fld_hi[2])   , PARAM_SELECT(BND_FLD_PERIODIC,
							  bnd_fld_descr) },

  { "bnd_particle_x", VAR(bnd_part[0])     , PARAM_SELECT(BND_PART_PERIODIC,
							  bnd_part_descr) },
  { "bnd_particle_y", VAR(bnd_part[1])     , PARAM_SELECT(BND_PART_PERIODIC,
							  bnd_part_descr) },
  { "bnd_particle_z", VAR(bnd_part[2])     , PARAM_SELECT(BND_PART_PERIODIC,
							  bnd_part_descr) },
  { "nproc_x",        VAR(nproc[0])        , PARAM_INT(1)         },
  { "nproc_y",        VAR(nproc[1])        , PARAM_INT(1)         },
  { "nproc_z",        VAR(nproc[2])        , PARAM_INT(1)         },
  {},
};

#undef VAR

void
init_param_domain()
{
  struct psc_domain *domain = &psc.domain;

  params_parse_cmdline_nodefault(domain, psc_domain_descr, "PSC domain",
				 MPI_COMM_WORLD);
  params_print(domain, psc_domain_descr, "PSC domain", MPI_COMM_WORLD);

  psc.pml.thick = 10;
  psc.pml.cushion = psc.pml.thick / 3;
  psc.pml.size = psc.pml.thick + psc.pml.cushion;
  psc.pml.order = 3;

  SET_param_pml();

  for (int d = 0; d < 3; d++) {
    if (domain->ihi[d] - domain->ilo[d] == 1) {
      // if invariant in this direction:
      // can't domain decompose in, set bnd to periodic
      assert(domain->nproc[d] == 1);
      domain->bnd_fld_lo[d] = BND_FLD_PERIODIC;
      domain->bnd_fld_hi[d] = BND_FLD_PERIODIC;
      domain->bnd_part[d]   = BND_PART_PERIODIC;
    } else {
      // on every proc, need domain at least nghost wide
      assert((domain->ihi[d] - domain->ilo[d]) >= domain->nproc[d] * domain->nghost[d]);
    }
  }
}

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
  {},
};

#undef VAR

void
init_param_psc()
{
  params_parse_cmdline_nodefault(&psc.prm, psc_param_descr, "PSC parameters",
				 MPI_COMM_WORLD);
  params_print(&psc.prm, psc_param_descr, "PSC parameters", MPI_COMM_WORLD);
}

// ----------------------------------------------------------------------
// set up cases

static struct psc_case_ops *psc_case_ops_list[] = {
  &psc_case_ops_langmuir,
  &psc_case_ops_wakefield,
  &psc_case_ops_thinfoil,
  &psc_case_ops_singlepart,
  &psc_case_ops_harris,
  NULL,
};

static struct psc_case_ops *
psc_find_case(const char *name)
{
  for (int i = 0; psc_case_ops_list[i]; i++) {
    if (strcasecmp(psc_case_ops_list[i]->name, name) == 0)
      return psc_case_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_case '%s' not available.\n", name);
  abort();
}

struct psc_case *
psc_case_create(const char *case_name)
{
  struct psc_case *Case = malloc(sizeof(*Case));
  Case->ops = psc_find_case(case_name);

  size_t ctx_size = Case->ops->ctx_size;
  if (ctx_size) {
    Case->ctx = malloc(ctx_size);
    memset(Case->ctx, 0, ctx_size);
  }

  struct param *ctx_descr = Case->ops->ctx_descr;
  if (ctx_descr) {
    char cn[strlen(case_name) + 6];
    sprintf(cn, "case %s", case_name);
    params_parse_cmdline(Case->ctx, ctx_descr, cn, MPI_COMM_WORLD);
    params_print(Case->ctx, ctx_descr, cn, MPI_COMM_WORLD);
  }

  if (Case->ops->create) {
    Case->ops->create(Case);
  }
  
  return Case;
}

void
psc_case_destroy(struct psc_case *Case)
{
  if (Case->ops->destroy) {
    Case->ops->destroy(Case);
  }

  free(Case->ctx);
  free(Case);
}

struct par_case {
  char *case_name;
};

#define VAR(x) (void *)offsetof(struct par_case, x)

static struct param par_case_descr[] = {
  { "case"          , VAR(case_name)        , PARAM_STRING(NULL)   },
  {},
};

#undef VAR

void
init_case()
{
  struct par_case par;
  params_parse_cmdline(&par, par_case_descr, "PSC case", MPI_COMM_WORLD);
  params_print(&par, par_case_descr, "PSC case", MPI_COMM_WORLD);

  if (par.case_name) {
    psc.Case = psc_case_create(par.case_name);
  }

  if (psc.Case) {
    psc_case_init_param(psc.Case);
  }
}

static void
init_param_domain_default()
{
  INIT_param_domain();

#if 0 // if we ever want to throw out Fortran entirely...
  psc.domain.nghost[0] = 3;
  psc.domain.nghost[1] = 3;
  psc.domain.nghost[2] = 3;

  psc.domain.length[0] = 1.  * 1e-6;
  psc.domain.length[1] = 1.  * 1e-6;
  psc.domain.length[2] = 20. * 1e-6;

  psc.domain.itot[0] = 10;
  psc.domain.itot[1] = 10;
  psc.domain.itot[2] = 400;

  psc.domain.nproc[0] = 1;
  psc.domain.nproc[1] = 1;
  psc.domain.nproc[2] = 1;

  psc.domain.ilo[0] = 8;
  psc.domain.ihi[0] = 9;
  psc.domain.ilo[1] = 8;
  psc.domain.ihi[1] = 9;
  psc.domain.ilo[2] = 0;
  psc.domain.ihi[2] = 400;

  psc.domain.bnd_fld[0] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld[1] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld[2] = BND_FLD_PERIODIC;
  psc.domain.bnd_part[0] = BND_PART_PERIODIC;
  psc.domain.bnd_part[1] = BND_PART_PERIODIC;
  psc.domain.bnd_part[2] = BND_PART_PERIODIC;
#endif
}

static void
init_param_psc_default()
{
  INIT_param_psc();

#if 0
  psc.prm.qq = 1.6021e-19;
  psc.prm.mm = 9.1091e-31;
  psc.prm.tt = 1.6021e-16;
  psc.prm.cc = 3.0e8;
  psc.prm.eps0 = 8.8542e-12;

  psc.prm.nmax = 0;
  psc.prm.cpum = 100000.0;
  psc.prm.lw = 3.2*1.0e-6;
  psc.prm.i0 = 1.0e21;
  psc.prm.n0 = 1.0e26;
  psc.prm.e0 = 0;

  psc.prm.nicell = 200;
#endif
}

static void
init_param_coeff()
{
  assert(psc.prm.nicell > 0);
  psc.coeff.cori = 1. / psc.prm.nicell;
  psc.coeff.wl = 2. * M_PI * psc.prm.cc / psc.prm.lw;
  psc.coeff.ld = psc.prm.cc / psc.coeff.wl;
  if (psc.prm.e0 == 0.) {
    psc.prm.e0 = sqrt(2.0 * psc.prm.i0 / psc.prm.eps0 / psc.prm.cc) /
      psc.prm.lw / 1.0e6;
  }
  psc.prm.b0 = psc.prm.e0 / psc.prm.cc;
  psc.prm.rho0 = psc.prm.eps0 * psc.coeff.wl * psc.prm.b0;
  psc.prm.phi0 = psc.coeff.ld * psc.prm.e0;
  psc.prm.a0 = psc.prm.e0 / psc.coeff.wl;
  psc.coeff.vos = psc.prm.qq * psc.prm.e0 / (psc.prm.mm * psc.coeff.wl);
  psc.coeff.vt = sqrt(psc.prm.tt / psc.prm.mm);
  psc.coeff.wp = sqrt(sqr(psc.prm.qq) * psc.prm.n0 / psc.prm.eps0 / psc.prm.mm);
  psc.coeff.alpha = psc.coeff.wp / psc.coeff.wl;
  psc.coeff.beta = psc.coeff.vt / psc.prm.cc;
  psc.coeff.eta = psc.coeff.vos / psc.prm.cc;

  for (int d = 0; d < 3; d++) {
    psc.dx[d] = psc.domain.length[d] / psc.coeff.ld / psc.domain.itot[d];
  }
  psc.dt = .75 * sqrt(1./(1./sqr(psc.dx[0]) + 1./sqr(psc.dx[1]) + 1./sqr(psc.dx[2])));

  // adjust to match laser cycles FIXME, this isn't a good place,
  // and hardcoded params (2, 30.)
  psc.coeff.nnp = ceil(2.*M_PI / psc.dt);
  if (psc.prm.adjust_dt_to_cycles) {
    psc.dt = 2.*M_PI / psc.coeff.nnp;
  }
  psc.coeff.np = 2 * psc.coeff.nnp;
  if (psc.prm.nmax == 0) {
    psc.prm.nmax = 30. * psc.coeff.nnp;
  }

}

void
psc_init_param()
{
  init_param_domain_default();
  init_param_psc_default();
  init_case();
  init_param_domain();
  init_param_psc();
  init_param_coeff();
}

// ======================================================================
// Fortran glue

#define INIT_param_domain_F77   F77_FUNC_(init_param_domain, INIT_PARAM_DOMAIN)
#define INIT_param_psc_F77      F77_FUNC_(init_param_psc, INIT_PARAM_PSC)
#define GET_param_domain_F77    F77_FUNC_(get_param_domain, GET_PARAM_DOMAIN)
#define SET_param_domain_F77    F77_FUNC_(set_param_domain, SET_PARAM_DOMAIN)
#define GET_param_psc_F77       F77_FUNC_(get_param_psc, GET_PARAM_PSC)
#define SET_param_psc_F77       F77_FUNC_(set_param_psc, SET_PARAM_PSC)
#define SET_param_coeff_F77     F77_FUNC_(set_param_coeff, SET_PARAM_COEFF)
#define C_init_param_F77        F77_FUNC_(c_init_param, C_INIT_PARAM)

void INIT_param_domain_F77(void);
void INIT_param_psc_F77(void);
void GET_param_domain_F77(f_real *length, f_int *itot, f_int *in, f_int *ix,
			  f_int *bnd_fld_lo, f_int *bnd_fld_hi, f_int *bnd_part,
			  f_int *nproc, f_int *nghost);
void SET_param_domain_F77(f_real *length, f_int *itot, f_int *in, f_int *ix,
			  f_int *bnd_fld_lo, f_int *bnd_fld_hi, f_int *bnd_part,
			  f_int *nproc, f_int *nghost);
void GET_param_psc_F77(f_real *qq, f_real *mm, f_real *tt, f_real *cc, f_real *eps0,
		       f_int *nmax, f_real *cpum, f_real *lw, f_real *i0, f_real *n0,
		       f_real *e0, f_real *b0, f_real *j0, f_real *rho0, f_real *phi0,
		       f_real *a0, f_int *nicell);
void SET_param_psc_F77(f_real *qq, f_real *mm, f_real *tt, f_real *cc, f_real *eps0,
		       f_int *nmax, f_real *cpum, f_real *lw, f_real *i0, f_real *n0,
		       f_real *e0, f_real *b0, f_real *j0, f_real *rho0, f_real *phi0,
		       f_real *a0, f_int *nicell);
void SET_param_coeff_F77(f_real *cori, f_real *alpha, f_real *beta, f_real *eta,
			 f_real *wl, f_real *ld, f_real *vos, f_real *vt, f_real *wp,
			 f_real *dx, f_real *dt, f_int *np, f_int *nnp);

void
GET_param_domain()
{
  struct psc_domain *p = &psc.domain;
  int imax[3];

  GET_param_domain_F77(p->length, p->itot, p->ilo, imax,
		       p->bnd_fld_lo, p->bnd_fld_hi, p->bnd_part, p->nproc, p->nghost);
  for (int d = 0; d < 3; d++) {
    p->ihi[d] = imax[d] + 1;
  }
}

void
SET_param_domain()
{
  struct psc_domain *p = &psc.domain;
  int imax[3];

  for (int d = 0; d < 3; d++) {
    imax[d] = p->ihi[d] - 1;
  }
  SET_param_domain_F77(p->length, p->itot, p->ilo, imax, p->bnd_fld_lo, p->bnd_fld_hi,
		       p->bnd_part, p->nproc, p->nghost);
}

void
GET_param_psc()
{
  struct psc_param *p = &psc.prm;
  GET_param_psc_F77(&p->qq, &p->mm, &p->tt, &p->cc, &p->eps0,
		    &p->nmax, &p->cpum, &p->lw, &p->i0, &p->n0, &p->e0, &p->b0,
		    &p->j0, &p->rho0, &p->phi0, &p->a0,
		    &p->nicell);
}

void
SET_param_psc()
{
  struct psc_param *p = &psc.prm;
  SET_param_psc_F77(&p->qq, &p->mm, &p->tt, &p->cc, &p->eps0,
		    &p->nmax, &p->cpum, &p->lw, &p->i0, &p->n0, &p->e0, &p->b0,
		    &p->j0, &p->rho0, &p->phi0, &p->a0,
		    &p->nicell);
}

void
SET_param_coeff()
{
  struct psc_coeff *p = &psc.coeff;
  SET_param_coeff_F77(&p->cori, &p->alpha, &p->beta, &p->eta,
		      &p->wl, &p->ld, &p->vos, &p->vt, &p->wp,
		      psc.dx, &psc.dt, &p->np, &p->nnp);
}

void
INIT_param_domain()
{
  INIT_param_domain_F77();
  GET_param_domain();
}

void
INIT_param_psc()
{
  INIT_param_psc_F77();
  GET_param_psc();
}

void
C_init_param_F77()
{
  psc_init_param();

  SET_param_domain();
  SET_param_psc();
  SET_param_coeff();
}
