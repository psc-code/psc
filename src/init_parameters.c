
#include "psc.h"
#include "util/params.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

struct psc_cmdline {
  const char *mod_particle;
};

#define VAR(x) (void *)offsetof(struct psc_domain, x)

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
  // FIXME, we should use something other than magic numbers here
  { "bnd_field_x"   , VAR(bnd_fld[0])      , PARAM_INT(1)         },
  { "bnd_field_y"   , VAR(bnd_fld[1])      , PARAM_INT(1)         },
  { "bnd_field_z"   , VAR(bnd_fld[2])      , PARAM_INT(1)         },
  { "bnd_particle_x", VAR(bnd_part[0])     , PARAM_INT(1)         },
  { "bnd_particle_y", VAR(bnd_part[1])     , PARAM_INT(1)         },
  { "bnd_particle_z", VAR(bnd_part[2])     , PARAM_INT(1)         },

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

#ifdef USE_PML
  fprintf(stderr, "PML not handled in C version!\n");
  MPI_Abort();
#endif

  for (int d = 0; d < 3; d++) {
    if (domain->ihi[d] - domain->ilo[d] == 1) {
      // if invariant in this direction, don't domain decompose
      assert(domain->nproc[d] == 1);
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
    psc.case_ops = psc_find_case(par.case_name);
    if (psc.case_ops->create) {
      psc.case_ops->create();
    }
    psc.case_ops->init_param();
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

  psc.domain.bnd_fld[0] = 1;
  psc.domain.bnd_fld[1] = 1;
  psc.domain.bnd_fld[2] = 1;
  psc.domain.bnd_part[0] = 1;
  psc.domain.bnd_part[1] = 1;
  psc.domain.bnd_part[2] = 1;
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
			  f_int *bnd_fld, f_int *bnd_part, f_int *nproc,
			  f_int *nghost);
void SET_param_domain_F77(f_real *length, f_int *itot, f_int *in, f_int *ix,
			  f_int *bnd_fld, f_int *bnd_part, f_int *nproc,
			  f_int *nghost);
void GET_param_psc_F77(f_real *qq, f_real *mm, f_real *tt, f_real *cc, f_real *eps0,
		       f_int *nmax, f_real *cpum, f_real *lw, f_real *i0, f_real *n0,
		       f_real *e0, f_int *nicell);
void SET_param_psc_F77(f_real *qq, f_real *mm, f_real *tt, f_real *cc, f_real *eps0,
		       f_int *nmax, f_real *cpum, f_real *lw, f_real *i0, f_real *n0,
		       f_real *e0, f_int *nicell);
void SET_param_coeff_F77(f_real *cori);

void
get_param_domain()
{
  struct psc_domain *p = &psc.domain;
  int imax[3];

  GET_param_domain_F77(p->length, p->itot, p->ilo, imax,
		       p->bnd_fld, p->bnd_part, p->nproc, p->nghost);
  for (int d = 0; d < 3; d++) {
    p->ihi[d] = imax[d] + 1;
  }
}

void
set_param_domain()
{
  struct psc_domain *p = &psc.domain;
  int imax[3];

  for (int d = 0; d < 3; d++) {
    imax[d] = p->ihi[d] - 1;
  }
  SET_param_domain_F77(p->length, p->itot, p->ilo, imax,
		       p->bnd_fld, p->bnd_part, p->nproc, p->nghost);
}

void
get_param_psc()
{
  struct psc_param *p = &psc.prm;
  GET_param_psc_F77(&p->qq, &p->mm, &p->tt, &p->cc, &p->eps0,
		    &p->nmax, &p->cpum, &p->lw, &p->i0, &p->n0, &p->e0,
		    &p->nicell);
}

void
set_param_psc()
{
  struct psc_param *p = &psc.prm;
  SET_param_psc_F77(&p->qq, &p->mm, &p->tt, &p->cc, &p->eps0,
		    &p->nmax, &p->cpum, &p->lw, &p->i0, &p->n0, &p->e0,
		    &p->nicell);
}

void
set_param_coeff()
{
  struct psc_coeff *p = &psc.coeff;
  SET_param_coeff_F77(&p->cori);
}

void
INIT_param_domain()
{
  INIT_param_domain_F77();
  get_param_domain();
}

void
INIT_param_psc()
{
  INIT_param_psc_F77();
  get_param_psc();
}

void
C_init_param_F77()
{
  init_param_domain_default();
  init_param_psc_default();
  init_case();
  init_param_domain();
  init_param_psc();
  init_param_coeff();

  set_param_domain();
  set_param_psc();
  set_param_coeff();
}

