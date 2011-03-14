
#include "psc.h"
#include <mrc_params.h>

#include <math.h>
#include <assert.h>

struct psc_cmdline {
  const char *mod_particle;
};

#define VAR(x) (void *)offsetof(struct psc_domain, x)

static struct mrc_param_select bnd_fld_descr[] = {
  { .val = BND_FLD_OPEN        , .str = "open"        },
  { .val = BND_FLD_PERIODIC    , .str = "periodic"    },
  { .val = BND_FLD_UPML        , .str = "upml"        },
  { .val = BND_FLD_TIME        , .str = "time"        },
  {},
};

static struct mrc_param_select bnd_part_descr[] = {
  { .val = BND_PART_REFLECTING , .str = "reflecting"  },
  { .val = BND_PART_PERIODIC   , .str = "periodic"    },
  {},
};

static struct param psc_domain_descr[] = {
  { "length_x"      , VAR(length[0])       , PARAM_DOUBLE(1e-6)   },
  { "length_y"      , VAR(length[1])       , PARAM_DOUBLE(1e-6)   },
  { "length_z"      , VAR(length[2])       , PARAM_DOUBLE(20e-6)  },
  { "gdims_x"       , VAR(gdims[0])        , PARAM_INT(1)         },
  { "gdims_y"       , VAR(gdims[1])        , PARAM_INT(1)         },
  { "gdims_z"       , VAR(gdims[2])        , PARAM_INT(400)       },

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
  { "use_pml",        VAR(use_pml)         , PARAM_BOOL(false)    },
  {},
};

#undef VAR

void
psc_set_default_domain()
{
  mrc_params_set_default(&psc.domain, psc_domain_descr);
  for (int d = 0; d < 3; d++) {
    psc.ibn[d] = 3;
  }
}

void
psc_set_from_options_domain()
{
  struct psc_domain *domain = &psc.domain;

  mrc_params_parse_nodefault(domain, psc_domain_descr, "PSC domain",
			     MPI_COMM_WORLD);
  mrc_params_print(domain, psc_domain_descr, "PSC domain", MPI_COMM_WORLD);

  psc.pml.thick = 10;
  psc.pml.cushion = psc.pml.thick / 3;
  psc.pml.size = psc.pml.thick + psc.pml.cushion;
  psc.pml.order = 3;
}

void
psc_setup_domain()
{
  struct psc_domain *domain = &psc.domain;

  SET_param_pml();

  bool need_pml = false;
  for (int d = 0; d < 3; d++) {
    if (domain->gdims[d] == 1) {
      // if invariant in this direction:
      // set bnd to periodic (FIXME?)
      domain->bnd_fld_lo[d] = BND_FLD_PERIODIC;
      domain->bnd_fld_hi[d] = BND_FLD_PERIODIC;
      domain->bnd_part[d]   = BND_PART_PERIODIC;
    } else {
      if (domain->bnd_fld_lo[d] >= BND_FLD_UPML ||
	  domain->bnd_fld_hi[d] >= BND_FLD_UPML) {
	need_pml = true;
      }
    }
  }
  if (need_pml && !psc.domain.use_pml) {
    fprintf(stderr,
	    "WARNING: use_pml is disabled but pml boundary conditions requested.\n");
    fprintf(stderr,
	    "         I'm enabling use_pml.\n");
    psc.domain.use_pml = true;
  }

  // FIXME, should be split to create, set_from_options, setup time?
  psc.mrc_domain = mrc_domain_create(MPI_COMM_WORLD);
  // create a very simple domain decomposition
  int bc[3] = {};
  for (int d = 0; d < 3; d++) {
    if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC &&
	psc.domain.gdims[d] > 1) {
      bc[d] = BC_PERIODIC;
    }
  }

  mrc_domain_set_type(psc.mrc_domain, "multi");
  mrc_domain_set_param_int3(psc.mrc_domain, "m", psc.domain.gdims);
  mrc_domain_set_param_int(psc.mrc_domain, "bcx", bc[0]);
  mrc_domain_set_param_int(psc.mrc_domain, "bcy", bc[1]);
  mrc_domain_set_param_int(psc.mrc_domain, "bcz", bc[2]);

  struct mrc_crds *crds = mrc_domain_get_crds(psc.mrc_domain);
  mrc_crds_set_type(crds, "multi_uniform");
  mrc_crds_set_param_int(crds, "sw", 2);
  mrc_crds_set_param_float3(crds, "h",  (float[3]) { psc.domain.length[0],
	psc.domain.length[1], psc.domain.length[2] });

  mrc_domain_set_from_options(psc.mrc_domain);
  mrc_domain_setup(psc.mrc_domain);

  // set up index bounds,
  // sanity checks for the decomposed domain
  int gdims[3];
  mrc_domain_get_global_dims(psc.mrc_domain, gdims);
  struct mrc_patch *patches = mrc_domain_get_patches(psc.mrc_domain, &psc.nr_patches);
  psc.patch = calloc(psc.nr_patches, sizeof(*psc.patch));
  foreach_patch(p) {
    struct psc_patch *patch = &psc.patch[p];
    for (int d = 0; d < 3; d++) {
      patch->ldims[d] = patches[p].ldims[d];
      patch->off[d] = patches[p].off[d];
      patch->xb[d]  = patches[p].off[d] * psc.dx[d];
      
      int min_size = 1;
      if (patch->off[d] == 0 && // left-most patch in this dir
	  (psc.domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	   psc.domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
	min_size += psc.pml.size;
      }
      if (patch->off[d] + patch->ldims[d] == gdims[d] && // right-most patch in this dir
	  (psc.domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	   psc.domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
	min_size += psc.pml.size;
      }
      assert(psc.patch[p].ldims[d] >= min_size);
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
  { "wallclock_limit"
                    , VAR(wallclock_limit)    , PARAM_DOUBLE(0.) },
  { "from_checkpoint"
                    , VAR(from_checkpoint)    , PARAM_BOOL(false) },
  {},
};

#undef VAR

void
psc_set_default_psc()
{
  mrc_params_set_default(&psc.prm, psc_param_descr);
}

void
psc_set_from_options_psc()
{
  mrc_params_parse_nodefault(&psc.prm, psc_param_descr, "PSC parameters",
			     MPI_COMM_WORLD);
  mrc_params_print(&psc.prm, psc_param_descr, "PSC parameters", MPI_COMM_WORLD);
}

// ----------------------------------------------------------------------
// set up cases

static struct psc_case_ops *psc_case_ops_list[] = {
  &psc_case_ops_langmuir,
  &psc_case_ops_wakefield,
  &psc_case_ops_thinfoil,
  &psc_case_ops_foils,
  &psc_case_ops_curvedfoil,
  &psc_case_ops_singlepart,
  &psc_case_ops_harris,
  &psc_case_ops_harris_xy,
  &psc_case_ops_collisions,
  &psc_case_ops_test_xz,
  &psc_case_ops_test_yz,
  &psc_case_ops_cone,	
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

void
psc_setup_coeff()
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
    if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC){
      psc.dx[d] = psc.domain.length[d] / psc.coeff.ld / psc.domain.gdims[d];
    } else {
      psc.dx[d] = psc.domain.length[d] / psc.coeff.ld / (psc.domain.gdims[d] - 1);
    }
  }
  psc.dt = .75 * sqrt(1./(1./sqr(psc.dx[0]) + 1./sqr(psc.dx[1]) + 1./sqr(psc.dx[2])));
  mpi_printf(MPI_COMM_WORLD, "::: dt      = %g\n", psc.dt);
  mpi_printf(MPI_COMM_WORLD, "::: dx      = %g %g %g\n", psc.dx[0], psc.dx[1], psc.dx[2]);

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

