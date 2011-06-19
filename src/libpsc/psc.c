
#include "psc.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_bnd.h"
#include "psc_collision.h"
#include "psc_randomize.h"
#include "psc_sort.h"
#include "psc_output_fields.h"
#include "psc_output_particles.h"
#include "psc_moments.h"
#include "psc_event_generator.h"
#include "psc_balance.h"

#include <mrc_common.h>
#include <mrc_params.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

struct psc *ppsc;

#define VAR(x) (void *)offsetof(struct psc, x)

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

struct param psc_descr[] = {
  // psc_domain
  { "length_x"      , VAR(domain.length[0])       , PARAM_DOUBLE(1e-6)   },
  { "length_y"      , VAR(domain.length[1])       , PARAM_DOUBLE(1e-6)   },
  { "length_z"      , VAR(domain.length[2])       , PARAM_DOUBLE(20e-6)  },
  { "corner_x"      , VAR(domain.corner[0])       , PARAM_DOUBLE(0.)     },
  { "corner_y"      , VAR(domain.corner[1])       , PARAM_DOUBLE(0.)     },
  { "corner_z"      , VAR(domain.corner[2])       , PARAM_DOUBLE(0.)     },
  { "gdims_x"       , VAR(domain.gdims[0])        , PARAM_INT(1)         },
  { "gdims_y"       , VAR(domain.gdims[1])        , PARAM_INT(1)         },
  { "gdims_z"       , VAR(domain.gdims[2])        , PARAM_INT(400)       },

  { "bnd_field_lo_x", VAR(domain.bnd_fld_lo[0])   , PARAM_SELECT(BND_FLD_PERIODIC,
								 bnd_fld_descr) },
  { "bnd_field_lo_y", VAR(domain.bnd_fld_lo[1])   , PARAM_SELECT(BND_FLD_PERIODIC,
								 bnd_fld_descr) },
  { "bnd_field_lo_z", VAR(domain.bnd_fld_lo[2])   , PARAM_SELECT(BND_FLD_PERIODIC,
								 bnd_fld_descr) },
  { "bnd_field_hi_x", VAR(domain.bnd_fld_hi[0])   , PARAM_SELECT(BND_FLD_PERIODIC,
								 bnd_fld_descr) },
  { "bnd_field_hi_y", VAR(domain.bnd_fld_hi[1])   , PARAM_SELECT(BND_FLD_PERIODIC,
								 bnd_fld_descr) },
  { "bnd_field_hi_z", VAR(domain.bnd_fld_hi[2])   , PARAM_SELECT(BND_FLD_PERIODIC,
								 bnd_fld_descr) },

  { "bnd_particle_x", VAR(domain.bnd_part[0])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },
  { "bnd_particle_y", VAR(domain.bnd_part[1])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },
  { "bnd_particle_z", VAR(domain.bnd_part[2])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },
  { "use_pml",        VAR(domain.use_pml)         , PARAM_BOOL(false)    },

  // psc_params
  { "qq"            , VAR(prm.qq)              , PARAM_DOUBLE(1.6021e-19)   },
  { "mm"            , VAR(prm.mm)              , PARAM_DOUBLE(9.1091e-31)   },
  { "tt"            , VAR(prm.tt)              , PARAM_DOUBLE(1.6021e-16)   },
  { "cc"            , VAR(prm.cc)              , PARAM_DOUBLE(3.0e8)        },
  { "eps0"          , VAR(prm.eps0)            , PARAM_DOUBLE(8.8542e-12)   },
  { "nmax"          , VAR(prm.nmax)            , PARAM_INT(0)               },
  { "cpum"          , VAR(prm.cpum)            , PARAM_DOUBLE(100000.)      },
  { "lw"            , VAR(prm.lw)              , PARAM_DOUBLE(3.2e-6)       },
  { "i0"            , VAR(prm.i0)              , PARAM_DOUBLE(1e21)         },
  { "n0"            , VAR(prm.n0)              , PARAM_DOUBLE(1e26)         },
  { "e0"            , VAR(prm.e0)              , PARAM_DOUBLE(0.)           },
  { "nicell"        , VAR(prm.nicell)          , PARAM_INT(200)             },
  // by default, we put the # of particles per cell according to the
  // density, using the weights (~ 1) only to fine-tune to the
  // right density.
  // if this parameter is set, we always use nicell particles / cell,
  // and adjust to the right density via the weights.
  { "const_num_particles_per_cell"
                    , VAR(prm.const_num_particles_per_cell), PARAM_BOOL(0)  },
  // a hack which allows to set the particle weight equal the density,
  // even though the # of particles was already set according to it.
  // this is for compatibility testing between Fortran and C initial
  // conditions, but I believe it's incorrect and should go away, eventually.
  { "fortran_particle_weight_hack"
                    , VAR(prm.fortran_particle_weight_hack), PARAM_BOOL(0)  },
  // yet another hackish thing for compatibility
  // adjust dt so that the laser period is an integer multiple of dt
  // only useful when actually doing lasers.
  { "adjust_dt_to_cycles"
                    , VAR(prm.adjust_dt_to_cycles), PARAM_BOOL(0)  },
  { "wallclock_limit"
                    , VAR(prm.wallclock_limit)    , PARAM_DOUBLE(0.) },
  { "from_checkpoint"
                    , VAR(prm.from_checkpoint)    , PARAM_BOOL(false) },

  {},
};

#undef VAR

// ----------------------------------------------------------------------
// psc_create

static void
_psc_create(struct psc *psc)
{
  MPI_Comm comm = psc_comm(psc);

  psc->push_particles = psc_push_particles_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->push_particles);
  psc->push_fields = psc_push_fields_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->push_fields);
  psc->bnd = psc_bnd_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->bnd);
  psc->collision = psc_collision_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->collision);
  psc->randomize = psc_randomize_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->randomize);
  psc->sort = psc_sort_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->sort);
  psc->output_fields = psc_output_fields_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->output_fields);
  psc->output_particles = psc_output_particles_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->output_particles);
  psc->moments = psc_moments_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->moments);
  psc->event_generator = psc_event_generator_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->event_generator);
  psc->balance = psc_balance_create(comm);
  psc_add_child(psc, (struct mrc_obj *) psc->balance);

  psc->time_start = MPI_Wtime();
}

// ======================================================================
// psc_setup

// ----------------------------------------------------------------------
// psc_setup_mrc_coeff

static void
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

// ----------------------------------------------------------------------
// psc_setup_mrc_domain

struct mrc_domain *
psc_setup_mrc_domain(struct psc *psc, int nr_patches)
{
  // FIXME, should be split to create, set_from_options, setup time?
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  // create a very simple domain decomposition
  int bc[3] = {};
  for (int d = 0; d < 3; d++) {
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_PERIODIC &&
	psc->domain.gdims[d] > 1) {
      bc[d] = BC_PERIODIC;
    }
  }

  mrc_domain_set_type(domain, "multi");
  mrc_domain_set_param_int3(domain, "m", psc->domain.gdims);
  mrc_domain_set_param_int(domain, "bcx", bc[0]);
  mrc_domain_set_param_int(domain, "bcy", bc[1]);
  mrc_domain_set_param_int(domain, "bcz", bc[2]);
  mrc_domain_set_param_int(domain, "nr_patches", nr_patches);

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "multi_uniform");
  mrc_crds_set_param_int(crds, "sw", 2);
  mrc_crds_set_param_float3(crds, "l",  (float[3]) { psc->domain.corner[0],
	psc->domain.corner[1], psc->domain.corner[2] });
  mrc_crds_set_param_float3(crds, "h",  (float[3]) {
      psc->domain.corner[0] + psc->domain.length[0],
      psc->domain.corner[1] + psc->domain.length[1],
      psc->domain.corner[2] + psc->domain.length[2] });

  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);

  // set up index bounds,
  // sanity checks for the decomposed domain
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &psc->nr_patches);
  psc->patch = calloc(psc->nr_patches, sizeof(*psc->patch));
  psc_foreach_patch(psc, p) {
    struct psc_patch *patch = &psc->patch[p];
    for (int d = 0; d < 3; d++) {
      patch->ldims[d] = patches[p].ldims[d];
      patch->off[d] = patches[p].off[d];
      patch->xb[d]  = patches[p].off[d] * psc->dx[d] + psc->domain.corner[d];
      
      int min_size = 1;
      if (patch->off[d] == 0 && // left-most patch in this dir
	  (psc->domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	   psc->domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
	min_size += psc->pml.size;
      }
      if (patch->off[d] + patch->ldims[d] == gdims[d] && // right-most patch in this dir
	  (psc->domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	   psc->domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
	min_size += psc->pml.size;
      }
      assert(psc->patch[p].ldims[d] >= min_size);
    }
  }

  return domain;
}

// ----------------------------------------------------------------------
// psc_setup_domain

static void
psc_setup_domain(struct psc *psc)
{
  struct psc_domain *domain = &psc->domain;

  for (int d = 0; d < 3; d++) {
    psc->ibn[d] = 3;
  }
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
  if (need_pml && !psc->domain.use_pml) {
    fprintf(stderr,
	    "WARNING: use_pml is disabled but pml boundary conditions requested.\n");
    fprintf(stderr,
	    "         I'm enabling use_pml.\n");
    psc->domain.use_pml = true;
  }
  psc->pml.thick = 10;
  psc->pml.cushion = psc->pml.thick / 3;
  psc->pml.size = psc->pml.thick + psc->pml.cushion;
  psc->pml.order = 3;

  psc->mrc_domain = psc_setup_mrc_domain(psc, -1);
}

// ----------------------------------------------------------------------
// _psc_setup

static void
_psc_setup(struct psc *psc)
{
  psc_setup_coeff(psc);
  psc_setup_domain(psc); // needs to be done before setting up psc_bnd
}

// ----------------------------------------------------------------------
// psc_destroy

static void
_psc_destroy(struct psc *psc)
{
  mfields_base_destroy(psc->flds);
  mparticles_base_destroy(&psc->particles);
  mphotons_destroy(&psc->mphotons);

  mrc_domain_destroy(psc->mrc_domain);
}

// ----------------------------------------------------------------------
// _psc_write

static void
_psc_write(struct psc *psc, struct mrc_io *io)
{
  mrc_domain_write(psc->mrc_domain, io);
}

// ======================================================================
// psc class

struct mrc_class_psc mrc_class_psc = {
  .name             = "psc",
  .size             = sizeof(struct psc),
  .param_descr      = psc_descr,
  .create           = _psc_create,
  .setup            = _psc_setup,
  .destroy          = _psc_destroy,
  .write            = _psc_write,
};

// ======================================================================

const char *fldname[NR_FIELDS] = {
  [NE]  = "ne",
  [NI]  = "ni",
  [NN]  = "nn",
  [JXI] = "jx",
  [JYI] = "jy",
  [JZI] = "jz",
  [EX]  = "ex",
  [EY]  = "ey",
  [EZ]  = "ez",
  [HX]  = "hx",
  [HY]  = "hy",
  [HZ]  = "hz",
  [DX]  = "dx",
  [DY]  = "dy",
  [DZ]  = "dz",
  [BX]  = "bx",
  [BY]  = "by",
  [BZ]  = "bz",
  [EPS] = "eps",
  [MU]  = "mu",
};
