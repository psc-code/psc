
#include "psc.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_bnd.h"
#include "psc_bnd_particles.h"
#include "psc_bnd_photons.h"
#include "psc_bnd_fields.h"
#include "psc_collision.h"
#include "psc_randomize.h"
#include "psc_sort.h"
#include "psc_marder.h"
#include "psc_diag.h"
#include "psc_output_fields_collection.h"
#include "psc_output_particles.h"
#include "psc_output_photons.h"
#include "psc_event_generator.h"
#include "psc_balance.h"
#include "psc_checks.h"
#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include <mrc_common.h>
#include <mrc_params.h>
#include <mrc_io.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <time.h>

struct psc *ppsc;

#define VAR(x) (void *)offsetof(struct psc, x)

static struct mrc_param_select bnd_fld_descr[] = {
  { .val = BND_FLD_OPEN           , .str = "open"            },
  { .val = BND_FLD_PERIODIC       , .str = "periodic"        },
  { .val = BND_FLD_UPML           , .str = "upml"            },
  { .val = BND_FLD_TIME           , .str = "time"            },
  { .val = BND_FLD_CONDUCTING_WALL, .str = "conducting_wall" },
  {},
};

static struct mrc_param_select bnd_part_descr[] = {
  { .val = BND_PART_REFLECTING , .str = "reflecting"  },
  { .val = BND_PART_PERIODIC   , .str = "periodic"    },
  { .val = BND_PART_ABSORBING  , .str = "absorbing"    },
  { .val = BND_PART_OPEN  , .str = "open"    },
  {},
};

static struct param psc_descr[] = {
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
  { "np_x"	    , VAR(domain.np[0])	          , PARAM_INT(1)	 },
  { "np_y"	    , VAR(domain.np[1])	          , PARAM_INT(1)	 },
  { "np_z"	    , VAR(domain.np[2])	          , PARAM_INT(1)	 },

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

  { "bnd_particle_lo_x", VAR(domain.bnd_part_lo[0])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },
  { "bnd_particle_lo_y", VAR(domain.bnd_part_lo[1])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },
  { "bnd_particle_lo_z", VAR(domain.bnd_part_lo[2])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },

  { "bnd_particle_hi_x", VAR(domain.bnd_part_hi[0])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },
  { "bnd_particle_hi_y", VAR(domain.bnd_part_hi[1])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },
  { "bnd_particle_hi_z", VAR(domain.bnd_part_hi[2])     , PARAM_SELECT(BND_PART_PERIODIC,
								 bnd_part_descr) },

  { "use_pml",        VAR(domain.use_pml)         , PARAM_BOOL(false)    },

  // psc_params
  { "qq"            , VAR(prm.qq)              , PARAM_DOUBLE(1.6021e-19)   },
  { "mm"            , VAR(prm.mm)              , PARAM_DOUBLE(9.1091e-31)   },
  { "tt"            , VAR(prm.tt)              , PARAM_DOUBLE(1.6021e-16)   },
  { "cc"            , VAR(prm.cc)              , PARAM_DOUBLE(3.0e8)        },
  { "eps0"          , VAR(prm.eps0)            , PARAM_DOUBLE(8.8542e-12)   },
  { "nmax"          , VAR(prm.nmax)            , PARAM_INT(0)               },
  { "lw"            , VAR(prm.lw)              , PARAM_DOUBLE(3.2e-6)       },
  { "i0"            , VAR(prm.i0)              , PARAM_DOUBLE(1e21)         },
  { "n0"            , VAR(prm.n0)              , PARAM_DOUBLE(1e26)         },
  { "e0"            , VAR(prm.e0)              , PARAM_DOUBLE(0.)           },
  { "cfl"           , VAR(prm.cfl)             , PARAM_DOUBLE(.75)          },
  { "nicell"        , VAR(prm.nicell)          , PARAM_INT(200)             },
  { "nr_populations", VAR(prm.nr_populations)  , PARAM_INT(-1),
    .help = "number of particle populations in the initial condition. "
    "init_npt() will be called this many times. By default, nr_populations "
    "will be set the the number of particle kinds." },
  { "neutralizing_population", VAR(prm.neutralizing_population)  , PARAM_INT(-1),
    .help = "this population will get density set to achieve neutrality "
    "in a given cell." },
  { "seed_by_time"  , VAR(prm.seed_by_time)    , PARAM_BOOL(false)          },
  // by default, we put the # of particles per cell according to the
  // density, using the weights (~ 1) only to fine-tune to the
  // right density.
  // if this parameter is set, we always use nicell particles / cell,
  // and adjust to the right density via the weights.
  { "fractional_n_particles_per_cell"
                    , VAR(prm.fractional_n_particles_per_cell), PARAM_BOOL(0)  },
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
  { "initial_momentum_gamma_correction"
                    , VAR(prm.initial_momentum_gamma_correction), PARAM_BOOL(0),
    .help = "if set, interpret momenta as velocities and multiply by gamma to get "
    "relativistic momenta." },
  
  { "adjust_dt_to_cycles"
                    , VAR(prm.adjust_dt_to_cycles), PARAM_BOOL(0)  },
  { "gdims_in_terms_of_cells"
                    , VAR(prm.gdims_in_terms_of_cells), PARAM_BOOL(false),
    .help = "if set, the given global dimensions (gdims) will be used for "
    "setting the number of cells, rather than the number of nodes in the grid." }, 
  { "initial_particle_shift"
                    , VAR(prm.initial_particle_shift), PARAM_DOUBLE(0.) },
  { "wallclock_limit"
                    , VAR(prm.wallclock_limit)    , PARAM_DOUBLE(0.) },
  { "write_checkpoint"
                    , VAR(prm.write_checkpoint)   , PARAM_BOOL(false) },
  { "write_checkpoint_every_step"
                    , VAR(prm.write_checkpoint_every_step), PARAM_INT(-1) },

  { "fields_base"   , VAR(prm.fields_base)        , PARAM_STRING("c") },
  { "particles_base", VAR(prm.particles_base)     , PARAM_STRING("c") },
  { "particles_base_flags"
                    , VAR(prm.particles_base_flags)  , PARAM_INT(0) },
  { "stats_every"
                    , VAR(prm.stats_every)        , PARAM_INT(1),
    .help = "sets every how many steps we log timing and other stats." },
  { "detailed_profiling"
                    , VAR(prm.detailed_profiling) , PARAM_BOOL(false),
    .help = "output profiling information by MPI process rather than aggregated." },
  { "theta_xz"      , VAR(prm.theta_xz)           , PARAM_DOUBLE(0.),
    .help = "rotate initial particle shifted Maxwellian in x-z plane." },

  { "push_particles"          , VAR(push_particles)          , MRC_VAR_OBJ(psc_push_particles) },
  { "push_fields"             , VAR(push_fields)             , MRC_VAR_OBJ(psc_push_fields) },
  { "bnd"                     , VAR(bnd)                     , MRC_VAR_OBJ(psc_bnd) },
  { "bnd_particles"           , VAR(bnd_particles)           , MRC_VAR_OBJ(psc_bnd_particles) },
  { "bnd_photons"             , VAR(bnd_photons)             , MRC_VAR_OBJ(psc_bnd_photons) },
  { "collision"               , VAR(collision)               , MRC_VAR_OBJ(psc_collision) },
  { "randomize"               , VAR(randomize)               , MRC_VAR_OBJ(psc_randomize) },
  { "marder"                  , VAR(marder)                  , MRC_VAR_OBJ(psc_marder) },
  { "sort"                    , VAR(sort)                    , MRC_VAR_OBJ(psc_sort) },
  { "diag"                    , VAR(diag)                    , MRC_VAR_OBJ(psc_diag) },
  { "output_fields_collection", VAR(output_fields_collection), MRC_VAR_OBJ(psc_output_fields_collection) },
  { "output_particles"        , VAR(output_particles)        , MRC_VAR_OBJ(psc_output_particles) },
  { "output_photons"          , VAR(output_photons)          , MRC_VAR_OBJ(psc_output_photons) },
  { "event_generator"         , VAR(event_generator)         , MRC_VAR_OBJ(psc_event_generator) },
  { "balance"                 , VAR(balance)                 , MRC_VAR_OBJ(psc_balance) },
  { "checks"                  , VAR(checks)                  , MRC_VAR_OBJ(psc_checks) },

  {},
};

#undef VAR

// ----------------------------------------------------------------------
// psc_create

static void
_psc_create(struct psc *psc)
{
  assert(!ppsc);
  ppsc = psc;

  // default: 2 species (e-, i+)
  psc_set_kinds(psc, 2, NULL);

  psc_bnd_set_psc(psc->bnd, psc); // FIXME, do general parent interface?
  psc_bnd_particles_set_psc(psc->bnd_particles, psc);
  psc_bnd_photons_set_psc(psc->bnd_photons, psc);
  psc_output_fields_collection_set_psc(psc->output_fields_collection, psc);

  psc->time_start = MPI_Wtime();
}

// ----------------------------------------------------------------------
// psc_set_from_options

static void
_psc_set_from_options(struct psc *psc)
{
  if (psc->use_dynamic_patches) {
    psc_patchmanager_set_from_options(&psc->patchmanager);
  }

  // make comma separated list of current kinds
  char *s = calloc(100, 1);
  char *s_save = s;
  if (psc->nr_kinds > 0) {
    strcpy(s, psc->kinds[0].name);
    for (int k = 1; k < psc->nr_kinds; k++) {
      strcat(s, ",");
      strcat(s, psc->kinds[k].name);
    }
  }
  char **ps = &s;
  // allow user to change names, or even change number
  mrc_params_get_option_string_help("particle_kinds", (const char **) ps,
				    "names of particle kinds, separated by commas");
  // parse comma separated list back into names
  char *p, *ss;
  int k = 0;
  ss = strdup(*ps);
  char *ss_save;
  ss_save = ss;
  
  while ((p = strsep(&ss, ", "))) {
    k++;
  }
  free(ss);
  free(ss_save);
  psc_set_kinds(psc, k, NULL);

  k = 0;
  ss = strdup(*ps);
  ss_save = ss;
  while ((p = strsep(&ss, ", "))) {
    free(psc->kinds[k].name);
    psc->kinds[k].name = strdup(p);
    k++;
  }
  free(ss);
  free(ss_save);
  s = s_save;
    
  // allow setting of parameters for each kind
  for (int k = 0; k < psc->nr_kinds; k++) {
    struct psc_kind *kind = psc->kinds + k;
    assert(kind->name);
    sprintf(s, "particle_%s_q", kind->name);
    mrc_params_get_option_double_help(s, &kind->q, 
				      "charge of this particle kind");
    sprintf(s, "particle_%s_m", kind->name);
    mrc_params_get_option_double_help(s, &kind->m, 
				      "mass of this particle kind");
    sprintf(s, "particle_%s_n", kind->name);
    mrc_params_get_option_double_help(s, &kind->n, 
				      "default density of this particle kind");
    sprintf(s, "particle_%s_T", kind->name);
    mrc_params_get_option_double_help(s, &kind->T, 
				      "default temperature of this particle kind");
  }
  free(s);
}

// ======================================================================
// psc_setup

// ----------------------------------------------------------------------
// psc_setup_coeff

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
}

// ----------------------------------------------------------------------
// psc_setup_mrc_domain

struct mrc_domain *
psc_setup_mrc_domain(struct psc *psc, int nr_patches)
{
  if (psc_ops(psc)->setup_mrc_domain) {
    return psc_ops(psc)->setup_mrc_domain(psc, nr_patches);
  }

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
  if(psc->use_dynamic_patches) {
    psc_patchmanager_timestep(&psc->patchmanager);
    mrc_domain_set_param_ptr(domain, "activepatches", psc->patchmanager.activepatches);
  }
  mrc_domain_set_param_int3(domain, "m", psc->domain.gdims);
  mrc_domain_set_param_int(domain, "bcx", bc[0]);
  mrc_domain_set_param_int(domain, "bcy", bc[1]);
  mrc_domain_set_param_int(domain, "bcz", bc[2]);
  mrc_domain_set_param_int(domain, "nr_patches", nr_patches);
  mrc_domain_set_param_int3(domain, "np", psc->domain.np);

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", 2);
  mrc_crds_set_param_double3(crds, "l",  (double[3]) { psc->domain.corner[0],
	psc->domain.corner[1], psc->domain.corner[2] });
  mrc_crds_set_param_double3(crds, "h",  (double[3]) {
      psc->domain.corner[0] + psc->domain.length[0],
      psc->domain.corner[1] + psc->domain.length[1],
      psc->domain.corner[2] + psc->domain.length[2] });

  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);

  return domain;
}

// ----------------------------------------------------------------------
// psc_setup_patches

void
psc_setup_patches(struct psc *psc, struct mrc_domain *domain)
{
  double dx[3];
  for (int d = 0; d < 3; d++) {
    // FIXME, we should eventually settle on one or the other way of handling
    // non-periodic boundaries -- use gdims to count the number of cells in the
    // grid, or count the number of nodes (corners)
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_PERIODIC ||
	psc->domain.gdims[d] == 1 ||
	psc->prm.gdims_in_terms_of_cells) {
      dx[d] = psc->domain.length[d] / psc->coeff.ld / psc->domain.gdims[d];
    } else {
      dx[d] = psc->domain.length[d] / psc->coeff.ld / (psc->domain.gdims[d] - 1);
    }
  }

  double inv_sum = 0.;
  int nr_levels;
  mrc_domain_get_nr_levels(psc->mrc_domain, &nr_levels);
  for (int d=0;d<3;d++) {
    if (psc->domain.gdims[d] > 1) {
      inv_sum += 1. / sqr(dx[d] / (1 << (nr_levels - 1)));
    }
  }
  assert(inv_sum);  //Simulation has 0 dimensions
  psc->dt = psc->prm.cfl * sqrt(1./inv_sum);

  mpi_printf(MPI_COMM_WORLD, "::: dt      = %g\n", psc->dt);
  mpi_printf(MPI_COMM_WORLD, "::: dx      = %g %g %g\n", dx[0], dx[1], dx[2]);

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

  // set up index bounds,
  // sanity checks for the decomposed domain
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &psc->nr_patches);
  psc->patch = calloc(psc->nr_patches, sizeof(*psc->patch));
  psc_foreach_patch(psc, p) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(domain, p, &info);
    struct psc_patch *patch = &psc->patch[p];
    for (int d = 0; d < 3; d++) {
      patch->ldims[d] = patches[p].ldims[d];
      patch->off[d] = patches[p].off[d];
      patch->dx[d]  = dx[d] / (1 << info.level);
      patch->xb[d]  = patches[p].off[d] * patch->dx[d] + psc->domain.corner[d] / psc->coeff.ld;
      
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
}

// ----------------------------------------------------------------------
// psc_setup_domain

void
psc_setup_domain(struct psc *psc)
{
  struct psc_domain *domain = &psc->domain;

  for (int d = 0; d < 3; d++) {
    // FIXME, old-style particle pushers need 3 ghost points still
    psc->ibn[d] = 2;
  }
  bool need_pml = false;
  for (int d = 0; d < 3; d++) {
    if (domain->gdims[d] == 1) {
      // if invariant in this direction:
      // set bnd to periodic (FIXME?)
      domain->bnd_fld_lo[d] = BND_FLD_PERIODIC;
      domain->bnd_fld_hi[d] = BND_FLD_PERIODIC;
      domain->bnd_part_lo[d]   = BND_PART_PERIODIC;
      domain->bnd_part_hi[d]   = BND_PART_PERIODIC;
      // and no ghost points
      psc->ibn[d] = 0;
    } else {
      if ((domain->bnd_fld_lo[d] >= BND_FLD_UPML && domain->bnd_fld_lo[d] <= BND_FLD_TIME) ||
	  (domain->bnd_fld_hi[d] >= BND_FLD_UPML && domain->bnd_fld_hi[d] <= BND_FLD_TIME)) {
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

  //Needs to be done before setting up the domain
  if (psc->use_dynamic_patches) {
    psc_patchmanager_setup(&psc->patchmanager, psc->domain.np);
  }

  if (!psc->mrc_domain) {
    psc->mrc_domain = psc_setup_mrc_domain(psc, -1);
  }
  psc_setup_patches(psc, psc->mrc_domain);
}

// ----------------------------------------------------------------------
// psc_setup_partition_and_particles

void
psc_setup_partition_and_particles(struct psc *psc)
{
  // alloc / initialize particles
  int particle_label_offset;
  int *nr_particles_by_patch = calloc(psc->nr_patches, sizeof(*nr_particles_by_patch));
  psc_setup_partition(psc, nr_particles_by_patch, &particle_label_offset);
  psc_balance_initial(psc->balance, psc, &nr_particles_by_patch);

  psc->particles = 
    psc_mparticles_create(mrc_domain_comm(psc->mrc_domain));
  psc_mparticles_set_type(psc->particles, psc->prm.particles_base);
  psc_mparticles_set_name(psc->particles, "mparticles");
  psc_mparticles_set_domain(psc->particles, psc->mrc_domain);
  if (psc->prm.particles_base_flags == 0) {
    psc->prm.particles_base_flags = psc_push_particles_get_mp_flags(ppsc->push_particles);
  }
  psc_mparticles_set_param_int(psc->particles, "flags", psc->prm.particles_base_flags);
  psc_mparticles_setup(psc->particles);
  psc_mparticles_reserve(psc->particles, nr_particles_by_patch);

  psc_setup_particles(psc, nr_particles_by_patch, particle_label_offset);
  free(nr_particles_by_patch);
}

// ----------------------------------------------------------------------
// psc_setup_photons

static void
psc_setup_photons(struct psc *psc)
{
  psc->mphotons = psc_mphotons_create(psc_comm(psc));
  psc_mphotons_set_name(psc->mphotons, "mphotons");
  psc_mphotons_set_domain(psc->mphotons, psc->mrc_domain);
  psc_mphotons_setup(psc->mphotons);

  if (!psc_ops(psc)->init_photon_np) {
    // if photons aren't initialized, we'll just have zero of them
    return;
  }

  psc_foreach_patch(psc, p) {
    photons_t *photons = &psc->mphotons->p[p];

    int np = 0;
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
      struct psc_photon_np photon_np = {}; // init to all zero
      psc_ops(psc)->init_photon_np(psc, xx, &photon_np);
      np += photon_np.n_in_cell;
    } psc_foreach_3d_end;

    photons_alloc(photons, np);

    int i = 0;
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
      struct psc_photon_np photon_np = {}; // init to all zero
      psc_ops(psc)->init_photon_np(psc, xx, &photon_np);
	    
      for (int cnt = 0; cnt < photon_np.n_in_cell; cnt++) {
	photon_t *p = photons_get_one(photons, i++);
	      
	float ran1, ran2, ran3, ran4, ran5, ran6;
	do {
	  ran1 = random() / ((float) RAND_MAX + 1);
	  ran2 = random() / ((float) RAND_MAX + 1);
	  ran3 = random() / ((float) RAND_MAX + 1);
	  ran4 = random() / ((float) RAND_MAX + 1);
	  ran5 = random() / ((float) RAND_MAX + 1);
	  ran6 = random() / ((float) RAND_MAX + 1);
	} while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f ||
		 ran4 >= 1.f || ran5 >= 1.f || ran6 >= 1.f);
	
	float kx =
	  sqrtf(-2.f*photon_np.sigma_k[0]*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2)
	  + photon_np.k[0];
	float ky =
	  sqrtf(-2.f*photon_np.sigma_k[1]*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4)
	  + photon_np.k[1];
	float kz =
	  sqrtf(-2.f*photon_np.sigma_k[2]*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6)
	  + photon_np.k[2];

	p->x[0] = xx[0];
	p->x[1] = xx[1];
	p->x[2] = xx[2];
	p->p[0] = kx;
	p->p[1] = ky;
	p->p[2] = kz;
	p->wni  = photon_np.n / photon_np.n_in_cell;
      }
    } psc_foreach_3d_end;

    photons->nr = i;
  }
}

// ----------------------------------------------------------------------
// _psc_setup

static void
_psc_setup(struct psc *psc)
{
  if (!psc_ops(psc)) { // old-style: setup is handled by psc_case
    psc_setup_member_objs(psc);
    return;
  }

  psc_setup_coeff(psc);
  psc_setup_domain(psc);

  psc_setup_partition_and_particles(psc);
  psc_setup_fields(psc);
  psc_setup_photons(psc);
#ifdef USE_FORTRAN
  psc_setup_fortran(psc);
#endif

  psc_setup_member_objs(psc);
}

// ----------------------------------------------------------------------
// psc_destroy

static void
_psc_destroy(struct psc *psc)
{
  psc_mfields_list_del(&psc_mfields_base_list, &psc->flds);
  psc_mfields_destroy(psc->flds);
  psc_mparticles_destroy(psc->particles);
  psc_mphotons_destroy(psc->mphotons);

  mrc_domain_destroy(psc->mrc_domain);
  free(psc->patch);

  if (psc->use_dynamic_patches) {
    psc_patchmanager_destroy(&psc->patchmanager);
  }

  if (psc->kinds) {
    for (int k = 0; k < psc->nr_kinds; k++) {
      free(psc->kinds[k].name);
    }
    free(psc->kinds);
  }
  
  ppsc = NULL;
}

// ----------------------------------------------------------------------
// _psc_write

static void
_psc_write(struct psc *psc, struct mrc_io *io)
{
  mrc_io_write_int(io, psc, "timestep", psc->timestep);
  mrc_io_write_int(io, psc, "nr_kinds", psc->nr_kinds);

  for (int k = 0; k < psc->nr_kinds; k++) {
    char s[20];
    sprintf(s, "kind_q%d", k);
    mrc_io_write_double(io, psc, s, psc->kinds[k].q);
    sprintf(s, "kind_m%d", k);
    mrc_io_write_double(io, psc, s, psc->kinds[k].m);
    sprintf(s, "kind_n%d", k);
    mrc_io_write_double(io, psc, s, psc->kinds[k].n);
    sprintf(s, "kind_T%d", k);
    mrc_io_write_double(io, psc, s, psc->kinds[k].T);
    sprintf(s, "kind_name%d", k);
    mrc_io_write_string(io, psc, s, psc->kinds[k].name);
  }

  mrc_io_write_ref(io, psc, "mrc_domain", psc->mrc_domain);
  mrc_io_write_ref(io, psc, "mparticles", psc->particles);
  mrc_io_write_ref(io, psc, "mfields", psc->flds);
  mrc_io_write_ref(io, psc, "mphotons", psc->mphotons);
}

// ----------------------------------------------------------------------
// _psc_read

static void
_psc_read(struct psc *psc, struct mrc_io *io)
{
  assert(!ppsc);
  ppsc = psc;

  psc_setup_coeff(psc);

  mrc_io_read_int(io, psc, "timestep", &psc->timestep);
  mrc_io_read_int(io, psc, "nr_kinds", &psc->nr_kinds);

  psc->kinds = calloc(psc->nr_kinds, sizeof(*psc->kinds));
  for (int k = 0; k < psc->nr_kinds; k++) {
    char s[20];
    sprintf(s, "kind_q%d", k);
    mrc_io_read_double(io, psc, s, &psc->kinds[k].q);
    sprintf(s, "kind_m%d", k);
    mrc_io_read_double(io, psc, s, &psc->kinds[k].m);
    sprintf(s, "kind_n%d", k);
    mrc_io_read_double(io, psc, s, &psc->kinds[k].n);
    sprintf(s, "kind_T%d", k);
    mrc_io_read_double(io, psc, s, &psc->kinds[k].T);
    sprintf(s, "kind_name%d", k);
    mrc_io_read_string(io, psc, s, &psc->kinds[k].name);
  }

  psc->mrc_domain = mrc_io_read_ref(io, psc, "mrc_domain", mrc_domain);
  psc_setup_domain(psc);
#ifdef USE_FORTRAN
  psc_setup_fortran(psc);
#endif

  psc->particles = mrc_io_read_ref(io, psc, "mparticles", psc_mparticles);
  psc->flds = mrc_io_read_ref(io, psc, "mfields", psc_mfields);
  psc_mfields_list_add(&psc_mfields_base_list, &psc->flds);
  psc->mphotons = mrc_io_read_ref(io, psc, "mphotons", psc_mphotons);

  psc_read_member_objs(psc, io);

  psc->time_start = MPI_Wtime();
}

// ----------------------------------------------------------------------
// get_n_in_cell
//
// helper function for partition / particle setup

static inline int
get_n_in_cell(struct psc *psc, struct psc_particle_npt *npt)
{
  if (psc->prm.const_num_particles_per_cell) {
    return psc->prm.nicell;
  }
  if (npt->particles_per_cell) {
    return npt->n * npt->particles_per_cell + .5;
  }
  if (psc->prm.fractional_n_particles_per_cell) {
    int n_prts = npt->n / psc->coeff.cori;
    float rmndr = npt->n / psc->coeff.cori - n_prts;
    float ran = random() / ((float) RAND_MAX + 1);
    if (ran < rmndr) {
      n_prts++;
    }
    return n_prts;
  }
  return npt->n / psc->coeff.cori + .5;
}

// ----------------------------------------------------------------------
// pml_find_bounds
//
// helper function for partition / particle setup

// particles must not be placed in the pml regions
// check if pml bnds are set and restrict avoid particle placement inside pml regions

static void
pml_find_bounds(struct psc *psc, int p, int ilo[3], int ihi[3])
{
  struct psc_patch *patch = &psc->patch[p];
  for (int d = 0; d < 3; d++) {
    ilo[d] = 0;
    if (patch->off[d] == 0 && // left-most proc in this dir
	(psc->domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	 psc->domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
      ilo[d] += psc->pml.size+1;
    }
    ihi[d] = patch->ldims[d];
    if (ihi[d] + patch->off[d] == psc->domain.gdims[d] && // right-most proc in this dir
	(psc->domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	 psc->domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
      ihi[d] -= psc->pml.size+1;
    }
  }
}

// ----------------------------------------------------------------------
// psc_setup_partition

void
psc_setup_partition(struct psc *psc, int *nr_particles_by_patch,
		    int *particle_label_offset)
{
  if (psc_ops(psc)->setup_particles) {
    psc_ops(psc)->setup_particles(psc, nr_particles_by_patch, true);
    return;
  }
  if (!psc_ops(psc)->init_npt) {
    psc_foreach_patch(psc, p) {
      nr_particles_by_patch[p] = 0;
    }
    return;
  }

  if (psc->prm.nr_populations < 0) {
    psc->prm.nr_populations = psc->nr_kinds;
  }
  if (psc->prm.neutralizing_population < 0) {
    psc->prm.neutralizing_population = psc->prm.nr_populations - 1;
  }

  int np_total = 0;
  psc_foreach_patch(psc, p) {
    int ilo[3], ihi[3];
    pml_find_bounds(psc, p, ilo, ihi);

    int np = 0;
    int nr_pop = psc->prm.nr_populations;
    for (int kind = 0; kind < nr_pop; kind++) {
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	    double xx[3] = { .5 * (CRDX(p, jx) + CRDX(p, jx+1)),
			     .5 * (CRDY(p, jy) + CRDY(p, jy+1)),
			     .5 * (CRDZ(p, jz) + CRDZ(p, jz+1)) };
	    // FIXME, the issue really is that (2nd order) particle pushers
	    // don't handle the invariant dim right
	    if (psc->domain.gdims[0] == 1) xx[0] = CRDX(p, jx);
	    if (psc->domain.gdims[1] == 1) xx[1] = CRDY(p, jy);
	    if (psc->domain.gdims[2] == 1) xx[2] = CRDZ(p, jz);

	    struct psc_particle_npt npt = {};
	    if (kind < psc->nr_kinds) {
	      npt.kind = kind;
	      npt.q    = psc->kinds[kind].q;
	      npt.m    = psc->kinds[kind].m;
	      npt.n    = psc->kinds[kind].n;
	      npt.T[0] = psc->kinds[kind].T;
	      npt.T[1] = psc->kinds[kind].T;
	      npt.T[2] = psc->kinds[kind].T;
	    };
	    psc_ops(psc)->init_npt(psc, kind, xx, &npt);

	    int n_in_cell = get_n_in_cell(psc, &npt);
	    if (psc->prm.fractional_n_particles_per_cell) {
	      n_in_cell++; // we may get an extra particle
	    }
	    np += n_in_cell;
	  }
	}
      }
    }
    nr_particles_by_patch[p] = np;
    np_total += np;
  }

  // calculate global particle label offset for unique numbering
  // FIXME, won't work with fractional particles
  *particle_label_offset = 0; // necessary on proc 0
  MPI_Exscan(&np_total, particle_label_offset, 1, MPI_INT, MPI_SUM,
	     MPI_COMM_WORLD);
}

void
psc_setup_particle(struct psc *psc, particle_t *prt, struct psc_particle_npt *npt,
		   int p, double xx[3])
{
  double beta = psc->coeff.beta;

  float ran1, ran2, ran3, ran4, ran5, ran6;
  do {
    ran1 = random() / ((float) RAND_MAX + 1);
    ran2 = random() / ((float) RAND_MAX + 1);
    ran3 = random() / ((float) RAND_MAX + 1);
    ran4 = random() / ((float) RAND_MAX + 1);
    ran5 = random() / ((float) RAND_MAX + 1);
    ran6 = random() / ((float) RAND_MAX + 1);
  } while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f ||
	   ran4 >= 1.f || ran5 >= 1.f || ran6 >= 1.f);
	      
  double pxi = npt->p[0] +
    sqrtf(-2.f*npt->T[0]/npt->m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
  double pyi = npt->p[1] +
    sqrtf(-2.f*npt->T[1]/npt->m*sqr(beta)*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4);
  double pzi = npt->p[2] +
    sqrtf(-2.f*npt->T[2]/npt->m*sqr(beta)*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6);

  if (psc->prm.initial_momentum_gamma_correction) {
    double gam;
    if (sqr(pxi) + sqr(pyi) + sqr(pzi) < 1.) {
      gam = 1. / sqrt(1. - sqr(pxi) - sqr(pyi) - sqr(pzi));
      pxi *= gam;
      pyi *= gam;
      pzi *= gam;
    }
  }
  
  assert(npt->kind >= 0 && npt->kind < psc->nr_kinds);
  prt->kind = npt->kind;
  assert(npt->q == psc->kinds[prt->kind].q);
  assert(npt->m == psc->kinds[prt->kind].m);
  /* prt->qni = psc->kinds[prt->kind].q; */
  /* prt->mni = psc->kinds[prt->kind].m; */
  prt->xi = xx[0] - psc->patch[p].xb[0];
  prt->yi = xx[1] - psc->patch[p].xb[1];
  prt->zi = xx[2] - psc->patch[p].xb[2];
  prt->pxi = pxi * cos(psc->prm.theta_xz) + pzi * sin(psc->prm.theta_xz);
  prt->pyi = pyi;
  prt->pzi = - pxi * sin(psc->prm.theta_xz) + pzi * cos(psc->prm.theta_xz);
}	      

// ----------------------------------------------------------------------
// psc_setup_particles

void
psc_setup_particles(struct psc *psc, int *nr_particles_by_patch,
		    int particle_label_offset)
{
  if (psc_ops(psc)->setup_particles) {
    psc_ops(psc)->setup_particles(psc, nr_particles_by_patch, false);
    return;
  }
  if (!psc_ops(psc)->init_npt)
    return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // FIXME, why does this do anything to the random seed?
  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, PARTICLE_TYPE, MP_DONT_COPY);

  if (psc->prm.seed_by_time) {
    srandom(10*rank + time(NULL));
  } else {
    srandom(rank);
  }

  psc_foreach_patch(psc, p) {
    int ilo[3], ihi[3];
    pml_find_bounds(psc, p, ilo, ihi);
    particle_range_t prts = particle_range_mprts(mprts, p);
  
    int i = 0;
    int nr_pop = psc->prm.nr_populations;
    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	  double xx[3] = { .5 * (CRDX(p, jx) + CRDX(p, jx+1)),
			   .5 * (CRDY(p, jy) + CRDY(p, jy+1)),
			   .5 * (CRDZ(p, jz) + CRDZ(p, jz+1)) };
	  // FIXME, the issue really is that (2nd order) particle pushers
	  // don't handle the invariant dim right
	  if (psc->domain.gdims[0] == 1) xx[0] = CRDX(p, jx);
	  if (psc->domain.gdims[1] == 1) xx[1] = CRDY(p, jy);
	  if (psc->domain.gdims[2] == 1) xx[2] = CRDZ(p, jz);
	  
	  int n_q_in_cell = 0;
	  for (int kind = 0; kind < nr_pop; kind++) {
	    struct psc_particle_npt npt = {};
	    if (kind < psc->nr_kinds) {
	      npt.kind = kind;
	      npt.q    = psc->kinds[kind].q;
	      npt.m    = psc->kinds[kind].m;
	      npt.n    = psc->kinds[kind].n;
	      npt.T[0] = psc->kinds[kind].T;
	      npt.T[1] = psc->kinds[kind].T;
	      npt.T[2] = psc->kinds[kind].T;
	    };
	    psc_ops(psc)->init_npt(psc, kind, xx, &npt);

	    int n_in_cell;
	    if (kind != psc->prm.neutralizing_population) {
	      n_in_cell = get_n_in_cell(psc, &npt);
	      n_q_in_cell += npt.q * n_in_cell;
	    } else {
	      // FIXME, should handle the case where not the last population is neutralizing
	      assert(psc->prm.neutralizing_population == nr_pop - 1);
	      n_in_cell = -n_q_in_cell / npt.q;
	    }
	    for (int cnt = 0; cnt < n_in_cell; cnt++) {
	      particle_t *prt = particle_iter_at(prts.begin, i++);
	      
	      psc_setup_particle(psc, prt, &npt, p, xx);
	      //p->lni = particle_label_offset + 1;
	      if (psc->prm.fractional_n_particles_per_cell) {
		prt->qni_wni = psc->kinds[prt->kind].q;
	      } else if (psc->prm.fortran_particle_weight_hack) {
		prt->qni_wni = psc->kinds[prt->kind].q * npt.n;
	      } else {
		prt->qni_wni = psc->kinds[prt->kind].q * npt.n / (n_in_cell * psc->coeff.cori);
	      }
	    }
	  }
	}
      }
    }
    psc_mparticles_patch_resize(mprts, p, i);
    if (!psc->prm.fractional_n_particles_per_cell) {
      assert(i == nr_particles_by_patch[p]);
    }
  }
  psc_mparticles_put_as(mprts, psc->particles, 0);
}

// ----------------------------------------------------------------------
// psc_setup_fields_default

static void
psc_setup_fields_default(struct psc *psc)
{
  double (*init_field)(struct psc *psc, double x[3], int m);
  init_field = psc_ops(psc)->init_field;
  if (!init_field)
    return;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_t *pf_base = psc_mfields_get_patch(psc->flds, p);
    fields_t *pf = psc_fields_get_as(pf_base, "c", 0, 0);

    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dx = psc->patch[p].dx[0], dy = psc->patch[p].dx[1], dz = psc->patch[p].dx[2];
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

      F3(pf, HX, jx,jy,jz) +=
	init_field(psc, (double []) { xx        , yy + .5*dy, zz + .5*dz }, HX);
      F3(pf, HY, jx,jy,jz) +=
	init_field(psc, (double []) { xx + .5*dx, yy        , zz + .5*dz }, HY);
      F3(pf, HZ, jx,jy,jz) +=
	init_field(psc, (double []) { xx + .5*dx, yy + .5*dy, zz         }, HZ);

      F3(pf, EX, jx,jy,jz) +=
	init_field(psc, (double []) { xx + .5*dx, yy        , zz         }, EX);
      F3(pf, EY, jx,jy,jz) +=
	init_field(psc, (double []) { xx        , yy + .5*dy, zz         }, EY);
      F3(pf, EZ, jx,jy,jz) +=
	init_field(psc, (double []) { xx        , yy        , zz + .5*dz }, EZ);

      F3(pf, JXI, jx,jy,jz) +=
	init_field(psc, (double []) { xx + .5*dx, yy        , zz         }, JXI);
      F3(pf, JYI, jx,jy,jz) +=
	init_field(psc, (double []) { xx        , yy + .5*dy, zz         }, JYI);
      F3(pf, JZI, jx,jy,jz) +=
	init_field(psc, (double []) { xx        , yy        , zz + .5*dz }, JZI);

    } foreach_3d_g_end;
    psc_fields_put_as(pf, pf_base, JXI, HX + 3);
  }
}


// ----------------------------------------------------------------------
// psc_setup_field_pml helper

static void
psc_setup_field_pml(struct psc *psc)
{
  mfields_base_t *flds = psc->flds;

  psc_mfields_copy_comp(flds, DX, flds, EX);
  psc_mfields_copy_comp(flds, DY, flds, EY);
  psc_mfields_copy_comp(flds, DZ, flds, EZ);
  psc_mfields_copy_comp(flds, BX, flds, HX);
  psc_mfields_copy_comp(flds, BY, flds, HY);
  psc_mfields_copy_comp(flds, BZ, flds, HZ);
  psc_mfields_set_comp(flds, EPS, 1.);
  psc_mfields_set_comp(flds, MU, 1.);
}

// ----------------------------------------------------------------------
// psc_setup_fields

void
psc_setup_fields(struct psc *psc)
{
  // create fields
  psc->flds = psc_mfields_create(mrc_domain_comm(psc->mrc_domain));
  psc_mfields_list_add(&psc_mfields_base_list, &psc->flds);
  psc_mfields_set_type(psc->flds, psc->prm.fields_base);
  psc_mfields_set_name(psc->flds, "mfields");
  psc_mfields_set_domain(psc->flds, psc->mrc_domain);
  psc_mfields_set_param_int(psc->flds, "nr_fields", NR_FIELDS);
  psc_mfields_set_param_int3(psc->flds, "ibn", psc->ibn);
  psc_mfields_setup(psc->flds);

  // type-specific other initial condition
  if (psc_ops(psc)->setup_fields) {
    psc_ops(psc)->setup_fields(psc, psc->flds);
  } else {
    psc_setup_fields_default(psc);
  }

  if (psc->domain.use_pml) {
    psc_setup_field_pml(psc);
  }
}

// ----------------------------------------------------------------------
// _psc_view

static void
_psc_view(struct psc *psc)
{
  mrc_domain_view(psc->mrc_domain);

  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "%20s|\n", "particle kinds");
  for (int k = 0; k < psc->nr_kinds; k++) {
    mpi_printf(comm, "%19s | q = %g m = %g n = %g T = %g\n", 
	       psc->kinds[k].name, psc->kinds[k].q, psc->kinds[k].m,
	       psc->kinds[k].n, psc->kinds[k].T);
  }
}

// ----------------------------------------------------------------------
// psc_set_kinds

void
psc_set_kinds(struct psc *psc, int nr_kinds, const struct psc_kind *kinds)
{
  if (!kinds && nr_kinds == psc->nr_kinds) {
    return;
  }

  if (psc->kinds) {
    for (int k = 0; k < psc->nr_kinds; k++) {
      free(psc->kinds[k].name);
    }
    free(psc->kinds);
  }
    
  psc->nr_kinds = nr_kinds;
  psc->kinds = calloc(nr_kinds, sizeof(*psc->kinds));
  if (kinds) {
    for (int k = 0; k < nr_kinds; k++) {
      psc->kinds[k] = kinds[k];
      psc->kinds[k].name = strdup(kinds[k].name);
    }
  } else {
    // set defaults, one electron species, the rest ions
    if (nr_kinds > KIND_ELECTRON) {
      psc->kinds[KIND_ELECTRON].name = strdup("e");
      psc->kinds[KIND_ELECTRON].q = -1.;
      psc->kinds[KIND_ELECTRON].m = 1.;
    }
    for (int k = 1; k < nr_kinds; k++) {
      char s[10];
      if (k == KIND_ION) {
	sprintf(s, "i");
      } else {
	sprintf(s, "i%d", k);
      }
      psc->kinds[k].name = strdup(s);
      psc->kinds[k].q = 1.;
      psc->kinds[k].m = 100.;
    }
  }
}

// ======================================================================
// psc class

struct mrc_class_psc mrc_class_psc = {
  .name             = "psc",
  .size             = sizeof(struct psc),
  .param_descr      = psc_descr,
  .create           = _psc_create,
  .set_from_options = _psc_set_from_options,
  .setup            = _psc_setup,
  .view             = _psc_view,
  .destroy          = _psc_destroy,
  .write            = _psc_write,
  .read             = _psc_read,
};

// ======================================================================

const char *fldname[NR_FIELDS] = {
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

// ======================================================================
// helpers

// ----------------------------------------------------------------------
// psc_default_dimensionless
//
// sets up parameter defaults for dimensionless units

void
psc_default_dimensionless(struct psc *psc)
{
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;
}

