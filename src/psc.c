
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

#include <mrc_common.h>
#include <mrc_params.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

struct psc psc;

#define VAR(x) (void *)offsetof(struct psc_mod_config, x)

static struct param psc_mod_config_descr[] = {
  { "mod_particle"    , VAR(mod_particle)       , PARAM_STRING(NULL)  },
  { "mod_field"       , VAR(mod_field)          , PARAM_STRING(NULL)  },
  { "mod_randomize"   , VAR(mod_randomize)      , PARAM_STRING(NULL)  },
  { "mod_sort"        , VAR(mod_sort)           , PARAM_STRING(NULL)  },
  { "mod_collision"   , VAR(mod_collision)      , PARAM_STRING(NULL)  },
  { "mod_output"      , VAR(mod_output)         , PARAM_STRING(NULL)  },
  { "mod_bnd"         , VAR(mod_bnd)            , PARAM_STRING(NULL)  },
  { "mod_moment"      , VAR(mod_moment)         , PARAM_STRING(NULL)  },
  {},
};

#undef VAR

void
psc_set_conf(struct psc *psc, struct psc_mod_config *conf)
{
  mrc_params_parse_nodefault(conf, psc_mod_config_descr, "PSC", MPI_COMM_WORLD);
  mrc_params_print(conf, psc_mod_config_descr, "PSC", MPI_COMM_WORLD);

  if (conf->mod_particle) {
    psc_push_particles_set_type(psc->push_particles, conf->mod_particle);
  }
  if (conf->mod_field) {
    psc_push_fields_set_type(psc->push_fields, conf->mod_field);
  }
  if (conf->mod_bnd) {
    psc_bnd_set_type(psc->bnd, conf->mod_bnd);
  }
  if (conf->mod_collision) {
    psc_collision_set_type(psc->collision, conf->mod_collision);
  }
  if (conf->mod_randomize) {
    psc_randomize_set_type(psc->randomize, conf->mod_randomize);
  }
  if (conf->mod_sort) {
    psc_sort_set_type(psc->sort, conf->mod_sort);
  }
  if (conf->mod_output) {
    psc_output_fields_set_type(psc->output_fields, conf->mod_output);
  }
  if (conf->mod_moment) {
    psc_moments_set_type(psc->moments, conf->mod_moment);
  }
}

// ----------------------------------------------------------------------
// psc_create

struct psc *
psc_create()
{
  memset(&psc, 0, sizeof(psc));

  MPI_Comm comm = MPI_COMM_WORLD;
  psc.push_particles = psc_push_particles_create(comm);
  psc.push_fields = psc_push_fields_create(comm);
  psc.bnd = psc_bnd_create(comm);
  psc.collision = psc_collision_create(comm);
  psc.randomize = psc_randomize_create(comm);
  psc.sort = psc_sort_create(comm);
  psc.output_fields = psc_output_fields_create(comm);
  psc.output_particles = psc_output_particles_create(comm);
  psc.moments = psc_moments_create(comm);

  psc.time_start = MPI_Wtime();

  psc_set_default_domain(&psc);
  psc_set_default_psc(&psc);

  return &psc;
}

// ----------------------------------------------------------------------
// psc_set_from_options

void
psc_set_from_options(struct psc *psc)
{
  psc_push_particles_set_from_options(psc->push_particles);
  psc_push_fields_set_from_options(psc->push_fields);
  psc_bnd_set_from_options(psc->bnd);
  psc_collision_set_from_options(psc->collision);
  psc_randomize_set_from_options(psc->randomize);
  psc_sort_set_from_options(psc->sort);
  psc_output_fields_set_from_options(psc->output_fields);
  psc_output_particles_set_from_options(psc->output_particles);
  psc_moments_set_from_options(psc->moments);

  psc_set_from_options_domain(psc);
  psc_set_from_options_psc(psc);
}

// ----------------------------------------------------------------------
// psc_setup

void
psc_setup(struct psc *psc)
{
  psc_setup_domain(psc); // needs to be done before setting up psc_bnd
  psc_setup_coeff(psc);

  psc_push_particles_setup(psc->push_particles);
  psc_push_fields_setup(psc->push_fields);
  psc_bnd_setup(psc->bnd);
  psc_collision_setup(psc->collision);
  psc_randomize_setup(psc->randomize);
  psc_sort_setup(psc->sort);
  psc_output_fields_setup(psc->output_fields);
  psc_output_particles_setup(psc->output_particles);
  psc_moments_setup(psc->moments);
}

// ----------------------------------------------------------------------
// psc_view

void
psc_view(struct psc *psc)
{
  psc_view_psc(psc);
  psc_view_domain(psc);

  psc_push_particles_view(psc->push_particles);
  psc_push_fields_view(psc->push_fields);
  psc_bnd_view(psc->bnd);
  psc_collision_view(psc->collision);
  psc_randomize_view(psc->randomize);
  psc_sort_view(psc->sort);
  psc_output_fields_view(psc->output_fields);
  psc_output_particles_view(psc->output_particles);
  psc_moments_view(psc->moments);
}

// ----------------------------------------------------------------------
// psc_destroy

void
psc_destroy(struct psc *psc)
{
  mfields_base_destroy(&psc->flds);
  mparticles_base_destroy(&psc->particles);

  psc_push_particles_destroy(psc->push_particles);
  psc_push_fields_destroy(psc->push_fields);
  psc_bnd_destroy(psc->bnd);
  psc_collision_destroy(psc->collision);
  psc_randomize_destroy(psc->randomize);
  psc_sort_destroy(psc->sort);
  psc_output_fields_destroy(psc->output_fields);
  psc_output_particles_destroy(psc->output_particles);
  psc_moments_destroy(psc->moments);
}

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
