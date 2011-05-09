
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

  psc_set_default_domain(psc);
  psc_set_default_psc(psc);
}

// ----------------------------------------------------------------------
// psc_set_from_options

static void
_psc_set_from_options(struct psc *psc)
{
  psc_set_from_options_domain(psc);
  psc_set_from_options_psc(psc);
}

// ----------------------------------------------------------------------
// psc_setup

static void
_psc_setup(struct psc *psc)
{
  psc_setup_coeff(psc);
  psc_setup_domain(psc); // needs to be done before setting up psc_bnd
}

// ----------------------------------------------------------------------
// psc_view

static void
_psc_view(struct psc *psc)
{
  psc_view_psc(psc);
  psc_view_domain(psc);
}

// ----------------------------------------------------------------------
// psc_destroy

static void
_psc_destroy(struct psc *psc)
{
  mfields_base_destroy(psc->flds);
  mparticles_base_destroy(&psc->particles);
  mphotons_destroy(&psc->mphotons);

  psc_destroy_domain(psc);
}

// ======================================================================
// psc class

struct mrc_class_psc mrc_class_psc = {
  .name             = "psc",
  .size             = sizeof(struct psc),
  .create           = _psc_create,
  .set_from_options = _psc_set_from_options,
  .setup            = _psc_setup,
  .destroy          = _psc_destroy,
  .view             = _psc_view,
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
