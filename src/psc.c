
#include "psc.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_bnd.h"
#include "psc_collision.h"
#include "psc_randomize.h"
#include "psc_sort.h"
#include "psc_output_fields.h"
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
psc_create(struct psc_mod_config *conf)
{
  memset(&psc, 0, sizeof(psc));

  mrc_params_parse_nodefault(conf, psc_mod_config_descr, "PSC", MPI_COMM_WORLD);
  mrc_params_print(conf, psc_mod_config_descr, "PSC", MPI_COMM_WORLD);

  MPI_Comm comm = MPI_COMM_WORLD;

  psc.push_particles = psc_push_particles_create(comm);
  if (conf->mod_particle) {
    psc_push_particles_set_type(psc.push_particles, conf->mod_particle);
  }

  psc.push_fields = psc_push_fields_create(comm);
  if (conf->mod_field) {
    psc_push_fields_set_type(psc.push_fields, conf->mod_field);
  }

  psc.bnd = psc_bnd_create(comm);
  if (conf->mod_bnd) {
    psc_bnd_set_type(psc.bnd, conf->mod_bnd);
  }

  psc.collision = psc_collision_create(comm);
  if (conf->mod_collision) {
    psc_collision_set_type(psc.collision, conf->mod_collision);
  }

  psc.randomize = psc_randomize_create(comm);
  if (conf->mod_randomize) {
    psc_randomize_set_type(psc.randomize, conf->mod_randomize);
  }

  psc.sort = psc_sort_create(comm);
  if (conf->mod_sort) {
    psc_sort_set_type(psc.sort, conf->mod_sort);
  }

  psc.output_fields = psc_output_fields_create(comm);
  if (conf->mod_output) {
    psc_output_fields_set_type(psc.output_fields, conf->mod_output);
  }

  psc.moments = psc_moments_create(comm);
  if (conf->mod_moment) {
    psc_moments_set_type(psc.moments, conf->mod_moment);
  }

  psc.time_start = MPI_Wtime();
}

// ----------------------------------------------------------------------
// psc_set_from_options

void
psc_set_from_options(void)
{
  psc_push_particles_set_from_options(psc.push_particles);
  psc_push_fields_set_from_options(psc.push_fields);
  psc_bnd_set_from_options(psc.bnd);
  psc_collision_set_from_options(psc.collision);
  psc_randomize_set_from_options(psc.randomize);
  psc_sort_set_from_options(psc.sort);
  psc_output_fields_set_from_options(psc.output_fields);
  psc_moments_set_from_options(psc.moments);
}

// ----------------------------------------------------------------------
// psc_setup

void
psc_setup(void)
{
  psc_push_particles_setup(psc.push_particles);
  psc_push_fields_setup(psc.push_fields);
  psc_bnd_setup(psc.bnd);
  psc_collision_setup(psc.collision);
  psc_randomize_setup(psc.randomize);
  psc_sort_setup(psc.sort);
  psc_output_fields_setup(psc.output_fields);
  psc_moments_setup(psc.moments);
}

// ----------------------------------------------------------------------
// psc_view

void
psc_view(void)
{
  psc_push_particles_view(psc.push_particles);
  psc_push_fields_view(psc.push_fields);
  psc_bnd_view(psc.bnd);
  psc_collision_view(psc.collision);
  psc_randomize_view(psc.randomize);
  psc_sort_view(psc.sort);
  psc_output_fields_view(psc.output_fields);
  psc_moments_view(psc.moments);
}

// ----------------------------------------------------------------------
// psc_destroy

void
psc_destroy()
{
  mfields_base_destroy(&psc.flds);
  mparticles_base_destroy(&psc.particles);

  psc_push_particles_destroy(psc.push_particles);
  psc_push_fields_destroy(psc.push_fields);
  psc_bnd_destroy(psc.bnd);
  psc_collision_destroy(psc.collision);
  psc_randomize_destroy(psc.randomize);
  psc_sort_destroy(psc.sort);
  psc_output_fields_destroy(psc.output_fields);
  psc_moments_destroy(psc.moments);
}

// ----------------------------------------------------------------------
// psc_out_particles

void
psc_out_particles()
{
  assert(psc.nr_patches == 1);
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, &psc.particles);

  OUT_part(&particles.p[0]);
  
  particles_fortran_put(&particles, &psc.particles);
}

// ----------------------------------------------------------------------
// psc_init

#define INIT_param_fortran_F77 F77_FUNC_(init_param_fortran, INIT_PARAM_FORTRAN)

void INIT_param_fortran_F77(void);

void
psc_init(const char *case_name)
{
  psc_init_param(case_name);

  PSC_set_params();
  PSC_set_coeff();
  INIT_basic();
  INIT_param_fortran_F77();

  int particle_label_offset;
  psc_init_partition(&particle_label_offset);
  PSC_set_domain();
  PSC_set_patch(0);

  psc_init_particles(particle_label_offset);
  psc_set_from_options();
  psc_setup();
  psc_view();
  OUT_params_set();

  mfields_base_alloc(&psc.flds, NR_FIELDS);
  psc_init_field(&psc.flds);
}

// ----------------------------------------------------------------------
// psc_write_checkpoint

void
psc_read_checkpoint(void)
{
  mpi_printf(MPI_COMM_WORLD, "INFO: Reading checkpoint.\n");
  
  int n_part;
  SERV_read_1(&psc.timestep, &n_part);
  particles_base_realloc(&psc.particles.p[0], n_part);
  psc.particles.p[0].n_part = n_part;
  
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, &psc.particles);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, 0, 0, &psc.flds);

  SERV_read_2(&particles.p[0], &flds.f[0]);

  particles_fortran_put(&particles, &psc.particles);
  fields_fortran_put(&flds, JXI, HZ + 1, &psc.flds);
}

// ----------------------------------------------------------------------
// psc_write_checkpoint

void
psc_write_checkpoint(void)
{
  mpi_printf(MPI_COMM_WORLD, "INFO: Writing checkpoint.\n");
  
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, &psc.particles);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, JXI, HZ + 1, &psc.flds);

  SERV_write(&particles.p[0], &flds.f[0]);

  particles_fortran_put(&particles, &psc.particles);
  fields_fortran_put(&flds, 0, 0, &psc.flds);
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
