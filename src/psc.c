
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

void
psc_fields_destroy(mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_free(&flds->f[p]);
  }
  free(flds->f);
}

void
psc_particles_destroy(mparticles_base_t *particles)
{
  foreach_patch(p) {
    particles_base_free(&particles->p[p]);
  }
  free(particles->p);
}

// ----------------------------------------------------------------------
// psc_destroy

void
psc_destroy()
{
  psc_fields_destroy(&psc.flds);
  psc_particles_destroy(&psc.particles);

  psc_push_particles_destroy(psc.push_particles);
  psc_push_fields_destroy(psc.push_fields);
  psc_bnd_destroy(psc.bnd);
  psc_collision_destroy(psc.collision);
  psc_randomize_destroy(psc.randomize);
  psc_sort_destroy(psc.sort);
  psc_output_fields_destroy(psc.output_fields);
  psc_moments_destroy(psc.moments);
}

static void
ascii_dump_field(mfields_base_t *flds, int m, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  foreach_patch(p) {
    char *filename = malloc(strlen(fname) + 10);
    sprintf(filename, "%s-p%d-p%d.asc", fname, rank, p);
    mpi_printf(MPI_COMM_WORLD, "ascii_dump_field: '%s'\n", filename);

    fields_base_t *pf = &flds->f[p];
    FILE *file = fopen(filename, "w");
    free(filename);
    foreach_patch(patch) {
      for (int iz = -psc.ibn[2]; iz < psc.patch[patch].ldims[2] + psc.ibn[2]; iz++) {
	for (int iy = -psc.ibn[1]; iy < psc.patch[patch].ldims[1] + psc.ibn[1]; iy++) {
	  for (int ix = -psc.ibn[0]; ix < psc.patch[patch].ldims[0] +  psc.ibn[0]; ix++) {
	    fprintf(file, "%d %d %d %g\n", ix, iy, iz, F3_BASE(pf, m, ix,iy,iz));
	  }
	  fprintf(file, "\n");
	}
	fprintf(file, "\n");
      }
    }
    fclose(file);
  }
}

static void
ascii_dump_particles(mparticles_base_t *particles, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  foreach_patch(p) {
    particles_base_t *pp = &particles->p[p];
    char *filename = malloc(strlen(fname) + 10);
    sprintf(filename, "%s-p%d-p%d.asc", fname, rank, p);
    mpi_printf(MPI_COMM_WORLD, "ascii_dump_particles: '%s'\n", filename);
    
    FILE *file = fopen(filename, "w");
    fprintf(file, "i\txi\tyi\tzi\tpxi\tpyi\tpzi\tqni\tmni\twni\n");
    for (int i = 0; i < pp->n_part; i++) {
      particle_base_t *p = particles_base_get_one(pp, i);
      fprintf(file, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      i, p->xi, p->yi, p->zi,
	      p->pxi, p->pyi, p->pzi, p->qni, p->mni, p->wni);
    }
    fclose(file);
    free(filename);
  }
}

void
psc_dump_field(mfields_base_t *flds, int m, const char *fname)
{
  ascii_dump_field(flds, m, fname);
}

void
psc_dump_particles(mparticles_base_t *particles, const char *fname)
{
  ascii_dump_particles(particles, fname);
}

// ----------------------------------------------------------------------
// psc_out_particles

void
psc_out_particles()
{
  OUT_part();
}

// ----------------------------------------------------------------------
// psc_p_pulse_z1

real
psc_p_pulse_z1(real x, real y, real z, real t)
{
  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_p_z1) { // default to Fortran
    return PSC_p_pulse_z1(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_p_z1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z1

real
psc_s_pulse_z1(real x, real y, real z, real t)
{
  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_s_z1) { // default to Fortran
    return PSC_s_pulse_z1(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_s_z1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z2

real
psc_p_pulse_z2(real x, real y, real z, real t)
{
  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_p_z2) { // default to Fortran
    return PSC_p_pulse_z2(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_p_z2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z2

real
psc_s_pulse_z2(real x, real y, real z, real t)
{
  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_s_z2) { // default to Fortran
    return PSC_s_pulse_z2(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_s_z2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_init

#define INIT_param_fortran_F77 F77_FUNC_(init_param_fortran, INIT_PARAM_FORTRAN)

void INIT_param_fortran_F77(void);

void
psc_init(const char *case_name)
{
  psc_init_param(case_name);

  SET_param_psc();
  SET_param_coeff();
  INIT_basic();
  INIT_param_fortran_F77();

  int particle_label_offset;
  psc_init_partition(&particle_label_offset);
  SET_param_domain();
  SET_subdomain();

  psc_init_particles(particle_label_offset);
  psc_set_from_options();
  psc_setup();
  psc_view();

  mfields_base_alloc(&psc.flds, NR_FIELDS);
  psc_init_field(&psc.flds);
}

// ----------------------------------------------------------------------
// psc_set_n_particles

void
psc_set_n_particles(particles_base_t *pp, int n_part)
{
  pp->n_part = n_part;
  SET_niloc(n_part);
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
  psc_set_n_particles(0, n_part);
  
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

void
mfields_base_alloc(mfields_base_t *flds, int nr_fields)
{
  flds->f = calloc(psc.nr_patches, sizeof(*flds->f));
  foreach_patch(p) {
    int ilg[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };
    int ihg[3] = { psc.patch[p].ldims[0] + psc.ibn[0],
		   psc.patch[p].ldims[1] + psc.ibn[1],
		   psc.patch[p].ldims[2] + psc.ibn[2] };
    fields_base_alloc(&flds->f[p], ilg, ihg, nr_fields);
  }
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
