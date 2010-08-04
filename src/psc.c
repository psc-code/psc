
#include "psc.h"
#include "util/params.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

struct psc psc;

// ----------------------------------------------------------------------
// assert_equal
//
// make sure that the two values are almost equal.

#define assert_equal(x, y, thres) __assert_equal(x, y, #x, #y, thres)

void
__assert_equal(double x, double y, const char *xs, const char *ys, double thres)
{
  double max = fmax(fabs(x), fabs(y)) + 1e-10;
  double eps = fabs((x - y) / max);
  if (eps > thres) {
    fprintf(stderr, "assert_equal: fail %s = %g, %s = %g rel err = %g thres = %g\n",
	    xs, x, ys, y, eps, thres);
    abort();
  }
}

// ----------------------------------------------------------------------

static struct psc_ops *psc_ops_list[] = {
  &psc_ops_generic_c,
  &psc_ops_fortran,
#ifdef USE_CUDA
  &psc_ops_cuda,
#endif
#ifdef USE_SSE2
  &psc_ops_sse2,
#endif
  NULL,
};

static struct psc_push_field_ops *psc_push_field_ops_list[] = {
  &psc_push_field_ops_fortran,
  &psc_push_field_ops_c,
  NULL,
};

static struct psc_randomize_ops *psc_randomize_ops_list[] = {
  &psc_randomize_ops_fortran,
  &psc_randomize_ops_none,
  NULL,
};

static struct psc_sort_ops *psc_sort_ops_list[] = {
  &psc_sort_ops_fortran,
  &psc_sort_ops_qsort,
  &psc_sort_ops_countsort,
  &psc_sort_ops_countsort2,
  &psc_sort_ops_none,
  NULL,
};

static struct psc_collision_ops *psc_collision_ops_list[] = {
  &psc_collision_ops_fortran,
  &psc_collision_ops_none,
  NULL,
};

static struct psc_output_ops *psc_output_ops_list[] = {
  &psc_output_ops_fortran,
  &psc_output_ops_c,
  NULL,
};

static struct psc_bnd_ops *psc_bnd_ops_list[] = {
  &psc_bnd_ops_fortran,
  &psc_bnd_ops_c,
  NULL,
};

static struct psc_ops *
psc_find_ops(const char *ops_name)
{
  for (int i = 0; psc_ops_list[i]; i++) {
    if (strcasecmp(psc_ops_list[i]->name, ops_name) == 0)
      return psc_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_ops '%s' not available.\n", ops_name);
  abort();
}

static struct psc_push_field_ops *
psc_find_push_field_ops(const char *ops_name)
{
  for (int i = 0; psc_push_field_ops_list[i]; i++) {
    if (strcasecmp(psc_push_field_ops_list[i]->name, ops_name) == 0)
      return psc_push_field_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_push_field_ops '%s' not available.\n", ops_name);
  abort();
}

static struct psc_randomize_ops *
psc_find_randomize_ops(const char *ops_name)
{
  for (int i = 0; psc_randomize_ops_list[i]; i++) {
    if (strcasecmp(psc_randomize_ops_list[i]->name, ops_name) == 0)
      return psc_randomize_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_randomize_ops '%s' not available.\n", ops_name);
  abort();
}

static struct psc_sort_ops *
psc_find_sort_ops(const char *ops_name)
{
  for (int i = 0; psc_sort_ops_list[i]; i++) {
    if (strcasecmp(psc_sort_ops_list[i]->name, ops_name) == 0)
      return psc_sort_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_sort_ops '%s' not available.\n", ops_name);
  abort();
}

static struct psc_collision_ops *
psc_find_collision_ops(const char *ops_name)
{
  for (int i = 0; psc_collision_ops_list[i]; i++) {
    if (strcasecmp(psc_collision_ops_list[i]->name, ops_name) == 0)
      return psc_collision_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_collision_ops '%s' not available.\n", ops_name);
  abort();
}

static struct psc_output_ops *
psc_find_output_ops(const char *ops_name)
{
  for (int i = 0; psc_output_ops_list[i]; i++) {
    if (strcasecmp(psc_output_ops_list[i]->name, ops_name) == 0)
      return psc_output_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_output_ops '%s' not available.\n", ops_name);
  abort();
}

static struct psc_bnd_ops *
psc_find_bnd_ops(const char *ops_name)
{
  for (int i = 0; psc_bnd_ops_list[i]; i++) {
    if (strcasecmp(psc_bnd_ops_list[i]->name, ops_name) == 0)
      return psc_bnd_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_bnd_ops '%s' not available.\n", ops_name);
  abort();
}

#define VAR(x) (void *)offsetof(struct psc_mod_config, x)

static struct param psc_mod_config_descr[] = {
  { "mod_particle"    , VAR(mod_particle)       , PARAM_STRING(NULL)  },
  { "mod_field"       , VAR(mod_field)          , PARAM_STRING(NULL)  },
  { "mod_randomize"   , VAR(mod_randomize)      , PARAM_STRING(NULL)  },
  { "mod_sort"        , VAR(mod_sort)           , PARAM_STRING(NULL)  },
  { "mod_collision"   , VAR(mod_collision)      , PARAM_STRING(NULL)  },
  { "mod_output"      , VAR(mod_output)         , PARAM_STRING(NULL)  },
  { "mod_bnd"         , VAR(mod_bnd)            , PARAM_STRING(NULL)  },
  {},
};

#undef VAR

void
psc_create(struct psc_mod_config *conf)
{
  // defaults
  if (!conf->mod_particle)
    conf->mod_particle = "fortran";
  if (!conf->mod_field)
    conf->mod_field = "fortran";
  if (!conf->mod_randomize)
    conf->mod_randomize = "none";
  if (!conf->mod_sort)
    conf->mod_sort = "none";
  if (!conf->mod_collision)
    conf->mod_collision = "none";
  if (!conf->mod_output)
    conf->mod_output = "fortran";
  if (!conf->mod_bnd)
    conf->mod_bnd = "fortran";

  params_parse_cmdline_nodefault(conf, psc_mod_config_descr, "PSC", MPI_COMM_WORLD);
  params_print(conf, psc_mod_config_descr, "PSC", MPI_COMM_WORLD);

  memset(&psc, 0, sizeof(psc));

  psc.ops = psc_find_ops(conf->mod_particle);
  if (psc.ops->create) {
    psc.ops->create();
  }
  psc.push_field_ops = psc_find_push_field_ops(conf->mod_field);
  if (psc.push_field_ops->create) {
    psc.push_field_ops->create();
  }
  psc.randomize_ops = psc_find_randomize_ops(conf->mod_randomize);
  if (psc.randomize_ops->create) {
    psc.randomize_ops->create();
  }
  psc.sort_ops = psc_find_sort_ops(conf->mod_sort);
  if (psc.sort_ops->create) {
    psc.sort_ops->create();
  }
  psc.collision_ops = psc_find_collision_ops(conf->mod_collision);
  if (psc.collision_ops->create) {
    psc.collision_ops->create();
  }
  psc.output_ops = psc_find_output_ops(conf->mod_output);
  if (psc.output_ops->create) {
    psc.output_ops->create();
  }
  psc.bnd_ops = psc_find_bnd_ops(conf->mod_bnd);
  if (psc.bnd_ops->create) {
    psc.bnd_ops->create();
  }
}

// ----------------------------------------------------------------------

// FIXME, obsolete
void
psc_alloc(int ilo[3], int ihi[3], int ibn[3], int n_part)
{
  for (int d = 0; d < 3; d++) {
    psc.ilo[d] = ilo[d];
    psc.ilg[d] = ilo[d] - ibn[d];
    psc.ihi[d] = ihi[d];
    psc.ihg[d] = ihi[d] + ibn[d];
    psc.img[d] = ihi[d] - ilo[d] + 2 * ibn[d];
    psc.ibn[d] = ibn[d];
    // for now, local == global (tests running on 1 proc)
    psc.domain.ilo[d] = ilo[d];
    psc.domain.ihi[d] = ihi[d];
  }
  psc.fld_size = psc.img[0] * psc.img[1] * psc.img[2];
  psc_fields_base_alloc(&psc.pf);

  psc_particles_base_alloc(&psc.pp, n_part);
  psc_set_n_particles(n_part);
}

void
psc_destroy()
{
  if (psc.ops->destroy) {
    psc.ops->destroy();
  }

  psc_fields_base_free(&psc.pf);
  psc_particles_base_free(&psc.pp);
}

void
psc_setup_parameters()
{
  psc.coeff.cori = 2.;
  psc.coeff.eta = 3.;
  psc.coeff.alpha = 5.;
  psc.dt = 1.;
  psc.dx[0] = 1.;
  psc.dx[1] = 1.;
  psc.dx[2] = 1.;
}

void
psc_setup_fields_zero()
{
  for (int m = 0; m < NR_FIELDS; m++) {
    psc_fields_base_zero(&psc.pf, m);
  }
}

void
psc_setup_fields_1()
{
  psc_setup_fields_zero();
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	F3_BASE(EX, jx,jy,jz) = .1 * sin(.5 * jx) + .2 * sin(.4 * jy) + .3 * sin(.3 * jz);
	F3_BASE(EY, jx,jy,jz) = .2 * sin(.4 * jx) + .3 * sin(.3 * jy) + .4 * sin(.2 * jz);
	F3_BASE(EZ, jx,jy,jz) = .3 * sin(.3 * jx) + .4 * sin(.2 * jy) + .5 * sin(.1 * jz);
	F3_BASE(BX, jx,jy,jz) = .1 * cos(.5 * jx) + .2 * cos(.4 * jy) + .3 * cos(.3 * jz);
	F3_BASE(BY, jx,jy,jz) = .2 * cos(.4 * jx) + .3 * cos(.3 * jy) + .4 * cos(.2 * jz);
	F3_BASE(BZ, jx,jy,jz) = .3 * cos(.3 * jx) + .4 * cos(.2 * jy) + .5 * cos(.1 * jz);
      }
    }
  }
}

static void
ascii_dump_field(int m, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  char *filename = malloc(strlen(fname) + 10);
  sprintf(filename, "%s-p%d.asc", fname, rank);
  mpi_printf(MPI_COMM_WORLD, "ascii_dump_field: '%s'\n", filename);

  FILE *file = fopen(filename, "w");
  for (int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++) {
    for (int iy = psc.ilg[1]; iy < psc.ihg[1]; iy++) {
      for (int ix = psc.ilg[0]; ix < psc.ihg[0]; ix++) {
	fprintf(file, "%d %d %d %g\n", ix, iy, iz, F3_BASE(m, ix,iy,iz));
      }
      fprintf(file, "\n");
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

static void
ascii_dump_particles(const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  char *filename = malloc(strlen(fname) + 10);
  sprintf(filename, "%s-p%d.asc", fname, rank);
  mpi_printf(MPI_COMM_WORLD, "ascii_dump_particles: '%s'\n", filename);

  FILE *file = fopen(filename, "w");
  fprintf(file, "i\txi\tyi\tzi\tpxi\tpyi\tpzi\tqni\tmni\twni\n");
  for (int i = 0; i < psc.pp.n_part; i++) {
    particle_base_t *p = psc_particles_base_get_one(&psc.pp, i);
    fprintf(file, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	    i, p->xi, p->yi, p->zi,
	    p->pxi, p->pyi, p->pzi, p->qni, p->mni, p->wni);
  }
  fclose(file);

  free(filename);
}

void
psc_dump_field(int m, const char *fname)
{
  if (psc.output_ops->dump_field) {
    psc.output_ops->dump_field(m, fname);
  } else {
    ascii_dump_field(m, fname);
  }
}

void
psc_dump_particles(const char *fname)
{
  if (psc.output_ops->dump_particles) {
    psc.output_ops->dump_particles(fname);
  } else {
    ascii_dump_particles(fname);
  }
}

void
psc_setup_particles_1()
{
  int n = 0;
  int n_per_cell = psc.pp.n_part / 
    ((psc.ihi[0]-psc.ilo[0])*(psc.ihi[1]-psc.ilo[1])*(psc.ihi[2]-psc.ilo[2]));
  for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
    for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
      for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {
	for (int cnt = 0; cnt < n_per_cell; cnt++) {
	  psc.pp.particles[n].xi = ix;
	  psc.pp.particles[n].yi = iy;
	  psc.pp.particles[n].zi = iz;
	  psc.pp.particles[n].pxi = .03;
	  psc.pp.particles[n].pyi = .02;
	  psc.pp.particles[n].pzi = .01;
	  psc.pp.particles[n].qni = -1.;
	  psc.pp.particles[n].mni = 1.;
	  psc.pp.particles[n].lni = 0.;
	  psc.pp.particles[n].wni = 1.;
	  n++;
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_push_part_xz

void
psc_push_part_xz()
{
  assert(psc.ops->push_part_xz);
  if (psc.ops->particles_from_fortran)
    psc.ops->particles_from_fortran();
  if (psc.ops->fields_from_fortran)
    psc.ops->fields_from_fortran();

  psc.ops->push_part_xz();

  if (psc.ops->particles_to_fortran)
    psc.ops->particles_to_fortran();
  if (psc.ops->fields_to_fortran)
    psc.ops->fields_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_part_yz

void
psc_push_part_yz()
{
  assert(psc.ops->push_part_yz);
  if (psc.ops->particles_from_fortran)
    psc.ops->particles_from_fortran();
  if (psc.ops->fields_from_fortran)
    psc.ops->fields_from_fortran();
  
  psc.ops->push_part_yz();

  if (psc.ops->particles_to_fortran)
    psc.ops->particles_to_fortran();
  if (psc.ops->fields_to_fortran)
    psc.ops->fields_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_part_z

void
psc_push_part_z()
{
  assert(psc.ops->push_part_z);
  if (psc.ops->particles_from_fortran)
    psc.ops->particles_from_fortran();
  if (psc.ops->fields_from_fortran)
    psc.ops->fields_from_fortran();

  psc.ops->push_part_z();

  if (psc.ops->fields_to_fortran)
    psc.ops->particles_to_fortran();
  if (psc.ops->particles_to_fortran)
    psc.ops->fields_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_part_yz_a

void
psc_push_part_yz_a()
{
  assert(psc.ops->push_part_yz_a);
  if (psc.ops->particles_from_fortran)
    psc.ops->particles_from_fortran();
  psc.ops->push_part_yz_a();
  if (psc.ops->particles_to_fortran)
    psc.ops->particles_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_part_yz_b

void
psc_push_part_yz_b()
{
  assert(psc.ops->push_part_yz_b);
  if (psc.ops->particles_from_fortran)
    psc.ops->particles_from_fortran();
  if (psc.ops->fields_from_fortran)
    psc.ops->fields_from_fortran();

  psc.ops->push_part_yz_b();

  if (psc.ops->fields_to_fortran)
    psc.ops->fields_to_fortran();
  if (psc.ops->particles_to_fortran)
    psc.ops->particles_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_particles

void
psc_push_particles()
{
  int im[3] = {
    psc.domain.ihi[0] - psc.domain.ilo[0],
    psc.domain.ihi[1] - psc.domain.ilo[1],
    psc.domain.ihi[2] - psc.domain.ilo[2],
  };
  if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
    assert(0);
  } else if (im[0] > 1 && im[2] > 1) { // xz
    psc_push_part_xz();
  } else if (im[0] > 1 && im[1] > 1) { // xy
    assert(0);
  } else if (im[1] > 1 && im[2] > 1) { // yz
    psc_push_part_yz();
  } else if (im[2] > 1) { // z
    psc_push_part_z();
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_push_field_a

void
psc_push_field_a()
{
  assert(psc.push_field_ops->push_field_a);
  psc.push_field_ops->push_field_a();
}

// ----------------------------------------------------------------------
// psc_push_field_b

void
psc_push_field_b()
{
  assert(psc.push_field_ops->push_field_b);
  psc.push_field_ops->push_field_b();
}

// ----------------------------------------------------------------------
// psc_add_ghosts

void
psc_add_ghosts(int mb, int me)
{
  assert(psc.bnd_ops->add_ghosts);
  psc.bnd_ops->add_ghosts(mb, me);
}

// ----------------------------------------------------------------------
// psc_fill_ghosts

void
psc_fill_ghosts(int mb, int me)
{
  assert(psc.bnd_ops->fill_ghosts);
  psc.bnd_ops->fill_ghosts(mb, me);
}

// ----------------------------------------------------------------------
// psc_exchange_particles

void
psc_exchange_particles(void)
{
  assert(psc.bnd_ops->exchange_particles);
  psc.bnd_ops->exchange_particles();
}

// ----------------------------------------------------------------------
// psc_randomize

void
psc_randomize()
{
  assert(psc.randomize_ops->randomize);
  psc.randomize_ops->randomize();
}

// ----------------------------------------------------------------------
// psc_sort

void
psc_sort()
{
  PIC_find_cell_indices();
  assert(psc.sort_ops->sort);
  psc.sort_ops->sort();
}

// ----------------------------------------------------------------------
// psc_collision

void
psc_collision()
{
  assert(psc.collision_ops->collision);
  psc.collision_ops->collision();
}

// ----------------------------------------------------------------------
// psc_out_field

void
psc_out_field()
{
  assert(psc.output_ops->out_field);
  psc.output_ops->out_field();
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

// psc_init

#define INIT_param_fortran_F77 F77_FUNC_(init_param_fortran, INIT_PARAM_FORTRAN)

void INIT_param_fortran_F77(void);

void
psc_init(const char *case_name)
{
  psc_init_param(case_name);

  SET_param_domain();
  SET_param_psc();
  SET_param_coeff();
  INIT_basic();
  INIT_param_fortran_F77();

  int n_part, particle_label_offset;
  psc_init_partition(&n_part, &particle_label_offset);
  SET_subdomain();

  psc_particles_base_alloc(&psc.pp, n_part);
  psc_init_particles(particle_label_offset);

  psc_fields_base_alloc(&psc.pf);
  psc_init_field();
}

// ----------------------------------------------------------------------
// psc_set_n_particles

void
psc_set_n_particles(int n_part)
{
  psc.pp.n_part = n_part;
  SET_niloc(n_part);
}

// ======================================================================
// testing related stuff


static particle_base_t *particle_ref;
static f_real *field_ref[NR_FIELDS];

// ----------------------------------------------------------------------
// psc_save_particles_ref
//
// save current particle data as reference solution

void
psc_save_particles_ref()
{
  if (!particle_ref) {
    particle_ref = calloc(psc.pp.n_part, sizeof(*particle_ref));
  }
  for (int i = 0; i < psc.pp.n_part; i++) {
    particle_ref[i] = *psc_particles_base_get_one(&psc.pp, i);
  }
}

// ----------------------------------------------------------------------
// psc_save_fields_ref
//
// save current field data as reference solution

void
psc_save_fields_ref()
{
  if (!field_ref[EX]) { //FIXME this is bad mojo
    for (int m = 0; m < NR_FIELDS; m++) {
      field_ref[m] = calloc(psc.fld_size, sizeof(f_real));
    }
  }
  for (int m = 0; m < NR_FIELDS; m++) {
    for (int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++) {
      for (int iy = psc.ilg[1]; iy < psc.ihg[1]; iy++) {
	for (int ix = psc.ilg[0]; ix < psc.ihg[0]; ix++) {
	  _FF3(field_ref[m], ix,iy,iz) = F3_BASE(m, ix,iy,iz);
	}
      }
    }
  }
} 



// ----------------------------------------------------------------------
// psc_check_particles_ref
//
// check current particle data agains previously saved reference solution

void
psc_check_particles_ref(double thres)
{
  assert(particle_ref);
  for (int i = 0; i < psc.pp.n_part; i++) {
    particle_base_t *part = psc_particles_base_get_one(&psc.pp, i);
    //    printf("i = %d\n", i);
    assert_equal(part->xi , particle_ref[i].xi, thres);
    assert_equal(part->yi , particle_ref[i].yi, thres);
    assert_equal(part->zi , particle_ref[i].zi, thres);
    assert_equal(part->pxi, particle_ref[i].pxi, thres);
    assert_equal(part->pyi, particle_ref[i].pyi, thres);
    assert_equal(part->pzi, particle_ref[i].pzi, thres);
  }
}


// ----------------------------------------------------------------------
// psc_check_fields_ref
//
// check field data against previously saved reference solution

void
psc_check_fields_ref(int *flds, double thres)
{
  assert(field_ref[EX]);
  for (int i = 0; flds[i] >= 0; i++) {
    int m = flds[i];
    for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
      for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
	for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {
	  //printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	  assert_equal(F3_BASE(m, ix,iy,iz), _FF3(field_ref[m], ix,iy,iz), thres);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_check_currents_ref
//
// check current current density data agains previously saved reference solution

void
psc_check_currents_ref(double thres)
{
  assert(field_ref[JXI]);
  for (int m = JXI; m <= JZI; m++){
    for (int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++) {
      for (int iy = psc.ilg[1]; iy < psc.ihg[1]; iy++) {
	for (int ix = psc.ilg[0]; ix < psc.ihg[0]; ix++) {
	  //	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	  assert_equal(F3_BASE(m, ix,iy,iz), _FF3(field_ref[m], ix,iy,iz), thres);
	}
      }
    }
  }
}

void
psc_check_currents_ref_noghost(double thres)
{
  assert(field_ref[JXI]);
  for (int m = JXI; m <= JZI; m++){
    for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
      for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
	for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {
	  //	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	  assert_equal(F3_BASE(m, ix,iy,iz), _FF3(field_ref[m], ix,iy,iz), thres);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_check_particles_sorted
//
// checks particles are sorted by cell index

void
psc_check_particles_sorted()
{
  int last = INT_MIN;

  for (int i = 0; i < psc.pp.n_part; i++) {
    assert(psc.pp.particles[i].cni >= last);
    last = psc.pp.particles[i].cni;
  }
}

// ----------------------------------------------------------------------
// psc_create_test_1
//
// set up test case 1

void
psc_create_test_1(const char *ops_name)
{
  int ilo[3] = { 0,  0,  0 };
  int ihi[3] = { 1, 16, 16 };
  int ibn[3] = { 2,  2,  2 }; // FIXME?

  int n_part = 1e3 * (ihi[2] - ilo[2]) * (ihi[1] - ilo[1]);

  struct psc_mod_config conf = {
    .mod_particle = ops_name,
  };
  psc_create(&conf);
  psc_alloc(ilo, ihi, ibn, n_part);
  psc_setup_parameters();
  psc_setup_fields_1();
  psc_setup_particles_1();
}

// ----------------------------------------------------------------------
// psc_create_test_xz

void
psc_create_test_xz(struct psc_mod_config *conf)
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  psc_create(conf);
  psc_init("test_xz");
}

// ----------------------------------------------------------------------
// psc_create_test_yz

void
psc_create_test_yz(struct psc_mod_config *conf)
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  psc_create(conf);
  psc_init("test_yz");
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
  [BX]  = "bx",
  [BY]  = "by",
  [BZ]  = "bz",
};
