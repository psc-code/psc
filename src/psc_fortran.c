
#include "psc.h"
#include <mrc_profile.h>

static void
fortran_push_part_xy(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_xy", 1., 0, 0);
  }
  prof_start(pr);
  PIC_push_part_xy(&particles.p[0], &flds.f[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

static void
fortran_push_part_xz(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_xz", 1., 0, 0);
  }
  prof_start(pr);
  PIC_push_part_xz(&particles.p[0], &flds.f[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

static void
fortran_push_part_yz(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_yz", 1., 0, 0);
  }
  prof_start(pr);
  PIC_push_part_yz(&particles.p[0], &flds.f[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

static void
fortran_push_part_xyz(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_xyz", 1., 0, 0);
  }
  prof_start(pr);
  PIC_push_part_xyz(&particles, &flds.f[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

static void
fortran_push_part_z(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_z", 1., 0, 0);
  }
  prof_start(pr);
  PIC_push_part_z(&particles.p[0], &flds.f[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

static void
fortran_push_part_yz_a(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_yz_a", 1., 0, 0);
  }
  prof_start(pr);
  PIC_push_part_yz_a(&particles.p[0], &flds.f[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

static void
fortran_push_part_yz_b(mfields_base_t *flds_base, mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, EX, EX + 6, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_yz_b", 1., 0, 0);
  }
  prof_start(pr);
  PIC_push_part_yz_b(&particles.p[0], &flds.f[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
  fields_fortran_put(&flds, JXI, JXI + 3, flds_base);
}

struct psc_ops psc_ops_fortran = {
  .name = "fortran",
  .push_part_xy           = fortran_push_part_xy,
  .push_part_xz           = fortran_push_part_xz,
  .push_part_yz           = fortran_push_part_yz,
  .push_part_xyz          = fortran_push_part_xyz,
  .push_part_z            = fortran_push_part_z,
  .push_part_yz_a         = fortran_push_part_yz_a,
  .push_part_yz_b         = fortran_push_part_yz_b,
};

// ======================================================================
// fortran randomize

static void
fortran_randomize(mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_randomize", 1., 0, 0);
  }
  prof_start(pr);
  PIC_randomize(&particles.p[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
}

struct psc_randomize_ops psc_randomize_ops_fortran = {
  .name      = "fortran",
  .randomize = fortran_randomize,
};

// ======================================================================
// fortran sort

static void
fortran_sort(mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_sort", 1., 0, 0);
  }
  prof_start(pr);
  PIC_find_cell_indices(&particles.p[0]);
  PIC_sort(&particles.p[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
}

struct psc_sort_ops psc_sort_ops_fortran = {
  .name = "fortran",
  .sort = fortran_sort,
};

// ======================================================================
// fortran collision

static void
fortran_collision(mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_collision", 1., 0, 0);
  }
  prof_start(pr);
  PIC_bin_coll(&particles.p[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
}

struct psc_collision_ops psc_collision_ops_fortran = {
  .name      = "fortran",
  .collision = fortran_collision,
};

// ======================================================================
// fortran push field

static void
fortran_push_field_a(mfields_base_t *flds_base)
{
  if (psc.domain.use_pml) {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, MU + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_pml_a", 1., 0, 0);
    }
    prof_start(pr);
    PIC_pml_msa(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, BZ + 1, flds_base);
  } else {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_a", 1., 0, 0);
    }
    prof_start(pr);
    PIC_msa(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, HZ + 1, flds_base);
  }
}

static void
fortran_push_field_b(mfields_base_t *flds_base)
{
  if (psc.domain.use_pml) {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, MU + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_pml_b", 1., 0, 0);
    }
    prof_start(pr);
    PIC_pml_msb(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, BZ + 1, flds_base);
  } else {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_b", 1., 0, 0);
    }
    prof_start(pr);
    PIC_msb(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, HZ + 1, flds_base);
  }
}

struct psc_push_field_ops psc_push_field_ops_fortran = {
  .name         = "fortran",
  .push_field_a = fortran_push_field_a,
  .push_field_b = fortran_push_field_b,
};

// ======================================================================
// fortran output

static void
fortran_out_field()
{
  assert(FIELDS_BASE == FIELDS_FORTRAN);
  static int pr;
  if (!pr) {
    pr = prof_register("fort_out_field", 1., 0, 0);
  }
  prof_start(pr);
  OUT_field();
  prof_stop(pr);
}

struct psc_output_ops psc_output_ops_fortran = {
  .name      = "fortran",
  .out_field = fortran_out_field,
};

// ======================================================================
// fortran bnd

static void
fortran_add_ghosts(mfields_base_t *flds_base, int mb, int me)
{
  assert(psc.nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  mfields_fortran_t flds;
  fields_fortran_get(&flds, mb, me, flds_base);

  for (int m = mb; m < me; m++) {
    if (psc.domain.gdims[0] > 1) {
      PIC_fax(&flds.f[0], m);
    }
    if (psc.domain.gdims[1] > 1) {
      PIC_fay(&flds.f[0], m);
    }
    if (psc.domain.gdims[2] > 1) {
      PIC_faz(&flds.f[0], m);
    }
  }

  fields_fortran_put(&flds, mb, me, flds_base);

  prof_stop(pr);
}

static void
fortran_fill_ghosts(mfields_base_t *flds_base, int mb, int me)
{
  assert(psc.nr_patches == 1);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_fill_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  mfields_fortran_t flds;
  fields_fortran_get(&flds, mb, me, flds_base);

  for (int m = mb; m < me; m++) {
    if (psc.domain.gdims[0] > 1) {
      PIC_fex(&flds.f[0], m);
    }
    if (psc.domain.gdims[1] > 1) {
      PIC_fey(&flds.f[0], m);
    }
    if (psc.domain.gdims[2] > 1) {
      PIC_fez(&flds.f[0], m);
    }
  }

  fields_fortran_put(&flds, mb, me, flds_base);

  prof_stop(pr);
}

static void
fortran_exchange_particles(mparticles_base_t *particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_xchg_part", 1., 0, 0);
  }
  prof_start(pr);

  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);
  
  SET_param_coeff();
  SET_niloc(particles.p[0].n_part);

  if (psc.domain.gdims[0] > 1) {
    PIC_pex();
  }
  if (psc.domain.gdims[1] > 1) {
    PIC_pey();
  }
  if (psc.domain.gdims[2] > 1) {
    PIC_pez();
  }

  GET_niloc(&particles.p[0].n_part);
  // don't really reallocate, just get the new array pointer
  // if PIC_pe[xyz]() reallocated during the previous calls
  particles.p[0].particles = REALLOC_particles(particles.p[0].n_part);
  particles_fortran_put(&particles, particles_base);

  prof_stop(pr);
}

struct psc_bnd_ops psc_bnd_ops_fortran = {
  .name               = "fortran",
  .add_ghosts         = fortran_add_ghosts,
  .fill_ghosts        = fortran_fill_ghosts,
  .exchange_particles = fortran_exchange_particles,
};

// ======================================================================
// fortran moment

static void
fortran_calc_densities(mfields_base_t *flds_base, mparticles_base_t *particles_base,
		       mfields_base_t *f)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_densities", 1., 0, 0);
  }
  prof_start(pr);

  mparticles_fortran_t particles;
  particles_fortran_get(&particles, &particles_base);
  mfields_fortran_t flds_fortran;
  fields_fortran_get_from(&flds_fortran, 0, 0, f, 0);

  CALC_densities(&particles.p[0], &flds_fortran.f[0]);

  particles_fortran_put(&particles, &psc.particles);
  fields_fortran_put_to(&flds_fortran, NE, NE + 3, f, 0);

  prof_stop(pr);
}

struct psc_moment_ops psc_moment_ops_fortran = {
  .name               = "fortran",
  .calc_densities     = fortran_calc_densities,
};

