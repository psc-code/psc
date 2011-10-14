
#include "psc_push_particles_private.h"
#include "psc_cuda.h"
#include "psc_bnd.h"
#include "psc_sort.h"

#include <mrc_profile.h>
#include <math.h>


static void
psc_push_particles_cuda_push_yz_a(struct psc_push_particles *push,
				  mparticles_base_t *particles_base,
				  mfields_base_t *flds_base)
{
  mparticles_cuda_t *particles
    = psc_mparticles_get_cuda(particles_base, MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
  mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, EX, EX + 6);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_part_yz_a", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);
    yz_a_set_constants(pp, pf);
    __cuda_push_part_yz_a(pp, pf);
  }
  prof_stop(pr);

  psc_mfields_put_cuda(flds, flds_base, JXI, JXI + 3);
  psc_mparticles_put_cuda(particles, particles_base);
}

#define PUSH_PART_B 2

static void
psc_push_particles_cuda_push_yz_b(struct psc_push_particles *push,
				  mparticles_base_t *particles_base,
				  mfields_base_t *flds_base)
{
  mparticles_cuda_t *particles
    = psc_mparticles_get_cuda(particles_base, MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
  mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, EX, EX + 6);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_part_yz_b", 1., 0, 0);
  }

  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);
#if PUSH_PART_B == 2
    yz_set_constants(pp, pf);
#else
    yz_b_set_constants(pp, pf);
#endif

#if PUSH_PART_B == 1
    __cuda_push_part_yz_b(pp, pf);
#elif PUSH_PART_B == 2
    yz_cuda_push_part_p2(pp, pf);
#elif PUSH_PART_B == 3
    __cuda_push_part_yz_b3(pp, pf);
#endif
  }
  prof_stop(pr);

  psc_mfields_put_cuda(flds, flds_base, JXI, JXI + 3);
  psc_mparticles_put_cuda(particles, particles_base);
}

static void
cuda_push_part(mparticles_base_t *particles_base,
	       mfields_base_t *flds_base,
	       void (*set_constants)(particles_cuda_t *, fields_cuda_t *),
	       void (*push_part_p1)(particles_cuda_t *, fields_cuda_t *, real **),
	       void (*push_part_p2)(particles_cuda_t *, fields_cuda_t *),
	       void (*push_part_p3)(particles_cuda_t *, fields_cuda_t *, real *, int),
	       void (*push_part_p4)(particles_cuda_t *, fields_cuda_t *, real *),
	       void (*push_part_p5)(particles_cuda_t *, fields_cuda_t *, real *))
{
  mparticles_cuda_t *particles =
    psc_mparticles_get_cuda(particles_base, MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
  mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, EX, EX + 6);

  static int pr, pr1, pr2, pr3, pr4, pr5;
  if (!pr) {
    pr  = prof_register("cuda_part", 1., 0, 0);
    pr1 = prof_register("cuda_part_p1", 1., 0, 0);
    pr2 = prof_register("cuda_part_p2", 1., 0, 0);
    pr3 = prof_register("cuda_part_p3", 1., 0, 0);
    pr4 = prof_register("cuda_part_p4", 1., 0, 0);
    pr5 = prof_register("cuda_part_p5", 1., 0, 0);
  }
  prof_start(pr);

  // d_scratch needs to be per patch
  assert(ppsc->nr_patches == 1);
  real *d_scratch;
  
  psc_foreach_patch(ppsc, p) {
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);

    set_constants(pp, pf);

    prof_start(pr1);
    push_part_p1(pp, pf, &d_scratch);
    prof_stop(pr1);
    
    prof_start(pr2);
    push_part_p2(pp, pf);
    prof_stop(pr2);
  }

  // FIXME, doing this here doesn't jive well with integrate.c wanting to do it..
  psc_mparticles_put_cuda(particles, particles_base);
  psc_bnd_exchange_particles(ppsc->bnd, particles_base);
  psc_sort_run(ppsc->sort, particles_base);
  particles = psc_mparticles_get_cuda(particles_base,
				      MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);

  psc_foreach_patch(ppsc, p) {
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);

    set_constants(pp, pf);
    prof_start(pr3);
    push_part_p3(pp, pf, d_scratch, 1);
    prof_stop(pr3);

    prof_start(pr4);
    push_part_p4(pp, pf, d_scratch);
    prof_stop(pr4);
    
    prof_start(pr5);
    push_part_p5(pp, pf, d_scratch);
    prof_stop(pr5);
  }
  
  prof_stop(pr);

  psc_mfields_put_cuda(flds, flds_base, JXI, JXI + 3);
  psc_mparticles_put_cuda(particles, particles_base);
}

static void
cuda_push_partq(struct psc_push_particles *push,
		mparticles_base_t *particles_base,
		mfields_base_t *flds_base,
		void (*set_constants)(particles_cuda_t *, fields_cuda_t *),
		void (*push_part_p2)(particles_cuda_t *, fields_cuda_t *),
		void (*push_part_p3)(particles_cuda_t *, fields_cuda_t *, real *, int))
{
  const int block_stride = 4;
  unsigned int mp_flags = psc_push_particles_get_mp_flags(push);
  
  // FIXME distinguish alloc/calc_block_offsets?
  mparticles_cuda_t *particles = 
    psc_mparticles_get_cuda(particles_base, mp_flags);
  mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, EX, EX + 6);

  static int pr, pr2, pr3;
  if (!pr) {
    pr  = prof_register("cuda_part", 1., 0, 0);
    pr2 = prof_register("cuda_part_p2", 1., 0, 0);
    pr3 = prof_register("cuda_part_p3", 1., 0, 0);
  }
  prof_start(pr);

  // d_scratch needs to be per patch
  assert(ppsc->nr_patches == 1);
  
  psc_foreach_patch(ppsc, p) {
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);

    set_constants(pp, pf);

    prof_start(pr2);
    push_part_p2(pp, pf);
    prof_stop(pr2);
  }

  // FIXME, doing this here doesn't jive well with integrate.c wanting to do it..
  psc_mparticles_put_cuda(particles, particles_base);
  psc_bnd_exchange_particles(ppsc->bnd, particles_base);

  // block/cell offsets will be calculated by sort_patch, anyway
  particles = psc_mparticles_get_cuda(particles_base,
				      mp_flags &~(MP_NEED_CELL_OFFSETS | MP_NEED_BLOCK_OFFSETS));
  psc_foreach_patch(ppsc, p) {
    if (mp_flags & MP_NEED_CELL_OFFSETS) {
      cuda_sort_patch_by_cell(p, psc_mparticles_get_patch_cuda(particles, p));
    } else {
      cuda_sort_patch(p, psc_mparticles_get_patch_cuda(particles, p));
    }
  }

  psc_foreach_patch(ppsc, p) {
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);

    set_constants(pp, pf);
    prof_start(pr3);
    push_part_p3(pp, pf, NULL, block_stride);
    prof_stop(pr3);
  }
  
  prof_stop(pr);

  psc_mfields_put_cuda(flds, flds_base, JXI, JXI + 3);
  psc_mparticles_put_cuda(particles, particles_base);
}

// ----------------------------------------------------------------------
// push_z

static void __unused
psc_push_particles_cuda_push_z(struct psc_push_particles *push,
			       mparticles_base_t *particles_base,
			       mfields_base_t *flds_base)
{
  cuda_push_part(particles_base, flds_base,
		 z_set_constants,
		 z_cuda_push_part_p1,
		 z_cuda_push_part_p2,
		 z_cuda_push_part_p3,
		 z_cuda_push_part_p4,
		 z_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_z2(struct psc_push_particles *push,
				mparticles_base_t *particles_base,
				mfields_base_t *flds_base)
{
  cuda_push_part(particles_base, flds_base,
		 z2_set_constants,
		 z2_cuda_push_part_p1,
		 z2_cuda_push_part_p2,
		 z2_cuda_push_part_p3,
		 z2_cuda_push_part_p4,
		 z2_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_z3(struct psc_push_particles *push,
				mparticles_base_t *particles_base,
				mfields_base_t *flds_base)
{
  cuda_push_part(particles_base, flds_base,
		 z3_set_constants,
		 z3_cuda_push_part_p1,
		 z3_cuda_push_part_p2,
		 z3_cuda_push_part_p3,
		 z3_cuda_push_part_p4,
		 z3_cuda_push_part_p5);
}

// ----------------------------------------------------------------------
// push_yz

static void __unused
psc_push_particles_cuda_push_yz(struct psc_push_particles *push,
				mparticles_base_t *particles_base,
				mfields_base_t *flds_base)
{
  cuda_push_part(particles_base, flds_base,
		 yz_set_constants,
		 yz_cuda_push_part_p1,
		 yz_cuda_push_part_p2,
		 yz_cuda_push_part_p3,
		 yz_cuda_push_part_p4,
		 yz_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz2(struct psc_push_particles *push,
				 mparticles_base_t *particles_base,
				 mfields_base_t *flds_base)
{
  cuda_push_part(particles_base, flds_base,
		 yz2_set_constants,
		 yz2_cuda_push_part_p1,
		 yz2_cuda_push_part_p2,
		 yz2_cuda_push_part_p3,
		 yz2_cuda_push_part_p4,
		 yz2_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz3(struct psc_push_particles *push,
				 mparticles_base_t *particles_base,
				 mfields_base_t *flds_base)
{
  cuda_push_part(particles_base, flds_base,
		 yz3_set_constants,
		 yz3_cuda_push_part_p1,
		 yz3_cuda_push_part_p2,
		 yz3_cuda_push_part_p3,
		 yz3_cuda_push_part_p4,
		 yz3_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz4(struct psc_push_particles *push,
				 mparticles_base_t *particles_base,
				 mfields_base_t *flds_base)
{
  cuda_push_part(particles_base, flds_base,
		 yz4_set_constants,
		 yz4_cuda_push_part_p1,
		 yz4_cuda_push_part_p2,
		 yz4_cuda_push_part_p3,
		 yz4_cuda_push_part_p4,
		 yz4_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz5(struct psc_push_particles *push,
				 mparticles_base_t *particles_base,
				 mfields_base_t *flds_base)
{
  cuda_push_partq(push, particles_base, flds_base,
		  yz5_set_constants,
		  yz5_cuda_push_part_p2,
		  yz5_cuda_push_part_p3);
}

static void __unused
psc_push_particles_cuda_push_yz6(struct psc_push_particles *push,
				 mparticles_base_t *particles_base,
				 mfields_base_t *flds_base)
{
  cuda_push_partq(push, particles_base, flds_base,
		  yz6_set_constants,
		  yz6_cuda_push_part_p2,
		  yz6_cuda_push_part_p3);
}

// ======================================================================
// psc_push_particles: subclass "cuda"

struct psc_push_particles_ops psc_push_particles_cuda_ops = {
  .name                  = "cuda",
  .push_z                = psc_push_particles_cuda_push_z3,
  .push_yz               = psc_push_particles_cuda_push_yz6,
  .push_yz_a             = psc_push_particles_cuda_push_yz_a,
  .push_yz_b             = psc_push_particles_cuda_push_yz_b,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS,
};

// ======================================================================

DECLARE_CUDA(yz4x4_1st);

static void
psc_push_particles_cuda_1st_push_yz(struct psc_push_particles *push,
				    mparticles_base_t *particles_base,
				    mfields_base_t *flds_base)
{
  cuda_push_partq(push, particles_base, flds_base,
		  yz4x4_1st_set_constants,
		  yz4x4_1st_cuda_push_part_p2,
		  yz4x4_1st_cuda_push_part_p3);
}

// ======================================================================
// psc_push_particles: subclass "cuda_1st"

struct psc_push_particles_ops psc_push_particles_cuda_1st_ops = {
  .name                  = "cuda_1st",
  .push_yz               = psc_push_particles_cuda_1st_push_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS |
                           MP_BLOCKSIZE_4X4X4,
};

// ======================================================================

DECLARE_CUDA(yz4x4_1vb);

static void
psc_push_particles_cuda_1vb_push_yz4x4(struct psc_push_particles *push,
				       mparticles_base_t *particles_base,
				       mfields_base_t *flds_base)
{
  cuda_push_partq(push, particles_base, flds_base,
		  yz4x4_1vb_set_constants,
		  yz4x4_1vb_cuda_push_part_p2,
		  yz4x4_1vb_cuda_push_part_p3);
}

// ======================================================================
// psc_push_particles: subclass "cuda_1vb"

struct psc_push_particles_ops psc_push_particles_cuda_1vb_ops = {
  .name                  = "cuda_1vb",
  .push_yz               = psc_push_particles_cuda_1vb_push_yz4x4,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4,
};
