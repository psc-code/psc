
#include "psc_push_particles_private.h"
#include "psc_cuda.h"
#include "psc_bnd.h"
#include "psc_sort.h"
#include "psc_fields.h"

#include <mrc_profile.h>
#include <math.h>

EXTERN_C void yz4x4_1vb_cuda_push_mprts_a(struct psc_mparticles *mprts, struct psc_mfields *mflds);
EXTERN_C void yz4x4_1vb_cuda_push_mprts_b(struct psc_mparticles *mprts, struct psc_mfields *mflds);

static void
cuda_push_part(struct psc_particles *prts_base,
	       struct psc_fields *flds_base,
	       void (*set_constants)(struct psc_particles *, struct psc_fields *),
	       void (*push_part_p1)(struct psc_particles *, struct psc_fields *, real **),
	       void (*push_part_p2)(struct psc_particles *, struct psc_fields *),
	       void (*push_part_p3)(struct psc_particles *, struct psc_fields *, real *, int),
	       void (*push_part_p4)(struct psc_particles *, struct psc_fields *, real *),
	       void (*push_part_p5)(struct psc_particles *, struct psc_fields *, real *))
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, "cuda",
						    MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, EX + 6);

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
  real *d_scratch;
  
  set_constants(prts, flds);
  
  prof_start(pr1);
  push_part_p1(prts, flds, &d_scratch);
  prof_stop(pr1);
  
  prof_start(pr2);
  push_part_p2(prts, flds);
  prof_stop(pr2);

  assert(0);
#if 0
  // FIXME, doing this here doesn't jive well with integrate.c wanting to do it..
  psc_particles_put_as(particles, particles_base, 0);
  psc_bnd_exchange_particles(ppsc->bnd, particles_base);
  psc_sort_run(ppsc->sort, particles_base);
  prts = psc_particles_get_as(prts_base, "cuda",
			      MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
#endif

  set_constants(prts, flds);
  prof_start(pr3);
  push_part_p3(prts, flds, d_scratch, 1);
  prof_stop(pr3);
  
  prof_start(pr4);
  push_part_p4(prts, flds, d_scratch);
  prof_stop(pr4);
  
  prof_start(pr5);
  push_part_p5(prts, flds, d_scratch);
  prof_stop(pr5);
  
  prof_stop(pr);

  psc_fields_put_as(flds, flds_base, JXI, JXI + 3);
  psc_particles_put_as(prts, prts_base, 0);
}

static void
cuda_push_partq_a(struct psc_push_particles *push,
		  struct psc_particles *prts_base,
		  struct psc_fields *flds_base,
		  void (*push_part_p2)(struct psc_particles *, struct psc_fields *))
{
  unsigned int mp_flags = psc_push_particles_get_mp_flags(push);
  
  static int pr;
  if (!pr) {
    pr  = prof_register("cuda_part_a", 1., 0, 0);
  }

  struct psc_particles *prts = psc_particles_get_as(prts_base, "cuda", mp_flags);
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, EX + 6);
  
  prof_start(pr);
  push_part_p2(prts, flds);
  prof_stop(pr);
  
  psc_particles_put_as(prts, prts_base, 0);
  psc_fields_put_as(flds, flds_base, 0, 0);
}

static void
cuda_push_partq_b(struct psc_push_particles *push,
		  struct psc_particles *prts_base,
		  struct psc_fields *flds_base,
		  void (*push_part_p3)(struct psc_particles *, struct psc_fields *, real *, int))
{
  static int pr;
  if (!pr) {
    pr = prof_register("cuda_part_b", 1., 0, 0);
  }

  const int block_stride = 4;
  unsigned int mp_flags = psc_push_particles_get_mp_flags(push);

  struct psc_particles *prts = psc_particles_get_as(prts_base, "cuda", mp_flags);
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", 0, 0);

  prof_start(pr);
  push_part_p3(prts, flds, NULL, block_stride);
  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  psc_fields_put_as(flds, flds_base, JXI, JXI + 3);
}

// ----------------------------------------------------------------------
// push_z

static void __unused
psc_push_particles_cuda_push_z(struct psc_push_particles *push,
			       struct psc_particles *prts_base,
			       struct psc_fields *flds_base)
{
  cuda_push_part(prts_base, flds_base,
		 z_set_constants,
		 z_cuda_push_part_p1,
		 z_cuda_push_part_p2,
		 z_cuda_push_part_p3,
		 z_cuda_push_part_p4,
		 z_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_z2(struct psc_push_particles *push,
			       struct psc_particles *prts_base,
			       struct psc_fields *flds_base)
{
  cuda_push_part(prts_base, flds_base,
		 z2_set_constants,
		 z2_cuda_push_part_p1,
		 z2_cuda_push_part_p2,
		 z2_cuda_push_part_p3,
		 z2_cuda_push_part_p4,
		 z2_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_z3(struct psc_push_particles *push,
			       struct psc_particles *prts_base,
			       struct psc_fields *flds_base)
{
  cuda_push_part(prts_base, flds_base,
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
				struct psc_particles *prts_base,
				struct psc_fields *flds_base)
{
  cuda_push_part(prts_base, flds_base,
		 yz_set_constants,
		 yz_cuda_push_part_p1,
		 yz_cuda_push_part_p2,
		 yz_cuda_push_part_p3,
		 yz_cuda_push_part_p4,
		 yz_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz2(struct psc_push_particles *push,
				 struct psc_particles *prts_base,
				 struct psc_fields *flds_base)
{
  cuda_push_part(prts_base, flds_base,
		 yz2_set_constants,
		 yz2_cuda_push_part_p1,
		 yz2_cuda_push_part_p2,
		 yz2_cuda_push_part_p3,
		 yz2_cuda_push_part_p4,
		 yz2_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz3(struct psc_push_particles *push,
				 struct psc_particles *prts_base,
				 struct psc_fields *flds_base)
{
  cuda_push_part(prts_base, flds_base,
		 yz3_set_constants,
		 yz3_cuda_push_part_p1,
		 yz3_cuda_push_part_p2,
		 yz3_cuda_push_part_p3,
		 yz3_cuda_push_part_p4,
		 yz3_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz4(struct psc_push_particles *push,
				 struct psc_particles *prts_base,
				 struct psc_fields *flds_base)
{
  cuda_push_part(prts_base, flds_base,
		 yz4_set_constants,
		 yz4_cuda_push_part_p1,
		 yz4_cuda_push_part_p2,
		 yz4_cuda_push_part_p3,
		 yz4_cuda_push_part_p4,
		 yz4_cuda_push_part_p5);
}

static void __unused
psc_push_particles_cuda_push_yz5(struct psc_push_particles *push,
				 struct psc_particles *prts_base,
				 struct psc_fields *flds_base)
{
  assert(0);
  /* cuda_push_partq(push, prts_base, flds_base, */
  /* 		  yz5_set_constants, */
  /* 		  yz5_cuda_push_part_p2, */
  /* 		  yz5_cuda_push_part_p3); */
}

static void __unused
psc_push_particles_cuda_push_yz6(struct psc_push_particles *push,
				 mparticles_base_t *particles_base,
				 mfields_base_t *flds_base)
{
  assert(0);
  /* cuda_push_partq(push, particles_base, flds_base, */
  /* 		  yz6_set_constants, */
  /* 		  yz6_cuda_push_part_p2, */
  /* 		  yz6_cuda_push_part_p3); */
}

// ======================================================================
// psc_push_particles: subclass "cuda"

struct psc_push_particles_ops psc_push_particles_cuda_ops = {
  .name                  = "cuda",
  //  .push_z                = psc_push_particles_cuda_push_z3,
  //  .push_yz               = psc_push_particles_cuda_push_yz6,
  .mp_flags              = MP_BLOCKSIZE_4X4X4 |
                           MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS,
};

// ======================================================================

DECLARE_CUDA(yz4x4_1st);

static void
psc_push_particles_cuda_1st_push_a_yz(struct psc_push_particles *push,
				      struct psc_particles *prts_base,
				      struct psc_fields *flds_base)
{
  assert(0);
  /* cuda_push_partq_a(push, prts_base, flds_base, */
  /* 		    yz4x4_1st_set_constants, */
  /* 		    yz4x4_1st_cuda_push_part_p2, */
  /* 		    yz4x4_1st_cuda_push_part_p3); */
}

// ======================================================================
// psc_push_particles: subclass "cuda_1st"

struct psc_push_particles_ops psc_push_particles_cuda_1st_ops = {
  .name                  = "cuda_1st",
  .push_a_yz             = psc_push_particles_cuda_1st_push_a_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS |
                           MP_BLOCKSIZE_4X4X4,
};

// ======================================================================

DECLARE_CUDA(yz4x4_1vb);

static void
psc_push_particles_1vb_4x4_cuda_push_mprts_yz(struct psc_push_particles *push,
					      struct psc_mparticles *mprts,
					      struct psc_mfields *mflds_base)
{
  static int pr_A, pr_B, pr_D;
  if (!pr_A) {
    pr_A  = prof_register("cuda_part_a", 1., 0, 0);
    pr_B  = prof_register("cuda_part_b", 1., 0, 0);
    pr_D  = prof_register("xchg_REORDER", 1., 0, 0);
  }

  // it's difficult to convert mprts because of the ordering constraints (?)
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", EX, EX + 6);

  prof_start(pr_A);
  yz4x4_1vb_cuda_push_mprts_a(mprts, mflds);
  prof_stop(pr_A);
  
  prof_start(pr_B);
  yz4x4_1vb_cuda_push_mprts_b(mprts, mflds);
  prof_stop(pr_B);

  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
}

// ======================================================================
// psc_push_particles: subclass "1vb_4x4_cuda"

struct psc_push_particles_ops psc_push_particles_1vb_4x4_cuda_ops = {
  .name                  = "1vb_4x4_cuda",
  .push_mprts_yz         = psc_push_particles_1vb_4x4_cuda_push_mprts_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD,
};

// ======================================================================

DECLARE_CUDA(yz2x2_1vb);

static void
psc_push_particles_cuda_2x2_1vb_push_a_yz(struct psc_push_particles *push,
					  struct psc_particles *prts_base,
					  struct psc_fields *flds_base)
{
  cuda_push_partq_a(push, prts_base, flds_base, yz2x2_1vb_cuda_push_part_p2);
  cuda_push_partq_b(push, prts_base, flds_base, yz2x2_1vb_cuda_push_part_p3);
}

// ======================================================================
// psc_push_particles: subclass "1vb_2x2_cuda"

struct psc_push_particles_ops psc_push_particles_1vb_2x2_cuda_ops = {
  .name                  = "1vb_2x2_cuda",
  .push_a_yz             = psc_push_particles_cuda_2x2_1vb_push_a_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_2X2X2 | MP_NO_CHECKERBOARD,
};

// ======================================================================

DECLARE_CUDA(yz8x8_1vb);

static void
psc_push_particles_cuda_8x8_1vb_push_a_yz(struct psc_push_particles *push,
					  struct psc_particles *prts_base,
					  struct psc_fields *flds_base)
{
  cuda_push_partq_a(push, prts_base, flds_base, yz8x8_1vb_cuda_push_part_p2);
  cuda_push_partq_b(push, prts_base, flds_base, yz8x8_1vb_cuda_push_part_p3);
}

// ======================================================================
// psc_push_particles: subclass "1vb_8x8_cuda"

struct psc_push_particles_ops psc_push_particles_1vb_8x8_cuda_ops = {
  .name                  = "1vb_8x8_cuda",
  .push_a_yz             = psc_push_particles_cuda_8x8_1vb_push_a_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_8X8X8 | MP_NO_CHECKERBOARD,
};
