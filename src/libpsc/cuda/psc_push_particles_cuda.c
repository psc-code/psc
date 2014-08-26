
#include "psc_push_particles_private.h"
#include "psc_cuda.h"
#include "psc_fields.h"

#include <mrc_profile.h>
#include <math.h>

// ======================================================================
// psc_push_particles: "1vb_4x4_cuda"

static void
psc_push_particles_1vb_4x4_cuda_push_mprts_yz(struct psc_push_particles *push,
					      struct psc_mparticles *mprts,
					      struct psc_mfields *mflds_base)
{
  // it's difficult to convert mprts because of the ordering constraints (?)
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", EX, EX + 6);
  yz4x4_1vb_cuda_push_mprts(mprts, mflds);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vb_4x4_cuda"

struct psc_push_particles_ops psc_push_particles_1vb_4x4_cuda_ops = {
  .name                  = "1vb_4x4_cuda",
  .push_mprts_yz         = psc_push_particles_1vb_4x4_cuda_push_mprts_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD,
};

// ======================================================================
// psc_push_particles: subclass "1vbec3d_2x2_cuda"

static void
psc_push_particles_1vbec3d_2x2_cuda_push_mprts_yz(struct psc_push_particles *push,
						  struct psc_mparticles *mprts,
						  struct psc_mfields *mflds_base)
{
  // it's difficult to convert mprts because of the ordering constraints (?)
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", EX, EX + 6);
  yz2x2_1vbec3d_cuda_push_mprts(mprts, mflds);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vb_2x2_cuda"

struct psc_push_particles_ops psc_push_particles_1vbec3d_2x2_cuda_ops = {
  .name                  = "1vbec3d_2x2_cuda",
  .push_mprts_yz         = psc_push_particles_1vbec3d_2x2_cuda_push_mprts_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_2X2X2 | MP_NO_CHECKERBOARD,
};

// ======================================================================
// psc_push_particles: subclass "1vbec3d_4x4_cuda"

static void
psc_push_particles_1vbec3d_4x4_cuda_push_mprts_yz(struct psc_push_particles *push,
						  struct psc_mparticles *mprts,
						  struct psc_mfields *mflds_base)
{
  // it's difficult to convert mprts because of the ordering constraints (?)
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", EX, EX + 6);
  yz4x4_1vbec3d_cuda_push_mprts(mprts, mflds);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vb_4x4_cuda"

struct psc_push_particles_ops psc_push_particles_1vbec3d_4x4_cuda_ops = {
  .name                  = "1vbec3d_4x4_cuda",
  .push_mprts_yz         = psc_push_particles_1vbec3d_4x4_cuda_push_mprts_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD,
};

// ======================================================================
// psc_push_particles: subclass "1vbec3d_8x8_cuda"

static void
psc_push_particles_1vbec3d_8x8_cuda_push_mprts_yz(struct psc_push_particles *push,
						  struct psc_mparticles *mprts,
						  struct psc_mfields *mflds_base)
{
  // it's difficult to convert mprts because of the ordering constraints (?)
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", EX, EX + 6);
  yz8x8_1vbec3d_cuda_push_mprts(mprts, mflds);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vb_8x8_cuda"

struct psc_push_particles_ops psc_push_particles_1vbec3d_8x8_cuda_ops = {
  .name                  = "1vbec3d_8x8_cuda",
  .push_mprts_yz         = psc_push_particles_1vbec3d_8x8_cuda_push_mprts_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_8X8X8 | MP_NO_CHECKERBOARD,
};

