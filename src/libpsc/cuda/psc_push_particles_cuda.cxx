
#include "cuda_iface.h"
#include "psc_push_particles_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#include <string.h>

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
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;
  int bs[3] = { 1, 4, 4 };
  cuda_push_mprts_yz(cmprts, cmflds, bs, false, false, false);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vb_4x4_cuda"

struct psc_push_particles_ops psc_push_particles_1vb_4x4_cuda_ops = {
  .name                  = "1vb_4x4_cuda",
  .push_mprts_yz         = psc_push_particles_1vb_4x4_cuda_push_mprts_yz,
  .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD,
};

#define MAKE_1VBEC3D_YZ(BY, BZ, MP_BS, CURRMEM_GLOBAL, MEM)		\
									\
  /* =============================================================== */ \
  /* psc_push_particles: subclass "1vbec3d_BYxBZ_MEM_cuda"    	     */	\
  									\
  static void								\
  psc_push_particles_1vbec3d_ ##BY## x ##BZ## MEM## _cuda_push_mprts_yz \
  (struct psc_push_particles *push, struct psc_mparticles *mprts,	\
   struct psc_mfields *mflds_base)					\
  {									\
  /* it's difficult to convert mprts due to ordering constraints (?) */	\
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);		\
									\
  struct psc_mfields *mflds =						\
    psc_mfields_get_as(mflds_base, "cuda", EX, EX + 6);			\
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;	\
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;	\
  int bs[3] = { 1, BY, BZ };						\
  cuda_push_mprts_yz(cmprts, cmflds, bs, true, true, CURRMEM_GLOBAL);	\
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);			\
  }									\
									\
  /* --------------------------------------------------------------- */	\
  /* psc_push_particles: subclass "1vbec3d_BYxBZ_MEM_cuda"           */ \
									\
  struct psc_push_particles_ops						\
  psc_push_particles_1vbec3d_ ##BY## x ##BZ## MEM## _cuda_ops = {	\
    .name                  = "1vbec3d_" #BY "x" #BZ #MEM "_cuda",	\
    .push_mprts_yz         = psc_push_particles_1vbec3d_ ##BY## x ##BZ## MEM## _cuda_push_mprts_yz, \
    .mp_flags              = MP_NEED_BLOCK_OFFSETS | MP_BS | MP_NO_CHECKERBOARD, \
  };

MAKE_1VBEC3D_YZ(2, 2, MP_BLOCKSIZE_2X2X2, false, );
MAKE_1VBEC3D_YZ(4, 4, MP_BLOCKSIZE_4X4X4, false, );
MAKE_1VBEC3D_YZ(8, 8, MP_BLOCKSIZE_8X8X8, false, );
MAKE_1VBEC3D_YZ(2, 2, MP_BLOCKSIZE_2X2X2, true, _gmem);
MAKE_1VBEC3D_YZ(4, 4, MP_BLOCKSIZE_4X4X4, true, _gmem);
MAKE_1VBEC3D_YZ(8, 8, MP_BLOCKSIZE_8X8X8, true, _gmem);

