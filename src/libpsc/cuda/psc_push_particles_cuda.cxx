
#include "cuda_iface.h"
#include "psc_push_particles_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#include <string.h>

struct CurrentGlobal : std::true_type {};
struct CurrentShared : std::false_type {};

// ======================================================================
// psc_push_particles: "1vb_4x4_cuda"

static void
psc_push_particles_1vb_4x4_cuda_push_mprts_yz(struct psc_push_particles *push,
					      struct psc_mparticles *mprts,
					      struct psc_mfields *mflds_base)
{
  // it's difficult to convert mprts because of the ordering constraints (?)
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);

  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(EX, EX + 6);
  struct cuda_mparticles *cmprts = mparticles_cuda_t(mprts)->cmprts();
  int bs[3] = { BS::x::value, BS::y::value, BS::z::value };
  cuda_push_mprts_yz(cmprts, mf->cmflds, bs, false, false, false);
  mf.put_as(mflds_base, JXI, JXI + 3);
}

template<typename Current>
static void psc_push_particles_1vbec3d_cuda_push_mprts_yz(struct psc_push_particles *push,
							  struct psc_mparticles *mprts,
							  struct psc_mfields *mflds_base) \
{
  /* it's difficult to convert mprts due to ordering constraints (?) */
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);

  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(EX, EX + 6);
  struct cuda_mparticles *cmprts = mparticles_cuda_t(mprts)->cmprts();
  struct cuda_mfields *cmflds = mf->cmflds;
  int bs[3] = { BS::x::value, BS::y::value, BS::z::value };
  cuda_push_mprts_yz(cmprts, cmflds, bs, true, true, Current::value);
  mf.put_as(mflds_base, JXI, JXI + 3);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vb_cuda"

struct psc_push_particles_ops_cuda : psc_push_particles_ops {
  psc_push_particles_ops_cuda() {
    name                  = "1vb_cuda";
    push_mprts_yz         = psc_push_particles_1vb_4x4_cuda_push_mprts_yz;
    mp_flags              = MP_BLOCKSIZE_4X4X4;
  }
} psc_push_particles_1vb_cuda_ops;

#define MAKE_1VBEC3D_YZ(MP_BS, MEM, CURRENT)		\
									\
  /* --------------------------------------------------------------- */	\
  /* psc_push_particles: subclass "1vbec3d_BYxBZ_MEM_cuda"           */ \
									\
  struct psc_push_particles_ops_1vbec3d ##MEM## _cuda :			\
    psc_push_particles_ops {						\
    psc_push_particles_ops_1vbec3d ##MEM## _cuda() {			\
      name                  = "1vbec" #MEM "_cuda";			\
      push_mprts_yz         = psc_push_particles_1vbec3d_cuda_push_mprts_yz<CURRENT>; \
      mp_flags              = MP_BS;					\
    }									\
  } psc_push_particles_1vbec3d ##MEM## _cuda_ops;


MAKE_1VBEC3D_YZ(MP_BLOCKSIZE_4X4X4, , CurrentShared);
MAKE_1VBEC3D_YZ(MP_BLOCKSIZE_4X4X4, _gmem, CurrentGlobal);

