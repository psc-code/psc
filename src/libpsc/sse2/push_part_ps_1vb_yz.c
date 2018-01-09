
#include "psc_push_particles_ps.h"

#include "psc_fields_single.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================

static void
fields_ip_alloc(fields_ip_t *pf, int ib[3], int ie[3], int nr_comp, int first_comp)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->first_comp = first_comp;
  pf->flds = calloc(nr_comp * size, sizeof(*pf->flds));
}

static void
fields_ip_free(fields_ip_t *pf)
{
  free(pf->flds);
}

#if SIMD_BITS == 0

static void
ip_fields_to_j(int p, fields_ip_t *fld, fields_t pf)
{
  struct psc_patch *patch = ppsc->patch + p;
  
  // FIXME, goes out of bounds
  for (int iz = -2; iz < patch->ldims[2] + 1; iz++) {
    for (int iy = -2; iy < patch->ldims[1] + 1; iy++) {
      fields_ip_real_t jx = F3_IP(fld, JXI, 0,iy,iz);
      _F3(pf, JXI, 0,iy  ,iz  ) += jx.f00;
      _F3(pf, JXI, 0,iy  ,iz+1) += jx.f01;
      _F3(pf, JXI, 0,iy+1,iz  ) += jx.f10;
      _F3(pf, JXI, 0,iy+1,iz+1) += jx.f11;

      fields_ip_real_t jyz = F3_IP(fld, JYI, 0,iy,iz);
      _F3(pf, JYI, 0,iy  ,iz  ) += jyz.f00;
      _F3(pf, JYI, 0,iy  ,iz+1) += jyz.f01;
      _F3(pf, JZI, 0,iy  ,iz  ) += jyz.f10;
      _F3(pf, JZI, 0,iy+1,iz  ) += jyz.f11;
    }
  }
}

#elif SIMD_BITS == 2

#define LOAD_JX_INC(fld, iy,iz) do {					\
    f00 = _mm_load_ps((float *) &F3_IP(fld, JXI, 0,iy  ,iz));		\
    f01 = _mm_load_ps((float *) &F3_IP(fld, JXI, 0,iy+1,iz));		\
    f10 = _mm_load_ps((float *) &F3_IP(fld, JXI, 0,iy+2,iz));		\
    f11 = _mm_load_ps((float *) &F3_IP(fld, JXI, 0,iy+3,iz));		\
    _MM_TRANSPOSE4_PS(f00, f01, f10, f11);				\
									\
    f10 = _mm_shuffle_ps(f10, f10, 0x93); /* rotate: 10 01 00 11 */	\
    f11 = _mm_shuffle_ps(f11, f11, 0x93);				\
    /* OPT, there may be a better way, e.g. mask with and, andnot */	\
    f00x = _mm_insert_ps(f10, f10, _MM_MK_INSERTPS_NDX(0, 0, 1));	\
    f40x = _mm_insert_ps(f10, f10, _MM_MK_INSERTPS_NDX(0, 0, 14));	\
    f01x = _mm_insert_ps(f11, f11, _MM_MK_INSERTPS_NDX(0, 0, 1));	\
    f41x = _mm_insert_ps(f11, f11, _MM_MK_INSERTPS_NDX(0, 0, 14));	\
  } while (0)

#define LOAD_FLD(jx00, pf, JXI, iy, iz) do {				\
    jx00 = (v4s) { _F3(pf, JXI, 0,iy  ,iz  ), _F3(pf, JXI, 0,iy+1,iz  ),	\
		   _F3(pf, JXI, 0,iy+2,iz  ), _F3(pf, JXI, 0,iy+3,iz  ) }; \
  } while (0)

#define SAVE_FLD(jx00, pf, JXI, iy, iz) do {			\
    for (int iiy = 0; iiy < 4; iiy++) {				\
      F3(pf, JXI, 0,iy+iiy  ,iz  ) = v4s_extract(jx00, iiy);	\
    }								\
  } while (0)							\
      
static void
ip_fields_to_j(int p, fields_ip_t *fld, fields_t pf)
{
  struct psc_patch *patch = ppsc->patch + p;

  // FIXME, goes out of bounds
  assert((patch->ldims[1] & 3) == 0);
  for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
    v4s f00, f01, f10, f11;
    v4s f00x, f40x, f01x, f41x;
    v4s jx00, jx01;

    int iy;
    iy = -2; {
      LOAD_FLD(jx00, pf, JXI, iy, iz  );
      LOAD_FLD(jx01, pf, JXI, iy, iz+1);
    }
    for (; iy < patch->ldims[1] - 2; iy += 4) {
      LOAD_JX_INC(fld, iy, iz);
      jx00 += f00 + f00x;
      jx01 += f01 + f01x;
      SAVE_FLD(jx00, pf, JXI, iy, iz  );
      SAVE_FLD(jx01, pf, JXI, iy, iz+1);
      LOAD_FLD(jx00, pf, JXI, iy+4,iz  );
      LOAD_FLD(jx01, pf, JXI, iy+4,iz+1);
      jx00 += f40x;
      jx01 += f41x;
    }
    for (; iy < patch->ldims[1] + 2; iy += 4) {
      LOAD_JX_INC(fld, iy, iz);
      jx00 += f00 + f00x;
      jx01 += f01 + f01x;
      SAVE_FLD(jx00, pf, JXI, iy, iz  );
      SAVE_FLD(jx01, pf, JXI, iy, iz+1);
    }
  }

  // FIXME, -> SSE2
  for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
    for (int iy = -2; iy < patch->ldims[1] + 2; iy++) {
      fields_ip_real_t jyz = F3_IP(fld, JYI, 0,iy,iz);
      _F3(pf, JYI, 0,iy  ,iz  ) += jyz.f00;
      _F3(pf, JYI, 0,iy  ,iz+1) += jyz.f01;
      _F3(pf, JZI, 0,iy  ,iz  ) += jyz.f10;
      _xF3(pf, JZI, 0,iy+1,iz  ) += jyz.f11;
    }
  }
}

#endif

// ======================================================================

#if SIMD_BITS == 0

static void _mrc_unused
ip_fields_from_em(int p, fields_ip_t *fld, fields_t pf)
{
  struct psc_patch *patch = ppsc->patch + p;

  for (int m = EX; m <= HZ; m++) {
    for (int iz = -2; iz < patch->ldims[2] + 1; iz++) {
      for (int iy = -2; iy < patch->ldims[1] + 1; iy++) {
	float f00 = _F3(pf, m, 0,iy  ,iz  );
	float f01 = _F3(pf, m, 0,iy  ,iz+1);
	float f10 = _F3(pf, m, 0,iy+1,iz  );
	float f11 = _F3(pf, m, 0,iy+1,iz+1);
	F3_IP(fld, m, 0,iy,iz).f00 = f00;
	F3_IP(fld, m, 0,iy,iz).f01 = f01 - f00;
	F3_IP(fld, m, 0,iy,iz).f10 = f10 - f00;
	F3_IP(fld, m, 0,iy,iz).f11 = f00 + f11 - f01 - f10;
      }
    }
  }
}

#elif SIMD_BITS == 2

static void _mrc_unused
ip_fields_from_em(int p, fields_ip_t *fld, fields_t pf)
{
  struct psc_patch *patch = ppsc->patch + p;

  assert((patch->ldims[1] & 3) == 0);
  for (int m = EX; m <= HZ; m++) {
    for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
      for (int iy = -2; iy < patch->ldims[1] + 2; iy += 4) {
	v4s f00 = { _F3(pf, m, 0,iy  ,iz  ), _F3(pf, m, 0,iy+1,iz  ), _F3(pf, m, 0,iy+2,iz  ), _F3(pf, m, 0,iy+3,iz  ) };
	v4s f01 = { _F3(pf, m, 0,iy  ,iz+1), _F3(pf, m, 0,iy+1,iz+1), _F3(pf, m, 0,iy+2,iz+1), _F3(pf, m, 0,iy+3,iz+1) };
	v4s f10 = { _F3(pf, m, 0,iy+1,iz  ), _F3(pf, m, 0,iy+2,iz  ), _F3(pf, m, 0,iy+3,iz  ), _F3(pf, m, 0,iy+4,iz  ) };
	v4s f11 = { _F3(pf, m, 0,iy+1,iz+1), _F3(pf, m, 0,iy+2,iz+1), _F3(pf, m, 0,iy+3,iz+1), _F3(pf, m, 0,iy+4,iz+1) };
	
	f11 -= f01 + f10 - f00;
	f01 -= f00;
	f10 -= f00;
	_MM_TRANSPOSE4_PS(f00, f01, f10, f11);
	_mm_store_ps((float *) &F3_IP(fld, m, 0,iy  ,iz), f00);
	_mm_store_ps((float *) &F3_IP(fld, m, 0,iy+1,iz), f01);
	_mm_store_ps((float *) &F3_IP(fld, m, 0,iy+2,iz), f10);
	_mm_store_ps((float *) &F3_IP(fld, m, 0,iy+3,iz), f11);
      }
    }
  }
}

#endif

// ======================================================================

void
psc_push_particles_1vb_ps_push_a_yz(struct psc_push_particles *push,
				    struct psc_mparticles *_mprts,
				    struct psc_mfields *mflds,
				    int p)
{
  static int pr;
  if (!pr) {
    pr = prof_register("ps_1vb_push_yz", 1., 0, 0);
  }

  mparticles_t mprts = mparticles_t(_mprts);

  struct psc_particles_single *sngl = psc_particles_single(_prts);
  assert(!sngl->need_reorder);

  prof_start(pr);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  struct psc_patch *patch = ppsc->patch + p;
  fields_ip_t fld_ip;
  // FIXME, can do -1 .. 1?
  int ib[3] = { 0, -2, -2 };
  int ie[3] = { 1, patch->ldims[1] + 2, patch->ldims[2] + 2 };

  fields_ip_alloc(&fld_ip, ib, ie, 9, 0); // JXI .. HZ
  ip_fields_from_em(p, &fld_ip, flds);

  particle_range_t prts = mprts[p].range();
  unsigned int n_prts = prts.size();
#if 1
  sb2_ps_1vb_yz_pxx_jxyz(p, &fld_ip, prts, 0);
  sb0_ps_1vb_yz_pxx_jxyz(p, &fld_ip, prts, n_prts & ~3);
#else
  sb2_ps_1vb_yz_p(p, &fld_ip, prts, 0);
  sb0_ps_1vb_yz_p(p, &fld_ip, prts, n_prts & ~3);

  sb0_ps_1vb_yz_xx_jxyz(p, &fld_ip, prts, 0);
#endif

  ip_fields_to_j(p, &fld_ip, flds);
  fields_ip_free(&fld_ip);
  prof_stop(pr);
}

// ======================================================================

#define SIMD_BITS 0
#include "ps_common.c"

void
psc_push_particles_1vb_ps2_push_a_yz(struct psc_push_particles *push,
				     struct psc_mparticles *mprts,
				     struct psc_mfields *mflds, int p)
{
  static int pr;
  if (!pr) {
    pr = prof_register("ps_1vb_push_yz", 1., 0, 0);
  }

  struct psc_particles_single *sngl = psc_particles_single(prts);

  if (!sngl->need_reorder) {
    MHERE;
    for (int n = 0; n < prts->n_part; n++) {
      sngl->b_ids[n] = n;
    }
  }
  sngl->need_reorder = false;

  // swap in alt array
  particle_single_t *tmp = sngl->particles;
  sngl->particles = sngl->particles_alt;
  sngl->particles_alt = tmp;
  // actual reorder happens while pushing
  
  prof_start(pr);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  struct psc_patch *patch = ppsc->patch + prts->p;
  fields_ip_t fld_ip;
  // FIXME, can do -1 .. 1?
  int ib[3] = { 0, -2, -2 };
  int ie[3] = { 1, patch->ldims[1] + 2, patch->ldims[2] + 2 };

  fields_ip_alloc(&fld_ip, ib, ie, 9, 0); // JXI .. HZ
  ip_fields_from_em(prts->p, &fld_ip, flds);

  sb2_ps2_1vb_yz_pxx_jxyz(prts->p, &fld_ip, prts, 0);
  sb0_ps2_1vb_yz_pxx_jxyz(prts->p, &fld_ip, prts, prts->n_part & ~3);

  ip_fields_to_j(prts->p, &fld_ip, flds);
  fields_ip_free(&fld_ip);
  prof_stop(pr);
}

