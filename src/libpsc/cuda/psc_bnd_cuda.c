
#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "psc_bnd_private.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>
#include <string.h>

// OPT opportunities:
// general: make fields really 2D in CUDA and in bnd xchg
// reduce packed buffers even more, e.g., add only on cuda side finally

// FIXME
EXTERN_C void __fields_cuda_from_device_inside(struct psc_fields *pf, struct cuda_fields_ctx *cf, int mb, int me);
EXTERN_C void __fields_cuda_to_device_outside(struct psc_fields *pf, struct cuda_fields_ctx *cf, int mb, int me);
EXTERN_C void __fields_cuda_to_device_inside(struct psc_fields *pf, struct cuda_fields_ctx *cf, int mb, int me);
void psc_bnd_cuda_xchg_setup(struct psc_bnd *bnd);
void psc_bnd_cuda_xchg_unsetup(struct psc_bnd *bnd);
void psc_bnd_cuda_xchg_exchange_particles(struct psc_bnd *bnd,
					  mparticles_base_t *particles_base);

static void
cuda_mfields_ctx_init(struct cuda_mfields_ctx *ctx, mfields_cuda_t *flds, int nr_fields)
{
  // FIXME, could be done once and cached
  ctx->cf = malloc(ppsc->nr_patches * sizeof(*ctx->cf));
  psc_foreach_patch(ppsc, p) {
    struct psc_fields *pf = psc_mfields_get_patch(flds, p);
    struct cuda_fields_ctx *cf = &ctx->cf[p];
    int sz = 1;
    for (int d = 0; d < 3; d++) {
      if (pf->im[d] == 1 - 2 * pf->ib[d]) { // only 1 non-ghost point
	cf->im[d] = 1;
	cf->ib[d] = 0;
      } else {
	cf->im[d] = pf->im[d];
	cf->ib[d] = pf->ib[d];
      }
      sz *= cf->im[d];
    }
    cf->arr = malloc(nr_fields * sz * sizeof(*cf->arr));
    cf->arr_off = ctx->cf[p].arr 
      - ((cf->ib[2] * cf->im[1] + cf->ib[1]) * cf->im[0] + cf->ib[0]);
  }
}

static void
cuda_mfields_ctx_free(struct cuda_mfields_ctx *ctx)
{
  psc_foreach_patch(ppsc, p) {
    free(ctx->cf[p].arr);
  }
  free(ctx->cf);
}

// ======================================================================
// ddc funcs

static void
copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct cuda_mfields_ctx *ctx = _ctx;
  struct cuda_fields_ctx *cf = &ctx->cf[p];
  fields_cuda_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf, m - mb, ix,iy,iz) = F3_CF(cf, m, ix,iy,iz);
	}
      }
    }
  }
}

static void
add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct cuda_mfields_ctx *ctx = _ctx;
  struct cuda_fields_ctx *cf = &ctx->cf[p];
  fields_cuda_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_CF(cf, m, ix,iy,iz) += MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

static void
copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct cuda_mfields_ctx *ctx = _ctx;
  struct cuda_fields_ctx *cf = &ctx->cf[p];
  fields_cuda_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_CF(cf, m, ix,iy,iz) = MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

static struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = copy_to_buf,
  .copy_from_buf = copy_from_buf,
  .add_from_buf  = add_from_buf,
};

// ----------------------------------------------------------------------
// psc_bnd_cuda_setup

static void
psc_bnd_cuda_setup(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);
  struct psc *psc = bnd->psc;

  bnd_cuda->ddc = mrc_domain_create_ddc(psc->mrc_domain);
  mrc_ddc_set_funcs(bnd_cuda->ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(bnd_cuda->ddc, "ibn", psc->ibn);
  mrc_ddc_set_param_int(bnd_cuda->ddc, "max_n_fields", 6);
  mrc_ddc_set_param_int(bnd_cuda->ddc, "size_of_type", sizeof(fields_cuda_real_t));
  mrc_ddc_setup(bnd_cuda->ddc);

  psc_bnd_cuda_xchg_setup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_unsetup

static void
psc_bnd_cuda_unsetup(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);

  mrc_ddc_destroy(bnd_cuda->ddc);
  psc_bnd_cuda_xchg_unsetup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_destroy

static void
psc_bnd_cuda_destroy(struct psc_bnd *bnd)
{
  psc_bnd_cuda_unsetup(bnd);
}

// ----------------------------------------------------------------------
// check_domain
//
// check if the underlying mrc_domain changed since setup(),
// which might happen, e.g., through rebalancing.
// In this case, do setup() over.

static void
check_domain(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);
  struct psc *psc = bnd->psc;

  struct mrc_domain *domain = mrc_ddc_get_domain(bnd_cuda->ddc);
  if (domain != psc->mrc_domain) {
    psc_bnd_cuda_unsetup(bnd);
    psc_bnd_setup(bnd);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_add_ghosts

static void
psc_bnd_cuda_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);
  check_domain(bnd);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);
  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_add_ghosts_periodic_yz(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_add_ghosts_periodic_z(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else {
    mfields_cuda_t *flds_cuda = psc_mfields_get_cuda(flds_base, mb, me);
    struct cuda_mfields_ctx ctx;
    cuda_mfields_ctx_init(&ctx, flds_cuda, me - mb);

    psc_foreach_patch(ppsc, p) {
      struct psc_fields *pf_cuda = psc_mfields_get_patch(flds_cuda, p);
      struct cuda_fields_ctx *cf = ctx.cf + p;
      __fields_cuda_from_device_inside(pf_cuda, cf, mb, me);
    }

    mrc_ddc_add_ghosts(bnd_cuda->ddc, 0, me - mb, &ctx);

    psc_foreach_patch(ppsc, p) {
      struct psc_fields *pf_cuda = psc_mfields_get_patch(flds_cuda, p);
      struct cuda_fields_ctx *cf = ctx.cf + p;
      __fields_cuda_to_device_inside(pf_cuda, cf, mb, me);
    }

    cuda_mfields_ctx_free(&ctx);
    psc_mfields_put_cuda(flds_cuda, flds_base, mb, me);
  }
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_fill_ghosts

static void
psc_bnd_cuda_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  struct psc_bnd_cuda *bnd_cuda = to_psc_bnd_cuda(bnd);
  check_domain(bnd);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_fill_ghosts", 1., 0, 0);
  }
  static int pr1, pr3, pr5;
  if (!pr1) {
    pr1 = prof_register("cuda_fill_ghosts_1", 1., 0, 0);
    pr3 = prof_register("cuda_fill_ghosts_3", 1., 0, 0);
    pr5 = prof_register("cuda_fill_ghosts_5", 1., 0, 0);
  }
  prof_start(pr);
  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_fill_ghosts_periodic_yz(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_fill_ghosts_periodic_z(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else {
    mfields_cuda_t *flds_cuda = psc_mfields_get_cuda(flds_base, mb, me);
    struct cuda_mfields_ctx ctx;
    cuda_mfields_ctx_init(&ctx, flds_cuda, me - mb);

    prof_start(pr1);
    psc_foreach_patch(ppsc, p) {
      struct psc_fields *pf_cuda = psc_mfields_get_patch(flds_cuda, p);
      struct cuda_fields_ctx *cf = ctx.cf + p;
      __fields_cuda_from_device_inside(pf_cuda, cf, mb, me);
    }
    prof_stop(pr1);

    prof_start(pr3);
    mrc_ddc_fill_ghosts(bnd_cuda->ddc, 0, me - mb, &ctx);
    prof_stop(pr3);

    prof_start(pr5);
    psc_foreach_patch(ppsc, p) {
      struct psc_fields *pf_cuda = psc_mfields_get_patch(flds_cuda, p);
      struct cuda_fields_ctx *cf = ctx.cf + p;
      __fields_cuda_to_device_outside(pf_cuda, cf, mb, me);
    }
    prof_stop(pr5);

    cuda_mfields_ctx_free(&ctx);
    psc_mfields_put_cuda(flds_cuda, flds_base, mb, me);
  }
  prof_stop(pr);
}

static void
psc_bnd_cuda_exchange_particles(struct psc_bnd *bnd,
				mparticles_base_t *particles_base)
{
  check_domain(bnd);
  psc_bnd_cuda_xchg_exchange_particles(bnd, particles_base);
}

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops psc_bnd_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_bnd_cuda),
  .setup                 = psc_bnd_cuda_setup,
  .destroy               = psc_bnd_cuda_destroy,
  .add_ghosts            = psc_bnd_cuda_add_ghosts,
  .fill_ghosts           = psc_bnd_cuda_fill_ghosts,
  .exchange_particles    = psc_bnd_cuda_exchange_particles,
};

