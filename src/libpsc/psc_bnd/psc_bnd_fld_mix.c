
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"
#include "psc_fields_single.h"
#include "psc_fields_cuda.h"
#include "psc.h"

#include <mrc_ddc.h>
#include <mrc_profile.h>

// ======================================================================
// ddc funcs

static void
copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct psc_mfields *flds = ctx;
  struct psc_fields *pf = psc_mfields_get_patch(flds, p);
  if (psc_fields_ops(pf) == &psc_fields_single_ops) {
    psc_bnd_fld_single_copy_to_buf(mb, me, p, ilo, ihi, _buf, ctx);
  } else if (psc_fields_ops(pf) == &psc_fields_cuda_ops) {
    psc_bnd_fld_cuda_copy_to_buf(mb, me, p, ilo, ihi, _buf, ctx);
  } else {
    assert(0);
  }
}

static void
add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct psc_mfields *flds = ctx;
  struct psc_fields *pf = psc_mfields_get_patch(flds, p);
  if (psc_fields_ops(pf) == &psc_fields_single_ops) {
    psc_bnd_fld_single_add_from_buf(mb, me, p, ilo, ihi, _buf, ctx);
  } else if (psc_fields_ops(pf) == &psc_fields_cuda_ops) {
    psc_bnd_fld_cuda_add_from_buf(mb, me, p, ilo, ihi, _buf, ctx);
  } else {
    assert(0);
  }
}

static void
copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct psc_mfields *flds = ctx;
  struct psc_fields *pf = psc_mfields_get_patch(flds, p);
  if (psc_fields_ops(pf) == &psc_fields_single_ops) {
    psc_bnd_fld_single_copy_from_buf(mb, me, p, ilo, ihi, _buf, ctx);
  } else if (psc_fields_ops(pf) == &psc_fields_cuda_ops) {
    psc_bnd_fld_cuda_copy_from_buf(mb, me, p, ilo, ihi, _buf, ctx);
  } else {
    assert(0);
  }
}

static struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = copy_to_buf,
  .copy_from_buf = copy_from_buf,
  .add_from_buf  = add_from_buf,
};

// ----------------------------------------------------------------------
// psc_bnd_fld_mix_create

void
psc_bnd_fld_mix_create(struct psc_bnd *bnd)
{
  struct mrc_ddc *ddc = mrc_domain_create_ddc(bnd->psc->mrc_domain);
  mrc_ddc_set_funcs(ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(ddc, "ibn", bnd->psc->ibn);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 12);
  assert(sizeof(fields_single_real_t) == sizeof(fields_cuda_real_t));
  mrc_ddc_set_param_int(ddc, "size_of_type", sizeof(fields_single_real_t));
  mrc_ddc_setup(ddc);
  bnd->ddc = ddc;
}

// ----------------------------------------------------------------------
// get_ops

static inline struct psc_bnd_ops *
get_ops(struct psc_fields *pf)
{
  if (psc_fields_ops(pf) == &psc_fields_single_ops) {
    return &psc_bnd_single_ops;
  } else if (psc_fields_ops(pf) == &psc_fields_cuda_ops) {
    return &psc_bnd_cuda_ops;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_fld_mix_add_ghosts

void
psc_bnd_fld_mix_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base, int mb, int me)
{
  for (int p = 0; p < flds_base->nr_patches; p++) {
    struct psc_fields *pf = psc_mfields_get_patch(flds_base, p);
    struct psc_bnd_ops *ops = get_ops(pf);
    if (ops->add_ghosts_prep) {
      ops->add_ghosts_prep(bnd, pf, mb, me);
    }
  }

  mrc_ddc_add_ghosts(bnd->ddc, mb, me, flds_base);

  for (int p = 0; p < flds_base->nr_patches; p++) {
    struct psc_fields *pf = psc_mfields_get_patch(flds_base, p);
    struct psc_bnd_ops *ops = get_ops(pf);
    if (ops->add_ghosts_post) {
      ops->add_ghosts_post(bnd, pf, mb, me);
    }
  }
}

// ----------------------------------------------------------------------
// psc_bnd_fld_mix_fill_ghosts

void
psc_bnd_fld_mix_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base, int mb, int me)
{
  for (int p = 0; p < flds_base->nr_patches; p++) {
    struct psc_fields *pf = psc_mfields_get_patch(flds_base, p);
    struct psc_bnd_ops *ops = get_ops(pf);
    if (ops->fill_ghosts_prep) {
      ops->fill_ghosts_prep(bnd, pf, mb, me);
    }
  }

  mrc_ddc_fill_ghosts(bnd->ddc, mb, me, flds_base);

  for (int p = 0; p < flds_base->nr_patches; p++) {
    struct psc_fields *pf = psc_mfields_get_patch(flds_base, p);
    struct psc_bnd_ops *ops = get_ops(pf);
    if (ops->fill_ghosts_post) {
      ops->fill_ghosts_post(bnd, pf, mb, me);
    }
  }
}
