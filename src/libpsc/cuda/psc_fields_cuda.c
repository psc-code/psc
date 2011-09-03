
#include "psc.h"
#include "psc_fields_cuda.h"

#include <mrc_profile.h>

static void
fields_cuda_alloc(fields_cuda_t *pf, int ib[3], int ie[3], int nr_comp)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->h_flds = calloc(nr_comp * size, sizeof(*pf->h_flds));
  pf->name = calloc(nr_comp, sizeof(*pf->name));
}

static void
fields_cuda_free(fields_cuda_t *pf)
{
  free(pf->h_flds);
  for (int m = 0; m < pf->nr_comp; m++) {
    free(pf->name[m]);
  }
  free(pf->name);
}

void
psc_mfields_cuda_get_from(mfields_cuda_t *flds, int mb, int me, void *_flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fields_cuda_get", 1., 0, 0);
  }
  prof_start(pr);

  mfields_base_t *flds_base = _flds_base;
  flds->f = calloc(ppsc->nr_patches, sizeof(*flds->f));
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf = &flds->f[p];
    struct psc_patch *patch = &ppsc->patch[p];
    int ilg[3] = { -ppsc->ibn[0], -ppsc->ibn[1], -ppsc->ibn[2] };
    int ihg[3] = { patch->ldims[0] + ppsc->ibn[0],
		   patch->ldims[1] + ppsc->ibn[1],
		   patch->ldims[2] + ppsc->ibn[2] };
    fields_cuda_alloc(pf, ilg, ihg, 12);

    fields_base_t *pf_base = &flds_base->f[p];
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_CUDA(pf, m, jx,jy,jz) = F3_BASE(pf_base, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
    __fields_cuda_get(pf);
  }

  prof_stop(pr);
}

void
psc_mfields_cuda_put_to(mfields_cuda_t *flds, int mb, int me, void *_flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fields_cuda_put", 1., 0, 0);
  }
  prof_start(pr);

  mfields_base_t *flds_base = _flds_base;
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf = &flds->f[p];
    fields_base_t *pf_base = &flds_base->f[p];
    __fields_cuda_put(pf);
  
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_BASE(pf_base, m, jx,jy,jz) = F3_CUDA(pf, m, jx,jy,jz);
      }
    } foreach_3d_g_end;

    fields_cuda_free(pf);
  }
  
  free(flds->f);
  flds->f = NULL;

  prof_stop(pr);
}

