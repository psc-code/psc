
#include "psc_debug.h"

#ifdef F3_CACHE

#include "psc_fields_single.h"

static struct psc_fields *
cache_fields_from_em(fields_t *pf)
{
  struct psc_fields *fld = psc_fields_create(psc_fields_comm(pf));
  psc_fields_set_type(fld, F3_CACHE_TYPE);
  // FIXME, can do -1 .. 1?
  psc_fields_set_param_int3(fld, "ib", (int[3]) { 0, -2, -2 });
  psc_fields_set_param_int3(fld, "im", (int[3]) { 1,
	pf->im[1] + 2 * pf->ib[1] + 4, pf->im[2] + 2 * pf->ib[2] + 4});
  psc_fields_set_param_int(fld, "nr_comp", 9); // JX .. HZ
  psc_fields_set_param_int(fld, "p", pf->p);
  psc_fields_setup(fld);
  for (int iz = fld->ib[2]; iz < fld->ib[2] + fld->im[2]; iz++) {
    for (int iy = fld->ib[1]; iy < fld->ib[1] + fld->im[1]; iy++) {
      F3_CACHE(fld, EX, 0,iy,iz) = F3(pf, EX, 0,iy,iz);
      F3_CACHE(fld, EY, 0,iy,iz) = F3(pf, EY, 0,iy,iz);
      F3_CACHE(fld, EZ, 0,iy,iz) = F3(pf, EZ, 0,iy,iz);
      F3_CACHE(fld, HX, 0,iy,iz) = F3(pf, HX, 0,iy,iz);
      F3_CACHE(fld, HY, 0,iy,iz) = F3(pf, HY, 0,iy,iz);
      F3_CACHE(fld, HZ, 0,iy,iz) = F3(pf, HZ, 0,iy,iz);
    }
  }
  return fld;
}

static void _mrc_unused
cache_fields_to_j(struct psc_fields *fld, fields_t *pf)
{
  for (int iz = fld->ib[2]; iz < fld->ib[2] + fld->im[2]; iz++) {
    for (int iy = fld->ib[1]; iy < fld->ib[1] + fld->im[1]; iy++) {
      F3(pf, JXI, 0,iy,iz) += F3_CACHE(fld, JXI, 0,iy,iz);
      F3(pf, JYI, 0,iy,iz) += F3_CACHE(fld, JYI, 0,iy,iz);
      F3(pf, JZI, 0,iy,iz) += F3_CACHE(fld, JZI, 0,iy,iz);
    }
  }
}

#endif

