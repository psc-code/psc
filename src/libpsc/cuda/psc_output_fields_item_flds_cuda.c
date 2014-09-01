
#include "psc_output_fields_item_private.h"

#include "psc_fields_as_single.h"

// ======================================================================

static void
calc_dive_nc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	     struct psc_particles *prts, struct psc_fields *f_base)
{
  assert(ppsc->domain.gdims[0] == 1);

  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 3);
  struct psc_fields *f = psc_fields_get_as(f_base, FIELDS_TYPE, 0, 0);

  float dy = ppsc->patch[f->p].dx[1];
  float dz = ppsc->patch[f->p].dx[2];

  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = 
      ((F3(flds, EY, ix,iy,iz) - F3(flds, EY, ix,iy-1,iz)) / dy +
       (F3(flds, EZ, ix,iy,iz) - F3(flds, EZ, ix,iy,iz-1)) / dz);
  } foreach_3d_end;

  psc_fields_put_as(flds, flds_base, 0, 0);
  psc_fields_put_as(f, f_base, 0, 1);
}

struct psc_output_fields_item_ops psc_output_fields_item_dive_cuda_ops = {
  .name      = "dive_cuda",
  .nr_comp   = 1,
  .fld_names = { "dive" },
  .run       = calc_dive_nc,
};

