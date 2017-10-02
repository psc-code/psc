
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"

#include <mrc_io.h>

// FIXME: checkpointing won't properly restore state
// FIXME: if they subclass creates objects, it'd be cleaner to have them
// be part of the subclass

// ----------------------------------------------------------------------
// fld_create
//
// FIXME, should be consolidated with psc_checks.c, and probably other places

static struct psc_mfields *
fld_create(struct psc *psc, const char *name)
{
  struct psc_mfields *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(fld, FIELDS_TYPE);
  psc_mfields_set_param_obj(fld, "domain", psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_set_param_int(fld, "nr_fields", 1);
  psc_mfields_setup(fld);
  psc_mfields_set_comp_name(fld, 0, name);

  return fld;
}

// ----------------------------------------------------------------------
// psc_marder_sub_setup

static void
psc_marder_sub_setup(struct psc_marder *marder)
{
  marder->div_e = fld_create(ppsc, "div_E");
  marder->rho = fld_create(ppsc, "rho");

  marder->bnd = psc_bnd_create(psc_marder_comm(marder));
  psc_bnd_set_name(marder->bnd, "marder_bnd");
  psc_bnd_set_type(marder->bnd, FIELDS_TYPE);
  psc_bnd_set_psc(marder->bnd, ppsc);
  psc_bnd_setup(marder->bnd);

  // FIXME, output_fields should be taking care of their own psc_bnd?
  marder->item_div_e = psc_output_fields_item_create(psc_comm(ppsc));
  psc_output_fields_item_set_type(marder->item_div_e, "dive_" FIELDS_TYPE);
  psc_output_fields_item_set_psc_bnd(marder->item_div_e, marder->bnd);
  psc_output_fields_item_setup(marder->item_div_e);

  marder->item_rho = psc_output_fields_item_create(psc_comm(ppsc));
  psc_output_fields_item_set_type(marder->item_rho, "rho_1st_nc_" PARTICLE_TYPE);
  psc_output_fields_item_set_psc_bnd(marder->item_rho, marder->bnd);
  psc_output_fields_item_setup(marder->item_rho);

  if (marder->dump) {
    struct mrc_io *io = mrc_io_create(psc_comm(ppsc));
    mrc_io_set_type(io, "xdmf_collective");
    mrc_io_set_name(io, "mrc_io_marder");
    mrc_io_set_param_string(io, "basename", "marder");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);

    marder->io = io;
  }
}

// ----------------------------------------------------------------------
// psc_marder_sub_destroy

static void
psc_marder_sub_destroy(struct psc_marder *marder)
{
  psc_mfields_destroy(marder->div_e);
  psc_mfields_destroy(marder->rho);

  psc_output_fields_item_destroy(marder->item_div_e);
  psc_output_fields_item_destroy(marder->item_rho);

  psc_bnd_destroy(marder->bnd);

  if (marder->dump) {
    mrc_io_destroy(marder->io);
  }
}

// ----------------------------------------------------------------------
// psc_marder_sub_correct_patch
//
// Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

#define psc_foreach_3d_more(psc, p, ix, iy, iz, l, r) {			\
  int __ilo[3] = { -l[0], -l[1], -l[2] };					\
  int __ihi[3] = { psc->patch[p].ldims[0] + r[0],				\
		   psc->patch[p].ldims[1] + r[1],				\
		   psc->patch[p].ldims[2] + r[2] };				\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {			\
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {			\
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define psc_foreach_3d_more_end				\
  } } }

static void
psc_marder_sub_correct_patch(struct psc_marder *marder, fields_t flds, fields_t f, int p)
{
  define_dxdydz(dx, dy, dz);

  // FIXME: how to choose diffusion parameter properly?
  //double deltax = ppsc->patch[p].dx[0];
  double deltay = ppsc->patch[p].dx[1]; // FIXME double/float
  double deltaz = ppsc->patch[p].dx[2];
  double inv_sum = 0.;
  int nr_levels;
  mrc_domain_get_nr_levels(ppsc->mrc_domain, &nr_levels);
  for (int d=0;d<3;d++) {
    if (ppsc->domain.gdims[d] > 1) {
      inv_sum += 1. / sqr(ppsc->patch[p].dx[d] / (1 << (nr_levels - 1)));
    }
  }
  double diffusion_max = 1. / 2. / (.5 * ppsc->dt) / inv_sum;
  double diffusion     = diffusion_max * marder->diffusion;

  int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
  int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
  for (int d = 0; d < 3; d++) {
   if (ppsc->domain.bnd_fld_lo[d] == BND_FLD_CONDUCTING_WALL && ppsc->patch[p].off[d] == 0) {
    l_cc[d] = -1;
    l_nc[d] = -1;
   }
   if (ppsc->domain.bnd_fld_hi[d] == BND_FLD_CONDUCTING_WALL && ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
    r_cc[d] = -1;
    r_nc[d] = 0;
   }
  }

#if 0
  psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
    // FIXME: F3 correct?
    _F3(flds, EX, ix,iy,iz) += 
      (_F3(f, DIVE_MARDER, ix+dx,iy,iz) - _F3(f, DIVE_MARDER, ix,iy,iz))
      * .5 * ppsc->dt * diffusion / deltax;
    _F3(flds, EY, ix,iy,iz) += 
      (_F3(f, DIVE_MARDER, ix,iy+dy,iz) - _F3(f, DIVE_MARDER, ix,iy,iz))
      * .5 * ppsc->dt * diffusion / deltay;
    _F3(flds, EZ, ix,iy,iz) += 
      (_F3(f, DIVE_MARDER, ix,iy,iz+dz) - _F3(f, DIVE_MARDER, ix,iy,iz))
      * .5 * ppsc->dt * diffusion / deltaz;
  } psc_foreach_3d_more_end;
#endif

  assert(ppsc->domain.gdims[0] == 1);

  {
    int l[3] = { l_nc[0], l_cc[1], l_nc[2] };
    int r[3] = { r_nc[0], r_cc[1], r_nc[2] };
    psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
      _F3(flds, EY, ix,iy,iz) += 
	(_F3(f, 0, ix,iy+dy,iz) - _F3(f, 0, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltay;
    } psc_foreach_3d_more_end;
  }

  {
    int l[3] = { l_nc[0], l_nc[1], l_cc[2] };
    int r[3] = { r_nc[0], r_nc[1], r_cc[2] };
    psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
      _F3(flds, EZ, ix,iy,iz) += 
	(_F3(f, 0, ix,iy,iz+dz) - _F3(f, 0, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltaz;
    } psc_foreach_3d_more_end;
  }
}

#undef psc_foreach_3d_more
#undef psc_foreach_3d_more_end

// ----------------------------------------------------------------------
// psc_marder_sub_correct

static void
psc_marder_sub_correct(struct psc_marder *marder, struct psc_mfields *mflds_base,
		       struct psc_mfields *div_e)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, EX + 3);

  for (int p = 0; p < div_e->nr_patches; p++) {
    psc_marder_sub_correct_patch(marder, fields_t_mflds(mflds, p),
				 fields_t_mflds(div_e, p), p);
  }

  psc_mfields_put_as(mflds, mflds_base, EX, EX + 3);
}

// ======================================================================
// psc_marder: subclass FIELDS_TYPE

struct psc_marder_ops psc_marder_sub_ops = {
  .name                  = FIELDS_TYPE,
  .setup                 = psc_marder_sub_setup,
  .destroy               = psc_marder_sub_destroy,
  .correct               = psc_marder_sub_correct,
};

