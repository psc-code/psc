
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_as_single.h"
#include "psc_particles_as_single.h"

#include <mrc_io.h>

#include <math.h>

// ----------------------------------------------------------------------
// marder_calc_aid_fields

static void
marder_calc_aid_fields(struct psc_marder *marder, 
		       struct psc_mfields *flds, struct psc_mparticles *particles,
		       struct psc_mfields *div_e, struct psc_mfields *rho)
{
  psc_output_fields_item_run(marder->item_div_e, flds, particles, div_e); // FIXME, should accept NULL for particles
  psc_output_fields_item_run(marder->item_rho, flds, particles, rho);
  
  if (marder->dump) {
    mrc_io_open(marder->io, "w", ppsc->timestep, ppsc->timestep * ppsc->dt);
    psc_mfields_write_as_mrc_fld(rho, marder->io);
    psc_mfields_write_as_mrc_fld(div_e, marder->io);
    mrc_io_close(marder->io);
  }

  psc_mfields_axpy_comp(div_e, 0, -1., rho, 0);
  // FIXME, why is this necessary?
  psc_bnd_fill_ghosts(marder->bnd, div_e, 0, 1);
}


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

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

// ======================================================================
// Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

static void
do_marder_correction(struct psc_marder *marder,
		     struct psc_fields *flds_base, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);

  // FIXME: how to choose diffusion parameter properly?
  //double deltax = ppsc->patch[f->p].dx[0];
  double deltay = ppsc->patch[f->p].dx[1];
  double deltaz = ppsc->patch[f->p].dx[2];
  double inv_sum = 0.;
  int nr_levels;
  mrc_domain_get_nr_levels(ppsc->mrc_domain, &nr_levels);
  for (int d=0;d<3;d++) {
    if (ppsc->domain.gdims[d] > 1) {
      inv_sum += 1. / sqr(ppsc->patch[f->p].dx[d] / (1 << (nr_levels - 1)));
    }
  }
  double diffusion_max = 1. / 2. / (.5 * ppsc->dt) / inv_sum;
  double diffusion     = diffusion_max * marder->diffusion;

  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 3);

  int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
  int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
  for (int d = 0; d < 3; d++) {
   if (ppsc->domain.bnd_fld_lo[d] == BND_FLD_CONDUCTING_WALL && ppsc->patch[flds->p].off[d] == 0) {
    l_cc[d] = -1;
    l_nc[d] = -1;
   }
   if (ppsc->domain.bnd_fld_hi[d] == BND_FLD_CONDUCTING_WALL && ppsc->patch[flds->p].off[d] + ppsc->patch[flds->p].ldims[d] == ppsc->domain.gdims[d]) {
    r_cc[d] = -1;
    r_nc[d] = 0;
   }
  }

#if 0
  psc_foreach_3d_more(ppsc, f->p, ix, iy, iz, l, r) {
    // FIXME: F3 correct?
    F3(flds, EX, ix,iy,iz) += 
      (F3(f, DIVE_MARDER, ix+dx,iy,iz) - F3(f, DIVE_MARDER, ix,iy,iz))
      * .5 * ppsc->dt * diffusion / deltax;
    F3(flds, EY, ix,iy,iz) += 
      (F3(f, DIVE_MARDER, ix,iy+dy,iz) - F3(f, DIVE_MARDER, ix,iy,iz))
      * .5 * ppsc->dt * diffusion / deltay;
    F3(flds, EZ, ix,iy,iz) += 
      (F3(f, DIVE_MARDER, ix,iy,iz+dz) - F3(f, DIVE_MARDER, ix,iy,iz))
      * .5 * ppsc->dt * diffusion / deltaz;
  } psc_foreach_3d_more_end;
#endif

  assert(ppsc->domain.gdims[0] == 1);

  {
    int l[3] = { l_nc[0], l_cc[1], l_nc[2] };
    int r[3] = { r_nc[0], r_cc[1], r_nc[2] };
    psc_foreach_3d_more(ppsc, f->p, ix, iy, iz, l, r) {
      F3(flds, EY, ix,iy,iz) += 
	(F3(f, 0, ix,iy+dy,iz) - F3(f, 0, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltay;
    } psc_foreach_3d_more_end;
  }

  {
    int l[3] = { l_nc[0], l_nc[1], l_cc[2] };
    int r[3] = { r_nc[0], r_nc[1], r_cc[2] };
    psc_foreach_3d_more(ppsc, f->p, ix, iy, iz, l, r) {
      F3(flds, EZ, ix,iy,iz) += 
	(F3(f, 0, ix,iy,iz+dz) - F3(f, 0, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltaz;
    } psc_foreach_3d_more_end;
  }

  psc_fields_put_as(flds, flds_base, EX, EX + 3);
}

static void
marder_correction_run(struct psc_marder *marder, struct psc_mfields *flds,
		      struct psc_mfields *f)
{
  for (int p = 0; p < f->nr_patches; p++) {
    do_marder_correction(marder, psc_mfields_get_patch(flds, p),
			 psc_mfields_get_patch(f, p));
  }
  psc_bnd_fill_ghosts(ppsc->bnd, flds, EX, EX+3);
}

// ----------------------------------------------------------------------
// psc_marder_run
//
// On ghost cells:
// It is possible (variant = 1) that ghost cells are set before this is called
// and the subsequent code expects ghost cells to still be set on return.
// We're calling fill_ghosts at the end of each iteration, so that's fine.
// However, for variant = 0, ghost cells aren't set on entry, and they're not
// expected to be set on return (though we do that, anyway.)

void
psc_marder_run(struct psc_marder *marder, 
	       struct psc_mfields *flds, struct psc_mparticles *particles)
{
  if (marder->every_step < 0 || ppsc->timestep % marder->every_step != 0) 
   return;

  // need to fill ghost cells first (should be unnecessary with only variant 1) FIXME
  psc_bnd_fill_ghosts(ppsc->bnd, flds, EX, EX+3);

  for (int i = 0; i < marder->loop; i++) {
    marder_calc_aid_fields(marder, flds, particles, marder->div_e, marder->rho);
    marder_correction_run(marder, flds, marder->div_e);
  }
}

// ----------------------------------------------------------------------
// psc_marder_init

static void
psc_marder_init()
{
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_single_ops);
}

// ----------------------------------------------------------------------
// psc_marder description

#define VAR(x) (void *)offsetof(struct psc_marder, x)
static struct param psc_marder_descr[] = {
  { "every_step"       , VAR(every_step)       , PARAM_INT(-1),     },
  { "diffusion"        , VAR(diffusion)        , PARAM_DOUBLE(0.9), },
  { "loop"             , VAR(loop)             , PARAM_INT(1),      },
  { "dump"             , VAR(dump)             , PARAM_BOOL(false), },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_marder class description

struct mrc_class_psc_marder mrc_class_psc_marder = {
  .name             = "psc_marder",
  .size             = sizeof(struct psc_marder),
  .param_descr      = psc_marder_descr,
  .init             = psc_marder_init,
};

#undef psc_foreach_3d_more
#undef psc_foreach_3d_more_end
