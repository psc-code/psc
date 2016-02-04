#include "psc_output_fields_item_private.h"

#include "psc_fields_as_c.h"

// FIXME, we're assuming that the result fields are "c" type

// ======================================================================

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

// ======================================================================

#define JX_NC(ix,iy,iz) (.5f * (F3(flds, JXI,ix,iy,iz) + F3(flds, JXI,ix-dx,iy,iz)))
#define JY_NC(ix,iy,iz) (.5f * (F3(flds, JYI,ix,iy,iz) + F3(flds, JYI,ix,iy-dy,iz)))
#define JZ_NC(ix,iy,iz) (.5f * (F3(flds, JZI,ix,iy,iz) + F3(flds, JZI,ix,iy,iz-dz)))

static void
calc_j_nc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	  struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, JXI + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = JX_NC(ix,iy,iz);
    F3(f, 1, ix,iy,iz) = JY_NC(ix,iy,iz);
    F3(f, 2, ix,iy,iz) = JZ_NC(ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_j_nc_ops = {
  .name      = "j_nc",
  .nr_comp   = 3,
  .fld_names = { "jx_nc", "jy_nc", "jz_nc" },
  .run       = calc_j_nc,
};

// ======================================================================

#define JX_CC(ix,iy,iz) (.25f * (F3(flds, JXI, ix,iy,iz   ) + F3(flds, JXI, ix,iy+dy,iz   ) + \
				 F3(flds, JXI, ix,iy,iz+dz) + F3(flds, JXI, ix,iy+dy,iz+dz)))
#define JY_CC(ix,iy,iz) (.25f * (F3(flds, JYI, ix,iy,iz   ) + F3(flds, JYI, ix+dx,iy,iz   ) + \
				 F3(flds, JYI, ix,iy,iz+dz) + F3(flds, JYI, ix+dx,iy,iz+dz)))
#define JZ_CC(ix,iy,iz) (.25f * (F3(flds, JZI, ix,iy,iz   ) + F3(flds, JZI, ix+dx,iy   ,iz) + \
				 F3(flds, JZI, ix,iy+dy,iz) + F3(flds, JZI, ix+dx,iy+dy,iz)))

static void
calc_j(struct psc_output_fields_item *item, struct psc_fields *flds_base,
       struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, JXI + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = JX_CC(ix,iy,iz);
    F3(f, 1, ix,iy,iz) = JY_CC(ix,iy,iz);
    F3(f, 2, ix,iy,iz) = JZ_CC(ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_j_ops = {
  .name      = "j",
  .nr_comp   = 3,
  .fld_names = { "jx", "jy", "jz" },
  .run       = calc_j,
};

// ======================================================================

static void
calc_j_ec(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	  struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, JXI + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = F3(flds, JXI, ix,iy,iz);
    F3(f, 1, ix,iy,iz) = F3(flds, JYI, ix,iy,iz);
    F3(f, 2, ix,iy,iz) = F3(flds, JZI, ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_j_ec_ops = {
  .name      = "j_ec",
  .nr_comp   = 3,
  .fld_names = { "jx_ec", "jy_ec", "jz_ec" },
  .run       = calc_j_ec,
};

// ======================================================================

#define EX_NC(ix,iy,iz) (.5f * (F3(flds,  EX,ix,iy,iz) + F3(flds,  EX,ix-dx,iy,iz)))
#define EY_NC(ix,iy,iz) (.5f * (F3(flds,  EY,ix,iy,iz) + F3(flds,  EY,ix,iy-dy,iz)))
#define EZ_NC(ix,iy,iz) (.5f * (F3(flds,  EZ,ix,iy,iz) + F3(flds,  EZ,ix,iy,iz-dz)))

static void
calc_E_nc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	  struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = EX_NC(ix,iy,iz);
    F3(f, 1, ix,iy,iz) = EY_NC(ix,iy,iz);
    F3(f, 2, ix,iy,iz) = EZ_NC(ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e_nc_ops = {
  .name      = "e_nc",
  .nr_comp   = 3,
  .fld_names = { "ex_nc", "ey_nc", "ez_nc" },
  .run       = calc_E_nc,
};

// ======================================================================

#define EX_CC(ix,iy,iz) (.25f * (F3(flds,  EX,ix,iy,iz   ) + F3(flds,  EX,ix,iy+dy,iz   ) + \
				 F3(flds,  EX,ix,iy,iz+dz) + F3(flds,  EX,ix,iy+dy,iz+dz)))
#define EY_CC(ix,iy,iz) (.25f * (F3(flds,  EY,ix,iy,iz   ) + F3(flds,  EY,ix+dx,iy,iz   ) + \
				 F3(flds,  EY,ix,iy,iz+dz) + F3(flds,  EY,ix+dx,iy,iz+dz)))
#define EZ_CC(ix,iy,iz) (.25f * (F3(flds,  EZ,ix,iy   ,iz) + F3(flds,  EZ,ix+dx,iy   ,iz) + \
				 F3(flds,  EZ,ix,iy+dy,iz) + F3(flds,  EZ,ix+dx,iy+dy,iz)))

static void
calc_E_cc(struct psc_output_fields_item *item, struct psc_fields *flds_base, struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = EX_CC(ix,iy,iz);
    F3(f, 1, ix,iy,iz) = EY_CC(ix,iy,iz);
    F3(f, 2, ix,iy,iz) = EZ_CC(ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e_ops = {
  .name      = "e",
  .nr_comp   = 3,
  .fld_names = { "ex", "ey", "ez" },
  .run       = calc_E_cc,
};

// ======================================================================

static void
calc_E_ec(struct psc_output_fields_item *item, struct psc_fields *flds_base, struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = F3(flds, EX, ix,iy,iz);
    F3(f, 1, ix,iy,iz) = F3(flds, EY, ix,iy,iz);
    F3(f, 2, ix,iy,iz) = F3(flds, EZ, ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e_ec_ops = {
  .name      = "e_ec",
  .nr_comp   = 3,
  .fld_names = { "ex_ec", "ey_ec", "ez_ec" },
  .run       = calc_E_ec,
};

// ======================================================================

#define HX_NC(ix,iy,iz) (.25f*(F3(flds, HX,ix,iy,iz   ) + F3(flds, HX,ix,iy-dy,iz   ) + \
			       F3(flds, HX,ix,iy,iz-dz) + F3(flds, HX,ix,iy-dy,iz-dz)))
#define HY_NC(ix,iy,iz) (.25f*(F3(flds, HY,ix,iy,iz   ) + F3(flds, HY,ix-dx,iy,iz   ) + \
			       F3(flds, HY,ix,iy,iz-dz) + F3(flds, HY,ix-dx,iy,iz-dz)))
#define HZ_NC(ix,iy,iz) (.25f*(F3(flds, HZ,ix,iy   ,iz) + F3(flds, HZ,ix-dx,iy   ,iz) + \
			       F3(flds, HZ,ix,iy-dy,iz) + F3(flds, HZ,ix-dx,iy-dy,iz)))

static void
calc_H_nc(struct psc_output_fields_item *item, struct psc_fields *flds_base, struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, HX, HX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = HX_NC(ix,iy,iz);
    F3(f, 1, ix,iy,iz) = HY_NC(ix,iy,iz);
    F3(f, 2, ix,iy,iz) = HZ_NC(ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h_nc_ops = {
  .name      = "h_nc",
  .nr_comp   = 3,
  .fld_names = { "hx_nc", "hy_nc", "hz_nc" },
  .run       = calc_H_nc,
};

// ======================================================================

#define HX_CC(ix,iy,iz) (.5f*(F3(flds, HX,ix,iy,iz) + F3(flds, HX,ix+dx,iy,iz)))
#define HY_CC(ix,iy,iz) (.5f*(F3(flds, HY,ix,iy,iz) + F3(flds, HY,ix,iy+dy,iz)))
#define HZ_CC(ix,iy,iz) (.5f*(F3(flds, HZ,ix,iy,iz) + F3(flds, HZ,ix,iy,iz+dz)))

static void
calc_H_cc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	  struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, HX, HX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = HX_CC(ix,iy,iz);
    F3(f, 1, ix,iy,iz) = HY_CC(ix,iy,iz);
    F3(f, 2, ix,iy,iz) = HZ_CC(ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h_ops = {
  .name      = "h",
  .nr_comp   = 3,
  .fld_names = { "hx", "hy", "hz" },
  .run       = calc_H_cc,
};

// ======================================================================

static void
calc_H_fc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	  struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, HX, HX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = F3(flds, HX, ix,iy,iz);
    F3(f, 1, ix,iy,iz) = F3(flds, HY, ix,iy,iz);
    F3(f, 2, ix,iy,iz) = F3(flds, HZ, ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h_fc_ops = {
  .name      = "h_fc",
  .nr_comp   = 3,
  .fld_names = { "hx_fc", "hy_fc", "hz_fc" },
  .run       = calc_H_fc,
};

// ======================================================================

static void
calc_jdote_cc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	      struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, EX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = JX_CC(ix,iy,iz) * EX_CC(ix,iy,iz);
    F3(f, 1, ix,iy,iz) = JY_CC(ix,iy,iz) * EY_CC(ix,iy,iz);
    F3(f, 2, ix,iy,iz) = JZ_CC(ix,iy,iz) * EZ_CC(ix,iy,iz);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_jdote_ops = {
  .name      = "jdote",
  .nr_comp   = 3,
  .fld_names = { "jxex", "jyey", "jzez" },
  .run       = calc_jdote_cc,
};

// ======================================================================

static void
calc_poyn_cc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	     struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, HX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = (EY_CC(ix,iy,iz) * HZ_CC(ix,iy,iz) - 
			  EZ_CC(ix,iy,iz) * HY_CC(ix,iy,iz));
    F3(f, 1, ix,iy,iz) = (EZ_CC(ix,iy,iz) * HX_CC(ix,iy,iz) -
			  EX_CC(ix,iy,iz) * HZ_CC(ix,iy,iz));
    F3(f, 2, ix,iy,iz) = (EX_CC(ix,iy,iz) * HY_CC(ix,iy,iz) -
			  EY_CC(ix,iy,iz) * HX_CC(ix,iy,iz));
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_poyn_ops = {
  .name      = "poyn",
  .nr_comp   = 3,
  .fld_names = { "poynx", "poyny", "poynz" },
  .run       = calc_poyn_cc,
};

// ======================================================================

static void
calc_E2_cc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	   struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = sqr(EX_CC(ix,iy,iz));
    F3(f, 1, ix,iy,iz) = sqr(EY_CC(ix,iy,iz));
    F3(f, 2, ix,iy,iz) = sqr(EZ_CC(ix,iy,iz));
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e2_ops = {
  .name      = "e2",
  .nr_comp   = 3,
  .fld_names = { "ex2", "ey2", "ez2" },
  .run       = calc_E2_cc,
};

// ======================================================================

static void
calc_H2_cc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	   struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, HX, HX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = sqr(HX_CC(ix,iy,iz));
    F3(f, 1, ix,iy,iz) = sqr(HY_CC(ix,iy,iz));
    F3(f, 2, ix,iy,iz) = sqr(HZ_CC(ix,iy,iz));
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h2_ops = {
  .name      = "h2",
  .nr_comp   = 3,
  .fld_names = { "hx2", "hy2", "hz2" },
  .run       = calc_H2_cc,
};

// ======================================================================

static void
calc_divb(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	   struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, HX, HX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = 
      ((F3(flds, HX, ix+dx,iy,iz) - F3(flds, HX, ix,iy,iz)) / ppsc->patch[f->p].dx[0] +
       (F3(flds, HY, ix,iy+dy,iz) - F3(flds, HY, ix,iy,iz)) / ppsc->patch[f->p].dx[1] +
       (F3(flds, HZ, ix,iy,iz+dz) - F3(flds, HZ, ix,iy,iz)) / ppsc->patch[f->p].dx[2]);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_divb_ops = {
  .name      = "divb",
  .nr_comp   = 1,
  .fld_names = { "divb" },
  .run       = calc_divb,
};

// ======================================================================

static void
calc_divj_nc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	     struct psc_particles *prts, struct psc_fields *f)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, JXI + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, 0, ix,iy,iz) = 
      ((F3(flds, JXI, ix,iy,iz) - F3(flds, JXI, ix-dx,iy,iz)) / ppsc->patch[f->p].dx[0] +
       (F3(flds, JYI, ix,iy,iz) - F3(flds, JYI, ix,iy-dy,iz)) / ppsc->patch[f->p].dx[1] +
       (F3(flds, JZI, ix,iy,iz) - F3(flds, JZI, ix,iy,iz-dz)) / ppsc->patch[f->p].dx[2]);
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_divj_ops = {
  .name      = "divj",
  .nr_comp   = 1,
  .fld_names = { "divj" },
  .run       = calc_divj_nc,
};

