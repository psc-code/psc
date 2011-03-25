
#include "psc_bnd_fields_private.h"
#include "psc_pulse.h"

// ----------------------------------------------------------------------
// psc_bnd_fields_create

static void
_psc_bnd_fields_create(struct psc_bnd_fields *bnd)
{
  bnd->pulse_x1 = psc_pulse_create(psc_bnd_fields_comm(bnd));
  psc_pulse_set_name(bnd->pulse_x1, "pulse_x1");
  bnd->pulse_x2 = psc_pulse_create(psc_bnd_fields_comm(bnd));
  psc_pulse_set_name(bnd->pulse_x2, "pulse_x2");
  bnd->pulse_y1 = psc_pulse_create(psc_bnd_fields_comm(bnd));
  psc_pulse_set_name(bnd->pulse_y1, "pulse_y1");
  bnd->pulse_y2 = psc_pulse_create(psc_bnd_fields_comm(bnd));
  psc_pulse_set_name(bnd->pulse_y2, "pulse_y2");
  bnd->pulse_z1 = psc_pulse_create(psc_bnd_fields_comm(bnd));
  psc_pulse_set_name(bnd->pulse_z1, "pulse_z1");
  bnd->pulse_z2 = psc_pulse_create(psc_bnd_fields_comm(bnd));
  psc_pulse_set_name(bnd->pulse_z2, "pulse_z2");
}

// ----------------------------------------------------------------------
// psc_bnd_fields_set_from_options

static void
_psc_bnd_fields_set_from_options(struct psc_bnd_fields *bnd)
{
  psc_pulse_set_from_options(bnd->pulse_x1);
  psc_pulse_set_from_options(bnd->pulse_x2);
  psc_pulse_set_from_options(bnd->pulse_y1);
  psc_pulse_set_from_options(bnd->pulse_y2);
  psc_pulse_set_from_options(bnd->pulse_z1);
  psc_pulse_set_from_options(bnd->pulse_z2);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_setup

static void
_psc_bnd_fields_setup(struct psc_bnd_fields *bnd)
{
  psc_pulse_setup(bnd->pulse_x1);
  psc_pulse_setup(bnd->pulse_x2);
  psc_pulse_setup(bnd->pulse_y1);
  psc_pulse_setup(bnd->pulse_y2);
  psc_pulse_setup(bnd->pulse_z1);
  psc_pulse_setup(bnd->pulse_z2);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_view

static void
_psc_bnd_fields_view(struct psc_bnd_fields *bnd)
{
  psc_pulse_view(bnd->pulse_x1);
  psc_pulse_view(bnd->pulse_x2);
  psc_pulse_view(bnd->pulse_y1);
  psc_pulse_view(bnd->pulse_y2);
  psc_pulse_view(bnd->pulse_z1);
  psc_pulse_view(bnd->pulse_z2);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_destroy

static void
_psc_bnd_fields_destroy(struct psc_bnd_fields *bnd)
{
  psc_pulse_destroy(bnd->pulse_x1);
  psc_pulse_destroy(bnd->pulse_x2);
  psc_pulse_destroy(bnd->pulse_y1);
  psc_pulse_destroy(bnd->pulse_y2);
  psc_pulse_destroy(bnd->pulse_z1);
  psc_pulse_destroy(bnd->pulse_z2);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_get_pulse_*

struct psc_pulse *
psc_bnd_fields_get_pulse_x1(struct psc_bnd_fields *bnd)
{
  return bnd->pulse_x1;
}

struct psc_pulse *
psc_bnd_fields_get_pulse_x2(struct psc_bnd_fields *bnd)
{
  return bnd->pulse_x2;
}

struct psc_pulse *
psc_bnd_fields_get_pulse_y1(struct psc_bnd_fields *bnd)
{
  return bnd->pulse_y1;
}

struct psc_pulse *
psc_bnd_fields_get_pulse_y2(struct psc_bnd_fields *bnd)
{
  return bnd->pulse_y2;
}

struct psc_pulse *
psc_bnd_fields_get_pulse_z1(struct psc_bnd_fields *bnd)
{
  return bnd->pulse_z1;
}

struct psc_pulse *
psc_bnd_fields_get_pulse_z2(struct psc_bnd_fields *bnd)
{
  return bnd->pulse_z2;
}

// ======================================================================
// forward to subclass

void
psc_bnd_fields_fill_ghosts_b_H(struct psc_bnd_fields *bnd, mfields_base_t *flds)
{
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bnd);
  assert(ops->fill_ghosts_b_H);
  ops->fill_ghosts_b_H(bnd, flds);
}

void
psc_bnd_fields_setup_fields(struct psc_bnd_fields *bnd_fields, mfields_base_t *mflds)
{
  struct psc *_psc = &psc;
  struct psc_pulse *pulse;

  psc_foreach_patch(&psc, p) {
    struct psc *psc = _psc; // FIXME, ugly hack, don't ask...
    fields_base_t *pf = &mflds->f[p];
    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dx = psc->dx[0], dy = psc->dx[1], dz = psc->dx[2];
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

      pulse = psc_bnd_fields_get_pulse_x1(bnd_fields);
      F3_BASE(pf, EZ, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx        , yy        , zz + .5*dz, 0.);
      F3_BASE(pf, HY, jx,jy,jz) += -psc_pulse_field_p(pulse, xx + .5*dx, yy        , zz + .5*dz, 0.);
      F3_BASE(pf, EY, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx        , yy + .5*dy, zz        , 0.);
      F3_BASE(pf, HZ, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx + .5*dx, yy + .5*dy, zz        , 0.);
      
      pulse = psc_bnd_fields_get_pulse_x2(bnd_fields);
      F3_BASE(pf, EZ, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx        , yy        , zz + .5*dz, 0.);
      F3_BASE(pf, HY, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx + .5*dx, yy        , zz + .5*dz, 0.);
      F3_BASE(pf, EY, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx        , yy + .5*dy, zz        , 0.);
      F3_BASE(pf, HZ, jx,jy,jz) += -psc_pulse_field_s(pulse, xx + .5*dx, yy + .5*dy, zz        , 0.);
      
      pulse = psc_bnd_fields_get_pulse_y1(bnd_fields);
      F3_BASE(pf, EX, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx + .5*dx, yy        , zz        , 0.);
      F3_BASE(pf, HZ, jx,jy,jz) += -psc_pulse_field_p(pulse, xx + .5*dx, yy + .5*dy, zz        , 0.);
      F3_BASE(pf, EZ, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx        , yy        , zz + .5*dz, 0.);
      F3_BASE(pf, HX, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx        , yy + .5*dy, zz + .5*dz, 0.);
      
      pulse = psc_bnd_fields_get_pulse_y2(bnd_fields);
      F3_BASE(pf, EX, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx + .5*dx, yy        , zz        , 0.);
      F3_BASE(pf, HZ, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx + .5*dx, yy + .5*dy, zz        , 0.);
      F3_BASE(pf, EZ, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx        , yy        , zz + .5*dz, 0.);
      F3_BASE(pf, HX, jx,jy,jz) += -psc_pulse_field_s(pulse, xx        , yy + .5*dy, zz + .5*dz, 0.);
      
      pulse = psc_bnd_fields_get_pulse_z1(bnd_fields);
      F3_BASE(pf, EY, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx        , yy + .5*dy, zz        , 0.);
      F3_BASE(pf, HX, jx,jy,jz) += -psc_pulse_field_p(pulse, xx        , yy + .5*dy, zz + .5*dz, 0.);
      F3_BASE(pf, EX, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx + .5*dx, yy        , zz        , 0.);
      F3_BASE(pf, HY, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx + .5*dx, yy        , zz + .5*dz, 0.);
      
      pulse = psc_bnd_fields_get_pulse_z2(bnd_fields);
      F3_BASE(pf, EY, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx        , yy + .5*dy, zz        , 0.);
      F3_BASE(pf, HX, jx,jy,jz) +=  psc_pulse_field_p(pulse, xx        , yy + .5*dy, zz + .5*dz, 0.);
      F3_BASE(pf, EX, jx,jy,jz) +=  psc_pulse_field_s(pulse, xx + .5*dx, yy        , zz        , 0.);
      F3_BASE(pf, HY, jx,jy,jz) += -psc_pulse_field_s(pulse, xx + .5*dx, yy        , zz + .5*dz, 0.);
    } psc_foreach_3d_g_end;
  }
}

// ======================================================================
// psc_bnd_fields_init

static void
psc_bnd_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_fortran_ops);
}

// ======================================================================
// psc_bnd_fields class

struct mrc_class_psc_bnd_fields mrc_class_psc_bnd_fields = {
  .name             = "psc_bnd_fields",
  .size             = sizeof(struct psc_bnd_fields),
  .init             = psc_bnd_fields_init,
  .create           = _psc_bnd_fields_create,
  .destroy          = _psc_bnd_fields_destroy,
  .set_from_options = _psc_bnd_fields_set_from_options,
  .setup            = _psc_bnd_fields_setup,
  .view             = _psc_bnd_fields_view,
};

