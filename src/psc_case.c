
#include "psc_case_private.h"

// ----------------------------------------------------------------------
// __psc_case_create

static void
__psc_case_create(struct _psc_case *_case)
{
  psc_create();
}

// ----------------------------------------------------------------------
// __psc_case_init_param

void
_psc_case_init_param(struct _psc_case *_case)
{
  psc_case_init_param(_psc_case->Case);
}

// ----------------------------------------------------------------------
// __psc_case_init_npt

void
_psc_case_init_npt(struct _psc_case *_case, int kind, double x[3],
		   struct psc_particle_npt *npt)
{
  _case->Case->ops->init_npt(_psc_case->Case, kind, x, npt);
}

// ----------------------------------------------------------------------
// __psc_case_init_field

void
_psc_case_init_field(struct _psc_case *_case, mfields_base_t *flds)
{
  _case->Case->ops->init_field(_case->Case, flds);
}

// ======================================================================
// _psc_case_init

static struct _psc_case_ops _psc_case_default_ops = {
  .name                  = "default",
};

static void
_psc_case_create_sub(struct _psc_case *_case)
{
  _case->Case = psc_case_create(_case->obj.ops->name);
}

static struct _psc_case_ops _psc_case_test_yz_ops = {
  .name                  = "test_yz",
  .create                = _psc_case_create_sub,
};

static struct _psc_case_ops _psc_case_test_xz_ops = {
  .name                  = "test_xz",
  .create                = _psc_case_create_sub,
};

static struct _psc_case_ops _psc_case_harris_xy_ops = {
  .name                  = "harris_xy",
  .create                = _psc_case_create_sub,
};

static void
_psc_case_init()
{
  mrc_class_register_subclass(&mrc_class__psc_case, &_psc_case_default_ops);
  mrc_class_register_subclass(&mrc_class__psc_case, &_psc_case_test_yz_ops);
  mrc_class_register_subclass(&mrc_class__psc_case, &_psc_case_test_xz_ops);
  mrc_class_register_subclass(&mrc_class__psc_case, &_psc_case_harris_xy_ops);
}

// ======================================================================
// _psc_case class

struct mrc_class__psc_case mrc_class__psc_case = {
  .name             = "_psc_case",
  .size             = sizeof(struct _psc_case),
  .init             = _psc_case_init,
  .create           = __psc_case_create,
};

