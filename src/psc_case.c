
#include "psc_case_private.h"

// ----------------------------------------------------------------------
// __psc_case_create

static void
__psc_case_create(struct _psc_case *_case)
{
  psc_create();
}

// ----------------------------------------------------------------------
// __psc_case_set_from_options

static void
__psc_case_set_from_options(struct _psc_case *_case)
{
  psc_set_from_options();
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
  if (_case->Case->ops->create) {
    _case->Case->ops->create(_case->Case);
  }
}

static void
_psc_case_sub_set_from_options(struct _psc_case *_case)
{
  struct param *ctx_descr = _case->Case->ops->ctx_descr;
  if (ctx_descr) {
    const char *case_name = _case->obj.ops->name;
    char cn[strlen(case_name) + 6];
    sprintf(cn, "case %s", case_name);
    mrc_params_parse(_case->Case->ctx, ctx_descr, cn, MPI_COMM_WORLD);
    mrc_params_print(_case->Case->ctx, ctx_descr, cn, MPI_COMM_WORLD);
  }

  // update psc-params accordingly
  if (_case->Case->ops->init_param) {
    _case->Case->ops->init_param(_case->Case);
  }
}

static struct _psc_case_ops _psc_case_test_yz_ops = {
  .name                  = "test_yz",
  .create                = _psc_case_create_sub,
  .set_from_options      = _psc_case_sub_set_from_options,
};

static struct _psc_case_ops _psc_case_test_xz_ops = {
  .name                  = "test_xz",
  .create                = _psc_case_create_sub,
  .set_from_options      = _psc_case_sub_set_from_options,
};

static struct _psc_case_ops _psc_case_harris_xy_ops = {
  .name                  = "harris_xy",
  .create                = _psc_case_create_sub,
  .set_from_options      = _psc_case_sub_set_from_options,
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
  .set_from_options = __psc_case_set_from_options,
};

