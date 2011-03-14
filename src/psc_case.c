
#include "psc_case_private.h"

// ----------------------------------------------------------------------
// __psc_case_set_from_options

static void
__psc_case_set_from_options(struct _psc_case *_case)
{
  // FIXME, should go away when we use type
  {
    struct mrc_obj *obj = &_case->obj;
    struct mrc_class *class = (struct mrc_class *) &mrc_class__psc_case;
    char *p = (char *) obj + class->param_offset;
    mrc_params_parse_nodefault(p, class->param_descr, mrc_obj_name(obj), obj->comm);
    mrc_params_parse_pfx(p, class->param_descr, mrc_obj_name(obj), obj->comm);
  }

  assert(_case->case_name);
  _case->Case = psc_case_create(_case->case_name);
}

// ======================================================================
// _psc_case_init

static void
_psc_case_init()
{
}

// ======================================================================
// _psc_case class

#define VAR(x) (void *)offsetof(struct _psc_case, x)
static struct param _psc_case_descr[] = {
  { "case"          , VAR(case_name)        , PARAM_STRING(NULL)   },
  {},
};
#undef VAR

struct mrc_class__psc_case mrc_class__psc_case = {
  .name             = "_psc_case",
  .size             = sizeof(struct _psc_case),
  .param_descr      = _psc_case_descr,
  .init             = _psc_case_init,
  .set_from_options = __psc_case_set_from_options,
};

