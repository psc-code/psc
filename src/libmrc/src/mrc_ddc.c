
#include "mrc_ddc.h"
#include "mrc_ddc_private.h"

#include <mrc_params.h>

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define mrc_ddc_ops(domain) ((struct mrc_ddc_ops *) ddc->obj.ops)

// ----------------------------------------------------------------------
// mrc_ddc_create

static void
_mrc_ddc_create(struct mrc_obj *obj)
{
  struct mrc_ddc *ddc = to_mrc_ddc(obj);

  MPI_Comm_rank(obj->comm, &ddc->rank);
  MPI_Comm_size(obj->comm, &ddc->size);
}

// ----------------------------------------------------------------------
// mrc_ddc_set_funcs

void
mrc_ddc_set_funcs(struct mrc_ddc *ddc, struct mrc_ddc_funcs *funcs)
{
  ddc->funcs = funcs;
}

// ----------------------------------------------------------------------
// mrc_ddc_fill_ghosts

void
mrc_ddc_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_ops *ops = mrc_ddc_ops(ddc);
  assert(ops->fill_ghosts);
  ops->fill_ghosts(ddc, mb, me, ctx);
}

// ----------------------------------------------------------------------
// mrc_ddc_add_ghosts

void
mrc_ddc_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  struct mrc_ddc_ops *ops = mrc_ddc_ops(ddc);
  assert(ops->add_ghosts);
  ops->add_ghosts(ddc, mb, me, ctx);
}

// ----------------------------------------------------------------------
// mrc_ddc_get_rank_nei

int
mrc_ddc_get_rank_nei(struct mrc_ddc *ddc, int dir[3])
{
  struct mrc_ddc_ops *ops = mrc_ddc_ops(ddc);
  assert(ops->get_rank_nei);
  return ops->get_rank_nei(ddc, dir);
}

// ======================================================================
// mrc_ddc class

static LIST_HEAD(mrc_ddc_subclasses);

void
libmrc_ddc_register(struct mrc_ddc_ops *ops)
{
  list_add_tail(&ops->list, &mrc_ddc_subclasses);
}

// ----------------------------------------------------------------------
// mrc_ddc_init

static void
mrc_ddc_init()
{
  libmrc_ddc_register_simple();
}

#define VAR(x) (void *)offsetof(struct mrc_ddc_params, x)
static struct param mrc_ddc_params_descr[] = {
  { "size_of_type"    , VAR(size_of_type) , PARAM_INT(0)           },
  { "max_n_fields"    , VAR(max_n_fields) , PARAM_INT(1)           },
  { "ibn"             , VAR(ibn)          , PARAM_INT3(0, 0, 0)    },
  {},
};
#undef VAR

struct mrc_class mrc_class_mrc_ddc = {
  .name             = "mrc_ddc",
  .size             = sizeof(struct mrc_ddc),
  .param_descr      = mrc_ddc_params_descr,
  .param_offset     = offsetof(struct mrc_ddc, prm),
  .subclasses       = &mrc_ddc_subclasses,
  .init             = mrc_ddc_init,
  .create           = _mrc_ddc_create,
};

