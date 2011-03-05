
#include "mrc_ddc.h"
#include "mrc_ddc_private.h"

#include <mrc_params.h>

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define mrc_ddc_ops(domain) ((struct mrc_ddc_ops *) ddc->obj.ops)

static int
get_rank(struct mrc_ddc *ddc, const int proc[3])
{
  for (int d = 0; d < 3; d++) {
    assert(proc[d] >= 0 && proc[d] < ddc->prm.n_proc[d]);
  }
  return (proc[2] * ddc->prm.n_proc[1] + proc[1]) * ddc->prm.n_proc[0] + proc[0];
}

int
mrc_ddc_get_rank_nei(struct mrc_ddc *ddc, int dir[3])
{
  int proc_nei[3];
  for (int d = 0; d < 3; d++) {
    proc_nei[d] = ddc->proc[d] + dir[d];
    if (ddc->prm.bc[d] == BC_PERIODIC) {
      if (proc_nei[d] < 0) {
	proc_nei[d] += ddc->prm.n_proc[d];
      }
      if (proc_nei[d] >= ddc->prm.n_proc[d]) {
	proc_nei[d] -= ddc->prm.n_proc[d];
      }
    }
    if (proc_nei[d] < 0 || proc_nei[d] >= ddc->prm.n_proc[d])
      return - 1;
  }
  return get_rank(ddc, proc_nei);
}

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

// ======================================================================
// mrc_ddc_ops_f3 for mrc_f3

#include <mrc_fld.h>

// FIXME, 0-based offsets and ghost points don't match well (not pretty anyway)

static void
mrc_f3_copy_to_buf(int mb, int me, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct mrc_f3 *fld = ctx;
  float *buf = _buf;
  int bnd = fld->sw;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf,m - mb, ix,iy,iz) = MRC_F3(fld,m, ix+bnd,iy+bnd,iz+bnd);
	}
      }
    }
  }
}

static void
mrc_f3_copy_from_buf(int mb, int me, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct mrc_f3 *fld = ctx;
  float *buf = _buf;
  int bnd = fld->sw;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_F3(fld,m, ix+bnd,iy+bnd,iz+bnd) = MRC_DDC_BUF3(buf,m - mb, ix,iy,iz);
	}
      }
    }
  }
}

struct mrc_ddc_funcs mrc_ddc_funcs_f3 = {
  .copy_to_buf   = mrc_f3_copy_to_buf,
  .copy_from_buf = mrc_f3_copy_from_buf,
};

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
  { "n_proc"          , VAR(n_proc)       , PARAM_INT3(0, 0, 0)    },
  { "ilo"             , VAR(ilo)          , PARAM_INT3(0, 0, 0)    },
  { "ihi"             , VAR(ihi)          , PARAM_INT3(0, 0, 0)    },
  { "ibn"             , VAR(ibn)          , PARAM_INT3(0, 0, 0)    },
  { "bc"              , VAR(bc)           , PARAM_INT3(0, 0, 0)    }, // FIXME, select
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

