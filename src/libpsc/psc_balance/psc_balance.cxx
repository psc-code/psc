
#include "psc_balance_private.h"
#include "psc_bnd_fields.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "particles.hxx"
#include "bnd_particles.hxx"
#include "bnd.hxx"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <stdlib.h>
#include <string.h>

LIST_HEAD(psc_mfields_base_list);

double *psc_balance_comp_time_by_patch;

int psc_balance_generation_cnt;

// ----------------------------------------------------------------------
// _psc_balance_read

static void
_psc_balance_read(struct psc_balance *bal, struct mrc_io *io)
{
  int nr_patches;
  mrc_domain_get_patches(ppsc->mrc_domain, &nr_patches);
  psc_balance_comp_time_by_patch = new double[nr_patches];
}

// ----------------------------------------------------------------------
// psc_balance_init

extern struct psc_balance_ops psc_balance_double_ops;
extern struct psc_balance_ops psc_balance_single_ops;

static void
psc_balance_init(void)
{
  mrc_class_register_subclass(&mrc_class_psc_balance, &psc_balance_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_balance, &psc_balance_single_ops);
}

// ======================================================================
// psc_balance class

#define VAR(x) (void *)offsetof(struct psc_balance, x)
static struct param psc_balance_descr[] = {
  { "every"            , VAR(every)            , PARAM_INT(0)            },
  { "factor_fields"    , VAR(factor_fields)    , PARAM_DOUBLE(1.)        },
  { "print_loads"      , VAR(print_loads)      , PARAM_BOOL(false)       },
  { "write_loads"      , VAR(write_loads)      , PARAM_BOOL(false)       },
  {},
};
#undef VAR

struct mrc_class_psc_balance_ : mrc_class_psc_balance
{
  mrc_class_psc_balance_()
  {
    name        = "psc_balance";
    size        = sizeof(struct psc_balance);
    param_descr = psc_balance_descr;
    init        = psc_balance_init;
    read        = _psc_balance_read;
  }
} mrc_class_psc_balance;
