
#include "psc_push_fields_private.h"

#include "psc_bnd.h"
#include "psc_bnd_fields.h"
#include <mrc_profile.h>
#include <mrc_params.h>

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

struct psc_bnd_fields *
psc_push_fields_get_bnd_fields(struct psc_push_fields *push)
{
  return push->bnd_fields;
}

// ======================================================================
// forward to subclass

void
psc_push_fields_push_E(struct psc_push_fields *push, struct psc_mfields *flds,
		       double dt_fac)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_E", 1., 0, 0);
  }
  
  psc_stats_start(st_time_field);
  prof_start(pr);
  prof_restart(pr_time_step_no_comm);

  assert(ops->push_mflds_E);
  ops->push_mflds_E(push, flds, dt_fac);

  prof_stop(pr_time_step_no_comm);
  prof_stop(pr);
  psc_stats_stop(st_time_field);
}

void
psc_push_fields_push_H(struct psc_push_fields *push, struct psc_mfields *flds,
		       double dt_fac)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_H", 1., 0, 0);
  }
  
  psc_stats_start(st_time_field);
  prof_start(pr);
  prof_restart(pr_time_step_no_comm);

  assert(ops->push_mflds_H);
  ops->push_mflds_H(push, flds, dt_fac);

  prof_stop(pr);
  prof_stop(pr_time_step_no_comm);
  psc_stats_stop(st_time_field);
}

// ======================================================================
// variant 0: default
// 
// That's the traditional way of how things have been done
//
// ======================================================================
// variant 1: opt
//
// This way does only the minimum amount of communication needed,
// and does all the communication (including particles) at the same time

// ======================================================================

// ======================================================================
// psc_push_fields_init

extern struct psc_push_fields_ops psc_push_fields_c_ops;
extern struct psc_push_fields_ops psc_push_fields_single_ops;
extern struct psc_push_fields_ops psc_push_fields_fortran_ops;
extern struct psc_push_fields_ops psc_push_fields_cbe_ops;
extern struct psc_push_fields_ops psc_push_fields_cuda_ops;
extern struct psc_push_fields_ops psc_push_fields_cuda2_ops;
extern struct psc_push_fields_ops psc_push_fields_acc_ops;
extern struct psc_push_fields_ops psc_push_fields_vpic_ops;
extern struct psc_push_fields_ops psc_push_fields_none_ops;

static void
psc_push_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_none_ops);
#ifdef USE_FORTRAN
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_fortran_ops);
#endif
#ifdef USE_CBE
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_cbe_ops);
#endif
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_cuda_ops);
#endif
#ifdef USE_CUDA2
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_cuda2_ops);
#endif
#ifdef USE_ACC
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_acc_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_vpic_ops);
#endif
}

// ======================================================================

#define VAR(x) (void *)offsetof(struct psc_push_fields, x)
static struct param psc_push_fields_descr[] = {
  { "variant"          , VAR(variant)          , PARAM_INT(0),
  .help = "selects different variants of the EM solver: "
  "(0): default (traditional way with 4 fill_ghosts()) "
  "(1): optimized with only 1 fill_ghosts(), 1st order "
  "particle shape only" },

  { "bnd_fields"       , VAR(bnd_fields)       , MRC_VAR_OBJ(psc_bnd_fields), },

  {},
};
#undef VAR

// ======================================================================
// psc_push_fields class

struct mrc_class_psc_push_fields_ : mrc_class_psc_push_fields {
  mrc_class_psc_push_fields_() {
    name             = "psc_push_fields";
    size             = sizeof(struct psc_push_fields);
    param_descr      = psc_push_fields_descr;
    init             = psc_push_fields_init;
  }
} mrc_class_psc_push_fields;

