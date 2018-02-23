
#pragma once

#include "psc_push_fields_private.h"
#include "fields3d.hxx"

#include "psc_bnd_fields.h"
#include "psc_bnd.h"

#include <mrc_profile.h>

extern int pr_time_step_no_comm; // FIXME

// ======================================================================
// PscPushFieldsBase

template<typename S>
struct PscPushFieldsBase
{
  using sub_t = S;

  explicit PscPushFieldsBase(psc_push_fields *pushf)
    : pushf_(pushf),
      sub_(mrc_to_subobj(pushf, sub_t))
  {}

  void advance_E(mfields_base_t mflds, double frac)
  {
    struct psc_push_fields_ops *ops = psc_push_fields_ops(pushf_);
    static int pr;
    if (!pr) {
      pr = prof_register("push_fields_E", 1., 0, 0);
    }
    
    psc_stats_start(st_time_field);
    prof_start(pr);
    prof_restart(pr_time_step_no_comm);
    
    assert(ops->push_mflds_E);
    ops->push_mflds_E(pushf_, mflds.mflds(), frac);
    
    prof_stop(pr_time_step_no_comm);
    prof_stop(pr);
    psc_stats_stop(st_time_field);
  }

  void advance_H(mfields_base_t mflds, double frac)
  {
    struct psc_push_fields_ops *ops = psc_push_fields_ops(pushf_);
    static int pr;
    if (!pr) {
      pr = prof_register("push_fields_H", 1., 0, 0);
    }
    
    psc_stats_start(st_time_field);
    prof_start(pr);
    prof_restart(pr_time_step_no_comm);
    
    assert(ops->push_mflds_H);
    ops->push_mflds_H(pushf_, mflds.mflds(), frac);
    
    prof_stop(pr);
    prof_stop(pr_time_step_no_comm);
    psc_stats_stop(st_time_field);
  }

  void advance_b2(mfields_base_t mflds)
  {
    // fill ghosts for H
    psc_bnd_fields_fill_ghosts_H(pushf_->bnd_fields, mflds.mflds());
    psc_bnd_fill_ghosts(ppsc->bnd, mflds.mflds(), HX, HX + 3);
    
    // add and fill ghost for J
    psc_bnd_fields_add_ghosts_J(pushf_->bnd_fields, mflds.mflds());
    psc_bnd_add_ghosts(ppsc->bnd, mflds.mflds(), JXI, JXI + 3);
    psc_bnd_fill_ghosts(ppsc->bnd, mflds.mflds(), JXI, JXI + 3);
    
    // push E
    advance_E(mflds, 1.);

    psc_bnd_fields_fill_ghosts_E(pushf_->bnd_fields, mflds.mflds());
    if (pushf_->variant == 0) {
      psc_bnd_fill_ghosts(ppsc->bnd, mflds.mflds(), EX, EX + 3);
    }
  }

  void advance_a(mfields_base_t mflds)
  {
    if (pushf_->variant == 0) {
      psc_bnd_fields_fill_ghosts_E(pushf_->bnd_fields, mflds.mflds());
      psc_bnd_fill_ghosts(ppsc->bnd, mflds.mflds(), EX, EX + 3);
    }
    
    // push H
    advance_H(mflds, .5);

    psc_bnd_fields_fill_ghosts_H(pushf_->bnd_fields, mflds.mflds());
    if (pushf_->variant == 0) {
      psc_bnd_fill_ghosts(ppsc->bnd, mflds.mflds(), HX, HX + 3);
    }
  }

  psc_push_fields *pushf() { return pushf_; }
  
  sub_t* operator->() { return sub_; }

  sub_t* sub() { return sub_; }

private:
  psc_push_fields *pushf_;
  sub_t *sub_;
};

