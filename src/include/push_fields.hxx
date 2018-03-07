
#pragma once

#include "psc_push_fields_private.h"
#include "fields3d.hxx"
#include "bnd.hxx"
#include "bnd_fields.hxx"

#include "psc_bnd_fields.h"

#include <mrc_profile.h>

extern int pr_time_step_no_comm; // FIXME

// ======================================================================
// PscPushFields

template<typename S>
struct PscPushFields
{
  using sub_t = S;

  explicit PscPushFields(psc_push_fields *pushf)
    : pushf_(pushf),
      sub_(mrc_to_subobj(pushf, sub_t))
  {}

  void advance_E(PscMfieldsBase mflds, double frac)
  {
    struct psc_push_fields_ops *ops = psc_push_fields_ops(pushf_);
    static int pr;
    if (!pr) {
      pr = prof_register("push_fields_E", 1., 0, 0);
    }
    
    psc_stats_start(st_time_field);
    prof_start(pr);
    prof_restart(pr_time_step_no_comm);
    
    sub()->push_E(mflds.mflds(), frac);
    
    prof_stop(pr_time_step_no_comm);
    prof_stop(pr);
    psc_stats_stop(st_time_field);
  }

  void advance_H(PscMfieldsBase mflds, double frac)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("push_fields_H", 1., 0, 0);
    }
    
    psc_stats_start(st_time_field);
    prof_start(pr);
    prof_restart(pr_time_step_no_comm);
    
    sub()->push_H(mflds.mflds(), frac);
    
    prof_stop(pr);
    prof_stop(pr_time_step_no_comm);
    psc_stats_stop(st_time_field);
  }

  void advance_b2(PscMfieldsBase mflds)
  {
    auto bnd = PscBndBase(ppsc->bnd);
    auto bndf = PscBndFieldsBase(pushf_->bnd_fields);

    // fill ghosts for H
    bndf.fill_ghosts_H(mflds);
    bnd.fill_ghosts(mflds, HX, HX + 3);
    
    // add and fill ghost for J
    bndf.add_ghosts_J(mflds);
    bnd.add_ghosts(mflds, JXI, JXI + 3);
    bnd.fill_ghosts(mflds, JXI, JXI + 3);
    
    // push E
    advance_E(mflds, 1.);

    bndf.fill_ghosts_E(mflds);
    if (pushf_->variant == 0) {
      bnd.fill_ghosts(mflds, EX, EX + 3);
    }
  }

  void advance_a(PscMfieldsBase mflds)
  {
    auto bnd = PscBndBase(ppsc->bnd);
    auto bndf = PscBndFieldsBase(pushf_->bnd_fields);

    if (pushf_->variant == 0) {
      bndf.fill_ghosts_E(mflds);
      bnd.fill_ghosts(mflds, EX, EX + 3);
    }
    
    // push H
    advance_H(mflds, .5);

    bndf.fill_ghosts_H(mflds);
    if (pushf_->variant == 0) {
      bnd.fill_ghosts(mflds, HX, HX + 3);
    }
  }

  psc_push_fields *pushf() { return pushf_; }
  
  sub_t* operator->() { return sub_; }

  sub_t* sub() { return sub_; }

private:
  psc_push_fields *pushf_;
  sub_t *sub_;
};

// ======================================================================
// class PushFieldsBase

class PushFieldsBase
{
public:
  virtual void push_E(struct psc_mfields *mflds_base, double dt_fac) = 0;
  virtual void push_H(struct psc_mfields *mflds_base, double dt_fac) = 0;
};

using PscPushFieldsBase = PscPushFields<PushFieldsBase>;
