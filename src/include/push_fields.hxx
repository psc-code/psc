
#pragma once

#include "fields3d.hxx"

extern int pr_time_step_no_comm; // FIXME

// ======================================================================
// class PushFieldsBase

class PushFieldsBase
{
public:
  virtual void push_E(MfieldsBase& mflds_base, double dt_fac) = 0;
  virtual void push_H(MfieldsBase& mflds_base, double dt_fac) = 0;
};

