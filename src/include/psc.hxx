
#pragma once

#include <psc_method.h>

struct Psc
{
  // ----------------------------------------------------------------------
  // ctor

  Psc(psc* psc)
    : psc_(psc)
  {}
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize(psc* psc_)
  {
    psc_method_initialize(psc_->method, psc_);
  }

protected:
  psc* psc_;
};
