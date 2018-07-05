
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
  // dtor

  ~Psc()
  {
    psc_destroy(psc_);
  }
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    psc_view(psc_);
    psc_mparticles_view(psc_->particles);
    psc_mfields_view(psc_->flds);
  
    psc_method_initialize(psc_->method, psc_);
  }

protected:
  psc* psc_;
};
