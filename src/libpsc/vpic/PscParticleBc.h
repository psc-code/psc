
#ifndef PSC_PARTICLE_BC_H
#define PSC_PARTICLE_BC_H

#include "VpicListBase.h"

#include <mrc_common.h>
#include <cassert>

#define IN_boundary
#include "boundary/boundary_private.h"

// ======================================================================
// PscParticleBc

struct PscParticleBc : particle_bc_t
{
};

// ======================================================================
// PscParticleBcList

struct PscParticleBcList : public VpicListBase<PscParticleBc>
{
  typedef PscParticleBc ParticleBc;
  
  ParticleBc* append(ParticleBc* pbc)
  {
    assert(0);
  }

  operator const particle_bc_t* () const
  {
    return head_;
  }
};


#endif
