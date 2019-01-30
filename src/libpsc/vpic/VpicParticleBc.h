
#ifndef VPIC_PARTICLE_BC_H
#define VPIC_PARTICLE_BC_H

#include "VpicListBase.h"

#include <mrc_common.h>
#include <cassert>

#define IN_boundary
#include "boundary/boundary_private.h"

// ======================================================================
// VpicParticleBc

struct VpicParticleBc : particle_bc_t
{
};

// ======================================================================
// VpicParticleBcList

struct VpicParticleBcList : public VpicListBase<VpicParticleBc>
{
  typedef VpicParticleBc ParticleBc;
  
  ParticleBc* append(ParticleBc* pbc)
  {
    return static_cast<ParticleBc*>(::append_particle_bc(pbc, reinterpret_cast<particle_bc_t**>(&head_)));
  }

  operator const particle_bc_t* () const
  {
    return head_;
  }
};


#endif
