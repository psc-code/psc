
#include "psc_sse2.h"
#include "psc_push_particles_private.h"

#include <stdlib.h>
#include <assert.h>

//---------------------------------------------
/// Fills the global constant numerical vectors used by sse2. 

/// **Must** be called before they are used
void 
init_vec_numbers(void) {		
  ones.r = pv_set1_real(1.0);			
  half.r = pv_set1_real(.5);			
  onepfive.r = pv_set1_real(1.5);		
  threefourths.r = pv_set1_real(.75);		
  third.r = pv_set1_real(1./3.);		
  ione.r = pv_set1_int(1);			
}

/// Pointers to functions optimized for SSE2

// ======================================================================
// psc_push_particles: subclass "sse2"

struct psc_push_particles_ops psc_push_particles_sse2_ops = {
  .name                  = "sse2",
  .push_yz               = psc_push_particles_sse2_push_yz,
  .push_yz_a             = psc_push_particles_sse2_push_yz_a,
  .push_yz_b             = psc_push_particles_sse2_push_yz_b,
};

/// \file psc_sse2.c Backend functions for SSE2 implementation.
