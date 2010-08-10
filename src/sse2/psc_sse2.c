
#include "psc_sse2.h"

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

struct psc_ops psc_ops_sse2 = {
  .name = "sse2",
  .push_part_yz_a         = sse2_push_part_yz_a,
  .push_part_yz_b         = sse2_push_part_yz_b,
  .push_part_yz           = sse2_push_part_yz,
  .push_part_xz           = sse2_push_part_xz,
}; 

/// \file psc_sse2.c Backend functions for SSE2 implementation.
