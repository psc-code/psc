#ifndef SIMD_PUSH_COMMON_H
#define SIMD_PUSH_COMMON_H
#include "psc_sse2.h"
#include "simd_wrap.h"

// I'm not sure about the exact rules in the standard, but these are actually declared
// inline where they are defined. This *should* mean that they are treated as inline functions
// when appropriate. There's a lot of hazyness about this, and it may cause problems
// further down the line.

void calc_vi(struct particle_vec *p, pvReal * restrict vxi, pvReal * restrict vyi,
		    pvReal * restrict vzi, pvReal * restrict root);
void find_index_Cround(pvReal *xi, pvReal *dxi, pvInt * restrict j, pvReal * restrict h);
void find_index_Iround(pvReal *xi, pvReal *dxi, pvInt * restrict j, pvReal * restrict h);
void find_index_minus_shift(pvReal *xi, pvReal *dxi, pvInt * restrict j, pvReal * restrict h, pvReal *shift);
void ip_to_grid_m(pvReal *h, pvReal * restrict xmx);
void ip_to_grid_O(pvReal *h, pvReal * restrict xOx);
void ip_to_grid_l(pvReal *h, pvReal * restrict xlx);
void form_factor_m(pvReal *h, pvReal * restrict xmx);
void form_factor_O(pvReal *h, pvReal * restrict xOx);
void form_factor_l(pvReal *h, pvReal * restrict xlx);
void push_pi_dt(struct particle_vec *p, pvReal *exq, pvReal *eyq, pvReal *ezq, pvReal *bxq, pvReal *byq, pvReal *bzq, pvReal *dqs);
#endif
