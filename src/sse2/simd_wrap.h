
#ifndef SSE2_CGEN_H
#define SSE2_CGEN_H

#include <xmmintrin.h>
#include <emmintrin.h>

#define VEC_OFF(foff,jy, jz) {				\
    foff.r = pv_sub_int(jz.r,ilg[2].r);					\
    itemp1.r = pv_sub_int(jy.r,ilg[1].r);				\
    foff.r = pv_mul_int(foff.r, img[1].r);				\
    foff.r = pv_add_int(foff.r, itemp1.r);				\
    foff.r = pv_mul_int(foff.r, img[0].r);				\
  }



#if SSE2_DOUBLE

#define VEC_SIZE 2
typedef double sse2_real;

union packed_vector{
  __m128d r;
  sse2_real v[2] __attribute__ ((aligned (32))); //FIXME : Might break for any non gcc
} ;

union packed_int{
  __m128i r;
  int v[2] __attribute__ ((aligned (16))); // this is a little strange, but only the lower
                                           // two values will ever be used.
} ; 

// real functions
#define pv_add_real(var1_r, var2_r) _mm_add_pd( var1_r, var2_r )
#define pv_sub_real(var1_r, var2_r) _mm_sub_pd( var1_r, var2_r )
#define pv_mul_real(var1_r, var2_r) _mm_mul_pd( var1_r, var2_r )
#define pv_div_real(var1_r, var2_r) _mm_div_pd( var1_r, var2_r )
#define pv_sqrt_real(var1_r) _mm_sqrt_pd( var1_r)
#define pv_set1_real(number) _mm_set1_pd( number )

// int functions
#define pv_add_int(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r )
#define pv_sub_int(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r )
#define pv_mul_int(var1_r, var2_r) _mm_mullo_epi16( var1_r, var2_r )
#define pv_set1_int(number) _mm_set1_epi32( number )

//conversion functions (round or pad)
#define pv_cvt_real_to_int(real_vec) _mm_cvtpd_epi32( real_vec )
#define pv_cvt_int_to_real(int_vec) _mm_cvtepi32_pd( int_vec )


#define INTERP_FIELD_YZ(F_ENUM, indx2, indx3, outer_coeff, inner_coeff, field) { \
    pvReal field##_in, field##tmp1, field##tmp2, field##tmp3;		\
    pvInt field##off1, field##off2, field##off3;			\
 									\
    VEC_OFF(field##off1,indx2##mns1,indx3##mns1);			\
    field##_in.v[0] = F_ENUM##point[field##off1.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off1.v[1]];	\
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##_in.r);	\
									\
    VEC_OFF(field##off2,indx2,indx3##mns1);				\
    field##_in.v[0] = F_ENUM##point[field##off2.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off2.v[1]];	\
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##_in.r);	\
									\
    VEC_OFF(field##off3,indx2##pls1,indx3##mns1);				\
    field##_in.v[0] = F_ENUM##point[field##off3.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off3.v[1]];	\
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_mul_real(outer_coeff##mz.r, field##tmp1.r);		\
    									\
    VEC_OFF(field##off1,indx2##mns1,indx3);				\
    field##_in.v[0] = F_ENUM##point[field##off1.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off1.v[1]];	\
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##_in.r);		\
    									\
    VEC_OFF(field##off2,indx2,indx3);				\
    field##_in.v[0] = F_ENUM##point[field##off2.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off2.v[1]];	\
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##_in.r);		\
    									\
    VEC_OFF(field##off3,indx2##pls1,indx3);				\
    field##_in.v[0] = F_ENUM##point[field##off3.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off3.v[1]];	\
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##Oz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
    									\
    VEC_OFF(field##off1,indx2##mns1,indx3##pls1);				\
    field##_in.v[0] = F_ENUM##point[field##off1.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off1.v[1]];	\
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##_in.r);		\
    									\
    VEC_OFF(field##off2,indx2,indx3##pls1);				\
    field##_in.v[0] = F_ENUM##point[field##off2.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off2.v[1]];	\
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##_in.r);		\
    									\
    VEC_OFF(field##off3,indx2##pls1,indx3##pls1);				\
    field##_in.v[0] = F_ENUM##point[field##off3.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off3.v[1]];	\
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##lz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
  }


#else
#define VEC_SIZE 4

typedef float sse2_real;

union packed_vector{
  __m128 r;
  sse2_real v[4] __attribute__ ((aligned (16))); //FIXME : Might break for any non gcc
} ;

union packed_int{
      __m128i r;
      int v[4] __attribute__ ((aligned (16)));
}; 

//THESE ARE ATOMICS!! Single variables only!!!!!
// ( Five exclamation marks, the sure sign of an insane mind. --Terry Pratchett)
// ( But I really mean it. --Steve)

// real functions
#define pv_add_real(var1_r, var2_r) _mm_add_ps( var1_r, var2_r )
#define pv_sub_real(var1_r, var2_r) _mm_sub_ps( var1_r, var2_r )
#define pv_mul_real(var1_r, var2_r) _mm_mul_ps( var1_r, var2_r )
#define pv_div_real(var1_r, var2_r) _mm_div_ps( var1_r, var2_r )
#define pv_sqrt_real(var1_r) _mm_sqrt_ps( var1_r)
#define pv_set1_real(number) _mm_set1_ps( number )
#define pv_load_real(float_ad) _mm_load_ps( float_ad )

// int functions
#define pv_add_int(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r )
#define pv_sub_int(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r )
#define pv_mul_int(var1_r, var2_r) _mm_mullo_epi16( var1_r, var2_r )
#define pv_set1_int(number) _mm_set1_epi32( number )

//conversion functions (round or pad)
#define pv_cvt_real_to_int(real_vec) _mm_cvtps_epi32( real_vec )
#define pv_cvt_int_to_real(int_vec) _mm_cvtepi32_ps( int_vec )


#define INTERP_FIELD_YZ(F_ENUM, indx2, indx3, outer_coeff, inner_coeff, field) { \
    pvReal field##_in, field##tmp1, field##tmp2, field##tmp3;		\
    pvInt field##off1, field##off2, field##off3;			\
    VEC_OFF(field##off1,indx2##mns1,indx3##mns1);			\
    field##_in.v[0] = F_ENUM##point[field##off1.v[0]];			\
    field##_in.v[1] = F_ENUM##point[field##off1.v[1]];			\
    field##_in.v[2] = F_ENUM##point[field##off1.v[2]];			\
    field##_in.v[3] = F_ENUM##point[field##off1.v[3]];			\
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##_in.r);	\
									\
    VEC_OFF(field##off2,indx2,indx3##mns1);				\
    field##_in.v[0] = F_ENUM##point[field##off2.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off2.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off2.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off2.v[3]];	\
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##_in.r);	\
									\
    VEC_OFF(field##off3,indx2##pls1,indx3##mns1);				\
    field##_in.v[0] = F_ENUM##point[field##off3.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off3.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off3.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off3.v[3]];	\
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_mul_real(outer_coeff##mz.r, field##tmp1.r);		\
    									\
    VEC_OFF(field##off1,indx2##mns1,indx3);				\
    field##_in.v[0] = F_ENUM##point[field##off1.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off1.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off1.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off1.v[3]];	\
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##_in.r);		\
    									\
    VEC_OFF(field##off2,indx2,indx3);				\
    field##_in.v[0] = F_ENUM##point[field##off2.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off2.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off2.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off2.v[3]];	\
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##_in.r);		\
    									\
    VEC_OFF(field##off3,indx2##pls1,indx3);				\
    field##_in.v[0] = F_ENUM##point[field##off3.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off3.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off3.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off3.v[3]];	\
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##Oz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
    									\
    VEC_OFF(field##off1,indx2##mns1,indx3##pls1);				\
    field##_in.v[0] = F_ENUM##point[field##off1.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off1.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off1.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off1.v[3]];	\
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##_in.r);		\
    									\
    VEC_OFF(field##off2,indx2,indx3##pls1);				\
    field##_in.v[0] = F_ENUM##point[field##off2.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off2.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off2.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off2.v[3]];	\
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##_in.r);		\
    									\
    VEC_OFF(field##off3,indx2##pls1,indx3##pls1);				\
    field##_in.v[0] = F_ENUM##point[field##off3.v[0]];	\
    field##_in.v[1] = F_ENUM##point[field##off3.v[1]];	\
    field##_in.v[2] = F_ENUM##point[field##off3.v[2]];	\
    field##_in.v[3] = F_ENUM##point[field##off3.v[3]];	\
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##lz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
  }

#endif

#endif
