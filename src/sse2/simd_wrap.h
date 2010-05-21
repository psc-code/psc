
#ifndef SIMD_WRAP_H
#define SIMD_WRAP_H

#include <xmmintrin.h>
#include <emmintrin.h>

//---------------------------------------------
// Calculates the offset for accessing pointers to 
// flattened (2-D only) fields using vector operations

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
  int v[2] __attribute__ ((aligned (16))); 
                                           
} ; 

// real functions
#define pv_add_real(var1_r, var2_r) _mm_add_pd( var1_r, var2_r )
#define pv_sub_real(var1_r, var2_r) _mm_sub_pd( var1_r, var2_r )
#define pv_mul_real(var1_r, var2_r) _mm_mul_pd( var1_r, var2_r )
#define pv_div_real(var1_r, var2_r) _mm_div_pd( var1_r, var2_r )
#define pv_sqrt_real(var1_r) _mm_sqrt_pd( var1_r)
#define pv_set1_real(number) _mm_set1_pd( number )
#define pv_load1_real(pointer) _mm_load_sd( pointer )
#define pv_unpacklo_real(var1_r, var2_r) _mm_unpacklo_pd( var1_r, var2_r )

// int functions
#define pv_add_int(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r )
#define pv_sub_int(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r )
#define pv_mul_int(var1_r, var2_r) _mm_mullo_epi16( var1_r, var2_r )
#define pv_set1_int(number) _mm_set1_epi32( number )

//conversion functions (round or pad)
#define pv_cvt_real_to_int(real_vec) _mm_cvtpd_epi32( real_vec )
#define pv_cvt_int_to_real(int_vec) _mm_cvtepi32_pd( int_vec )

//---------------------------------------------
// Field interpolation for true 2D yz pusher.
// This may be ugly, but it smokes the serial version.

#define INTERP_FIELD_YZ(F_ENUM, indx2, indx3, outer_coeff, inner_coeff, field) { \
    pvReal field##tmp1, field##tmp2, field##tmp3;			\
    pvReal field##_in_A1, field##_in_B1;	\
    pvReal field##_in_A2, field##_in_B2;	\
    pvReal field##_in_A3, field##_in_B3;	\
    pvInt field##off1, field##off2, field##off3;			\
    VEC_OFF(field##off1,indx2##mns1,indx3##mns1);			\
    VEC_OFF(field##off2,indx2,indx3##mns1);				\
    VEC_OFF(field##off3,indx2##pls1,indx3##mns1);			\
									\
    field##_in_A1.r = pv_load1_real(F_ENUM##point+field##off1.v[0]);		\
    field##_in_B1.r = pv_load1_real(F_ENUM##point+field##off1.v[1]);		\
    field##tmp1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_B1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##tmp1.r);	\
									\
    field##_in_A2.r = pv_load1_real(F_ENUM##point+field##off2.v[0]);		\
    field##_in_B2.r = pv_load1_real(F_ENUM##point+field##off2.v[1]);		\
    field##tmp2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_B2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##tmp2.r);	\
									\
    field##_in_A3.r = pv_load1_real(F_ENUM##point+field##off3.v[0]);		\
    field##_in_B3.r = pv_load1_real(F_ENUM##point+field##off3.v[1]);		\
    field##tmp3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_B3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##tmp3.r);	\
    									\
    VEC_OFF(field##off1,indx2##mns1,indx3);				\
    VEC_OFF(field##off2,indx2,indx3);					\
    VEC_OFF(field##off3,indx2##pls1,indx3);				\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_mul_real(outer_coeff##mz.r, field##tmp1.r);		\
    									\
    field##_in_A1.r = pv_load1_real(F_ENUM##point+field##off1.v[0]);		\
    field##_in_B1.r = pv_load1_real(F_ENUM##point+field##off1.v[1]);		\
    field##tmp1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_B1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##tmp1.r);		\
									\
    field##_in_A2.r = pv_load1_real(F_ENUM##point+field##off2.v[0]);		\
    field##_in_B2.r = pv_load1_real(F_ENUM##point+field##off2.v[1]);		\
    field##tmp2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_B2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##tmp2.r);		\
									\
    field##_in_A3.r = pv_load1_real(F_ENUM##point+field##off3.v[0]);		\
    field##_in_B3.r = pv_load1_real(F_ENUM##point+field##off3.v[1]);		\
    field##tmp3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_B3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##tmp3.r);	\
									\
    VEC_OFF(field##off1,indx2##mns1,indx3##pls1);			\
    VEC_OFF(field##off2,indx2,indx3##pls1);				\
    VEC_OFF(field##off3,indx2##pls1,indx3##pls1);				\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##Oz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
    									\
    field##_in_A1.r = pv_load1_real(F_ENUM##point+field##off1.v[0]);		\
    field##_in_B1.r = pv_load1_real(F_ENUM##point+field##off1.v[1]);		\
    field##tmp1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_B1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##tmp1.r);		\
    									\
    field##_in_A2.r = pv_load1_real(F_ENUM##point+field##off2.v[0]);		\
    field##_in_B2.r = pv_load1_real(F_ENUM##point+field##off2.v[1]);		\
    field##tmp2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_B2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##tmp2.r);	\
    									\
    field##_in_A3.r = pv_load1_real(F_ENUM##point+field##off3.v[0]);		\
    field##_in_B3.r = pv_load1_real(F_ENUM##point+field##off3.v[1]);		\
    field##tmp3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_B3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##tmp3.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##lz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
  }


#else // Single Precision
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

// real functions
#define pv_add_real(var1_r, var2_r) _mm_add_ps( var1_r, var2_r )
#define pv_sub_real(var1_r, var2_r) _mm_sub_ps( var1_r, var2_r )
#define pv_mul_real(var1_r, var2_r) _mm_mul_ps( var1_r, var2_r )
#define pv_div_real(var1_r, var2_r) _mm_div_ps( var1_r, var2_r )
#define pv_sqrt_real(var1_r) _mm_sqrt_ps( var1_r)
#define pv_set1_real(number) _mm_set1_ps( number )
#define pv_load_real(float_ad) _mm_load_ps( float_ad )
#define pv_load1_real(float_point) _mm_load_ss( float_point )
#define pv_unpacklo_real( var1_r, var2_r ) _mm_unpacklo_ps( var1_r, var2_r)

// int functions
#define pv_add_int(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r )
#define pv_sub_int(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r )
#define pv_mul_int(var1_r, var2_r) _mm_mullo_epi16( var1_r, var2_r )
#define pv_set1_int(number) _mm_set1_epi32( number )

//conversion functions (round or pad)
#define pv_cvt_real_to_int(real_vec) _mm_cvtps_epi32( real_vec )
#define pv_cvt_int_to_real(int_vec) _mm_cvtepi32_ps( int_vec )

//---------------------------------------------
// Field interpolation for true 2D yz pusher.
// This may be ugly, but it smokes the serial version.

#define INTERP_FIELD_YZ(F_ENUM, indx2, indx3, outer_coeff, inner_coeff, field) { \
    pvReal field##tmp1, field##tmp2, field##tmp3;		\
    pvReal field##shuffAC_1, field##shuffBD_1, field##shuffAC_2, field##shuffBD_2, field##shuffAC_3, field##shuffBD_3; \
    pvReal field##_in_A1, field##_in_B1, field##_in_C1, field##_in_D1;	\
    pvReal field##_in_A2, field##_in_B2, field##_in_C2, field##_in_D2;	\
    pvReal field##_in_A3, field##_in_B3, field##_in_C3, field##_in_D3;	\
    pvInt field##off1, field##off2, field##off3;			\
    VEC_OFF(field##off1,indx2##mns1,indx3##mns1);			\
    VEC_OFF(field##off2,indx2,indx3##mns1);				\
    VEC_OFF(field##off3,indx2##pls1,indx3##mns1);			\
									\
    field##_in_A1.r = pv_load1_real(F_ENUM##point+field##off1.v[0]);		\
    field##_in_B1.r = pv_load1_real(F_ENUM##point+field##off1.v[1]);		\
    field##_in_C1.r = pv_load1_real(F_ENUM##point+field##off1.v[2]);		\
    field##_in_D1.r = pv_load1_real(F_ENUM##point+field##off1.v[3]);		\
    field##shuffAC_1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_C1.r); \
    field##shuffBD_1.r = pv_unpacklo_real( field##_in_B1.r, field##_in_D1.r); \
    field##tmp1.r = pv_unpacklo_real( field##shuffAC_1.r, field##shuffBD_1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##tmp1.r);	\
									\
    field##_in_A2.r = pv_load1_real(F_ENUM##point+field##off2.v[0]);		\
    field##_in_B2.r = pv_load1_real(F_ENUM##point+field##off2.v[1]);		\
    field##_in_C2.r = pv_load1_real(F_ENUM##point+field##off2.v[2]);		\
    field##_in_D2.r = pv_load1_real(F_ENUM##point+field##off2.v[3]);		\
    field##shuffAC_2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_C2.r); \
    field##shuffBD_2.r = pv_unpacklo_real( field##_in_B2.r, field##_in_D2.r); \
    field##tmp2.r = pv_unpacklo_real( field##shuffAC_2.r, field##shuffBD_2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##tmp2.r);	\
									\
    field##_in_A3.r = pv_load1_real(F_ENUM##point+field##off3.v[0]);		\
    field##_in_B3.r = pv_load1_real(F_ENUM##point+field##off3.v[1]);		\
    field##_in_C3.r = pv_load1_real(F_ENUM##point+field##off3.v[2]);		\
    field##_in_D3.r = pv_load1_real(F_ENUM##point+field##off3.v[3]);		\
    field##shuffAC_3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_C3.r); \
    field##shuffBD_3.r = pv_unpacklo_real( field##_in_B3.r, field##_in_D3.r); \
    field##tmp3.r = pv_unpacklo_real( field##shuffAC_3.r, field##shuffBD_3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##tmp3.r);	\
    									\
    VEC_OFF(field##off1,indx2##mns1,indx3);				\
    VEC_OFF(field##off2,indx2,indx3);					\
    VEC_OFF(field##off3,indx2##pls1,indx3);				\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_mul_real(outer_coeff##mz.r, field##tmp1.r);		\
    									\
    field##_in_A1.r = pv_load1_real(F_ENUM##point+field##off1.v[0]);		\
    field##_in_B1.r = pv_load1_real(F_ENUM##point+field##off1.v[1]);		\
    field##_in_C1.r = pv_load1_real(F_ENUM##point+field##off1.v[2]);		\
    field##_in_D1.r = pv_load1_real(F_ENUM##point+field##off1.v[3]);		\
    field##shuffAC_1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_C1.r); \
    field##shuffBD_1.r = pv_unpacklo_real( field##_in_B1.r, field##_in_D1.r); \
    field##tmp1.r = pv_unpacklo_real( field##shuffAC_1.r, field##shuffBD_1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##tmp1.r);		\
									\
    field##_in_A2.r = pv_load1_real(F_ENUM##point+field##off2.v[0]);		\
    field##_in_B2.r = pv_load1_real(F_ENUM##point+field##off2.v[1]);		\
    field##_in_C2.r = pv_load1_real(F_ENUM##point+field##off2.v[2]);		\
    field##_in_D2.r = pv_load1_real(F_ENUM##point+field##off2.v[3]);		\
    field##shuffAC_2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_C2.r); \
    field##shuffBD_2.r = pv_unpacklo_real( field##_in_B2.r, field##_in_D2.r); \
    field##tmp2.r = pv_unpacklo_real( field##shuffAC_2.r, field##shuffBD_2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##tmp2.r);		\
									\
    field##_in_A3.r = pv_load1_real(F_ENUM##point+field##off3.v[0]);		\
    field##_in_B3.r = pv_load1_real(F_ENUM##point+field##off3.v[1]);		\
    field##_in_C3.r = pv_load1_real(F_ENUM##point+field##off3.v[2]);		\
    field##_in_D3.r = pv_load1_real(F_ENUM##point+field##off3.v[3]);		\
    field##shuffAC_3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_C3.r); \
    field##shuffBD_3.r = pv_unpacklo_real( field##_in_B3.r, field##_in_D3.r); \
    field##tmp3.r = pv_unpacklo_real( field##shuffAC_3.r, field##shuffBD_3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##tmp3.r);	\
									\
    VEC_OFF(field##off1,indx2##mns1,indx3##pls1);			\
    VEC_OFF(field##off2,indx2,indx3##pls1);				\
    VEC_OFF(field##off3,indx2##pls1,indx3##pls1);				\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##Oz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
    									\
    field##_in_A1.r = pv_load1_real(F_ENUM##point+field##off1.v[0]);		\
    field##_in_B1.r = pv_load1_real(F_ENUM##point+field##off1.v[1]);		\
    field##_in_C1.r = pv_load1_real(F_ENUM##point+field##off1.v[2]);		\
    field##_in_D1.r = pv_load1_real(F_ENUM##point+field##off1.v[3]);		\
    field##shuffAC_1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_C1.r); \
    field##shuffBD_1.r = pv_unpacklo_real( field##_in_B1.r, field##_in_D1.r); \
    field##tmp1.r = pv_unpacklo_real( field##shuffAC_1.r, field##shuffBD_1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##my.r, field##tmp1.r);		\
    									\
    field##_in_A2.r = pv_load1_real(F_ENUM##point+field##off2.v[0]);		\
    field##_in_B2.r = pv_load1_real(F_ENUM##point+field##off2.v[1]);		\
    field##_in_C2.r = pv_load1_real(F_ENUM##point+field##off2.v[2]);		\
    field##_in_D2.r = pv_load1_real(F_ENUM##point+field##off2.v[3]);		\
    field##shuffAC_2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_C2.r); \
    field##shuffBD_2.r = pv_unpacklo_real( field##_in_B2.r, field##_in_D2.r); \
    field##tmp2.r = pv_unpacklo_real( field##shuffAC_2.r, field##shuffBD_2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##Oy.r, field##tmp2.r);	\
    									\
    field##_in_A3.r = pv_load1_real(F_ENUM##point+field##off3.v[0]);		\
    field##_in_B3.r = pv_load1_real(F_ENUM##point+field##off3.v[1]);		\
    field##_in_C3.r = pv_load1_real(F_ENUM##point+field##off3.v[2]);		\
    field##_in_D3.r = pv_load1_real(F_ENUM##point+field##off3.v[3]);		\
    field##shuffAC_3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_C3.r); \
    field##shuffBD_3.r = pv_unpacklo_real( field##_in_B3.r, field##_in_D3.r); \
    field##tmp3.r = pv_unpacklo_real( field##shuffAC_3.r, field##shuffBD_3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##ly.r, field##tmp3.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##lz.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
  }

#endif

#endif
