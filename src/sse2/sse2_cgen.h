
#ifndef SSE2_CGEN_H
#define SSE2_CGEN_H

#include <xmmintrin.h>
#include <emmintrin.h>


#if SSE2_DOUBLE

#define VEC_SIZE 2
typedef double sse2_real;

union packed_vector{
  __m128d r;
  sse2_real v[2] __attribute__ ((aligned (16))); //FIXME : Might break for any non gcc
} ;

union packed_int{
  __m128i r;
  int v[2] __attribute__ ((aligned (16))); // this is a little strange, but only the lower
                                           // two values will ever be used.
} ; 

// real functions
#define pv_real_ADD(var1_r, var2_r) _mm_add_pd( var1_r, var2_r )
#define pv_real_SUB(var1_r, var2_r) _mm_sub_pd( var1_r, var2_r )
#define pv_real_MUL(var1_r, var2_r) _mm_mul_pd( var1_r, var2_r )
#define pv_real_DIV(var1_r, var2_r) _mm_div_pd( var1_r, var2_r )
#define pv_real_SQRT(var1_r) _mm_sqrt_pd( var1_r)
#define pv_real_SET1(number) _mm_set1_pd( number )

// int functions
#define pv_int_ADD(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r )
#define pv_int_SUB(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r )
#define pv_int_SET1(number) _mm_set1_epi32( number )

//conversion functions (round or pad)
#define pv_real_to_int_CVT(real_vec) _mm_cvtpd_epi32( real_vec )
#define pv_int_to_real_CVT(int_vec) _mm_cvtepi32_pd( int_vec )


// These are loop unwinding macros for some pretty nasty sections


#define LOAD_PART( source, loop_iter ) {		\
  							\
    p.xi.v[0] = source->part[loop_iter].xi;		\
    p.yi.v[0] = source->part[loop_iter].yi;		\
    p.zi.v[0] = source->part[loop_iter].zi;		\
    p.pxi.v[0] = source->part[loop_iter].pxi;		\
    p.pyi.v[0] = source->part[loop_iter].pyi;		\
    p.pzi.v[0] = source->part[loop_iter].pzi;		\
    p.qni.v[0] = source->part[loop_iter].qni;		\
    p.mni.v[0] = source->part[loop_iter].mni;		\
    p.wni.v[0] = source->part[loop_iter].wni;		\
							\
    p.xi.v[1] = source->part[loop_iter + 1].xi;		\
    p.yi.v[1] = source->part[loop_iter + 1].yi;		\
    p.zi.v[1] = source->part[loop_iter + 1].zi;		\
    p.pxi.v[1] = source->part[loop_iter + 1].pxi;		\
    p.pyi.v[1] = source->part[loop_iter + 1].pyi;		\
    p.pzi.v[1] = source->part[loop_iter + 1].pzi;		\
    p.qni.v[1] = source->part[loop_iter + 1].qni;		\
    p.mni.v[1] = source->part[loop_iter + 1].mni;		\
    p.wni.v[1] = source->part[loop_iter + 1].wni;		\
  }


#define STORE_PART_XP(out, loop_iter) {		\
  (out->part[loop_iter]).xi = p.xi.v[0];		\
  (out->part[loop_iter]).yi = p.yi.v[0];		\
  (out->part[loop_iter]).zi = p.zi.v[0];		\
  (out->part[loop_iter]).pxi = p.pxi.v[0];	\
  (out->part[loop_iter]).pyi = p.pyi.v[0];	\
  (out->part[loop_iter]).pzi = p.pzi.v[0];	\
						\
  (out->part[loop_iter+1]).xi = p.xi.v[1];	\
  (out->part[loop_iter+1]).yi = p.yi.v[1];	\
  (out->part[loop_iter+1]).zi = p.zi.v[1];	\
  (out->part[loop_iter+1]).pxi = p.pxi.v[1];	\
  (out->part[loop_iter+1]).pyi = p.pyi.v[1];	\
  (out->part[loop_iter+1]).pzi = p.pzi.v[1];	\
						\
  }

#define INTERP_FIELD_YZ(F_ENUM, indx1, indx2, indx3, outer_coeff, inner_coeff, field) { \
									\
    pvReal field##_in, field##tmp1, field##tmp2, field##tmp3;		\
									\
									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##mns1.v[0],indx3##mns1.v[0]);	\
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##mns1.v[1],indx3##mns1.v[1]);	\
    field##tmp1.r = pv_real_MUL(inner_coeff##my.r, field##_in.r);		\
									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2.v[0],indx3##mns1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2.v[1],indx3##mns1.v[1]); \
    field##tmp2.r = pv_real_MUL(inner_coeff##Oy.r, field##_in.r);	\
									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##pls1.v[0],indx3##mns1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##pls1.v[1],indx3##mns1.v[1]); \
    field##tmp3.r = pv_real_MUL(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_real_MUL(outer_coeff##mz.r, field##tmp1.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##mns1.v[0],indx3.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##mns1.v[1],indx3.v[1]); \
    field##tmp1.r = pv_real_MUL(inner_coeff##my.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2.v[0],indx3.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2.v[1],indx3.v[1]); \
    field##tmp2.r = pv_real_MUL(inner_coeff##Oy.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##pls1.v[0],indx3.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##pls1.v[1],indx3.v[1]); \
    field##tmp3.r = pv_real_MUL(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_real_MUL(outer_coeff##Oz.r, field##tmp1.r);	\
    field.r = pv_real_ADD(field.r, field##tmp1.r);			\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##mns1.v[0],indx3##pls1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##mns1.v[1],indx3##pls1.v[1]); \
    field##tmp1.r = pv_real_MUL(inner_coeff##my.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2.v[0],indx3##pls1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2.v[1],indx3##pls1.v[1]); \
    field##tmp2.r = pv_real_MUL(inner_coeff##Oy.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##pls1.v[0],indx3##pls1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##pls1.v[1],indx3##pls1.v[1]); \
    field##tmp3.r = pv_real_MUL(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_real_MUL(outer_coeff##lz.r, field##tmp1.r);	\
    field.r = pv_real_ADD(field.r, field##tmp1.r);			\
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
#define pv_real_ADD(var1_r, var2_r) _mm_add_ps( var1_r, var2_r )
#define pv_real_SUB(var1_r, var2_r) _mm_sub_ps( var1_r, var2_r )
#define pv_real_MUL(var1_r, var2_r) _mm_mul_ps( var1_r, var2_r )
#define pv_real_DIV(var1_r, var2_r) _mm_div_ps( var1_r, var2_r )
#define pv_real_SQRT(var1_r) _mm_sqrt_ps( var1_r)
#define pv_real_SET1(number) _mm_set1_ps( number )

// int functions
#define pv_int_ADD(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r )
#define pv_int_SUB(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r )
#define pv_int_SET1(number) _mm_set1_epi32( number )

//conversion functions (round or pad)
#define pv_real_to_int_CVT(real_vec) _mm_cvtps_epi32( real_vec )
#define pv_int_to_real_CVT(int_vec) _mm_cvtepi32_ps( int_vec )

// These are loop unwinding macros for some pretty nasty sections


#define LOAD_PART( source, loop_iter ) {		\
							\
    p.xi.v[0] = source->part[loop_iter].xi;		\
    p.yi.v[0] = source->part[loop_iter].yi;		\
    p.zi.v[0] = source->part[loop_iter].zi;		\
    p.pxi.v[0] = source->part[loop_iter].pxi;		\
    p.pyi.v[0] = source->part[loop_iter].pyi;		\
    p.pzi.v[0] = source->part[loop_iter].pzi;		\
    p.qni.v[0] = source->part[loop_iter].qni;		\
    p.mni.v[0] = source->part[loop_iter].mni;		\
    p.wni.v[0] = source->part[loop_iter].wni;		\
							\
    p.xi.v[1] = source->part[loop_iter + 1].xi;		\
    p.yi.v[1] = source->part[loop_iter + 1].yi;		\
    p.zi.v[1] = source->part[loop_iter + 1].zi;		\
    p.pxi.v[1] = source->part[loop_iter + 1].pxi;		\
    p.pyi.v[1] = source->part[loop_iter + 1].pyi;		\
    p.pzi.v[1] = source->part[loop_iter + 1].pzi;		\
    p.qni.v[1] = source->part[loop_iter + 1].qni;		\
    p.mni.v[1] = source->part[loop_iter + 1].mni;		\
    p.wni.v[1] = source->part[loop_iter + 1].wni;		\
							\
    p.xi.v[2] = source->part[loop_iter + 2].xi;		\
    p.yi.v[2] = source->part[loop_iter + 2].yi;		\
    p.zi.v[2] = source->part[loop_iter + 2].zi;		\
    p.pxi.v[2] = source->part[loop_iter + 2].pxi;		\
    p.pyi.v[2] = source->part[loop_iter + 2].pyi;		\
    p.pzi.v[2] = source->part[loop_iter + 2].pzi;		\
    p.qni.v[2] = source->part[loop_iter + 2].qni;		\
    p.mni.v[2] = source->part[loop_iter + 2].mni;		\
    p.wni.v[2] = source->part[loop_iter + 2].wni;		\
							\
    p.xi.v[3] = source->part[loop_iter + 3].xi;		\
    p.yi.v[3] = source->part[loop_iter + 3].yi;		\
    p.zi.v[3] = source->part[loop_iter + 3].zi;		\
    p.pxi.v[3] = source->part[loop_iter + 3].pxi;		\
    p.pyi.v[3] = source->part[loop_iter + 3].pyi;		\
    p.pzi.v[3] = source->part[loop_iter + 3].pzi;		\
    p.qni.v[3] = source->part[loop_iter + 3].qni;		\
    p.mni.v[3] = source->part[loop_iter + 3].mni;		\
    p.wni.v[3] = source->part[loop_iter + 3].wni;		\
  }

#define STORE_PART_XP(out, loop_iter) {		\
  (out->part[loop_iter]).xi = p.xi.v[0];		\
  (out->part[loop_iter]).yi = p.yi.v[0];		\
  (out->part[loop_iter]).zi = p.zi.v[0];		\
  (out->part[loop_iter]).pxi = p.pxi.v[0];	\
  (out->part[loop_iter]).pyi = p.pyi.v[0];	\
  (out->part[loop_iter]).pzi = p.pzi.v[0];	\
						\
  (out->part[loop_iter+1]).xi = p.xi.v[1];	\
  (out->part[loop_iter+1]).yi = p.yi.v[1];	\
  (out->part[loop_iter+1]).zi = p.zi.v[1];	\
  (out->part[loop_iter+1]).pxi = p.pxi.v[1];	\
  (out->part[loop_iter+1]).pyi = p.pyi.v[1];	\
  (out->part[loop_iter+1]).pzi = p.pzi.v[1];	\
						\
  (out->part[loop_iter+2]).xi = p.xi.v[2];	\
  (out->part[loop_iter+2]).yi = p.yi.v[2];	\
  (out->part[loop_iter+2]).zi = p.zi.v[2];	\
  (out->part[loop_iter+2]).pxi = p.pxi.v[2];	\
  (out->part[loop_iter+2]).pyi = p.pyi.v[2];	\
  (out->part[loop_iter+2]).pzi = p.pzi.v[2];	\
						\
  (out->part[loop_iter+3]).xi = p.xi.v[3];	\
  (out->part[loop_iter+3]).yi = p.yi.v[3];	\
  (out->part[loop_iter+3]).zi = p.zi.v[3];	\
  (out->part[loop_iter+3]).pxi = p.pxi.v[3];	\
  (out->part[loop_iter+3]).pyi = p.pyi.v[3];	\
  (out->part[loop_iter+3]).pzi = p.pzi.v[3];	\
  }

// The below is here because I'd like to try and generalize it to any 2D plane
#define INTERP_FIELD_YZ(F_ENUM, indx1, indx2, indx3, outer_coeff, inner_coeff, field) { \
  									\
    pvReal field##_in, field##tmp1, field##tmp2, field##tmp3;		\
									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##mns1.v[0],indx3##mns1.v[0]);	\
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##mns1.v[1],indx3##mns1.v[1]);	\
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2##mns1.v[2],indx3##mns1.v[2]);	\
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2##mns1.v[3],indx3##mns1.v[3]);	\
    field##tmp1.r = pv_real_MUL(inner_coeff##my.r, field##_in.r);		\
									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2.v[0],indx3##mns1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2.v[1],indx3##mns1.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2.v[2],indx3##mns1.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2.v[3],indx3##mns1.v[3]); \
    field##tmp2.r = pv_real_MUL(inner_coeff##Oy.r, field##_in.r);		\
									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##pls1.v[0],indx3##mns1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##pls1.v[1],indx3##mns1.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2##pls1.v[2],indx3##mns1.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2##pls1.v[3],indx3##mns1.v[3]);	\
    field##tmp3.r = pv_real_MUL(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_real_MUL(outer_coeff##mz.r, field##tmp1.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##mns1.v[0],indx3.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##mns1.v[1],indx3.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2##mns1.v[2],indx3.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2##mns1.v[3],indx3.v[3]); \
    field##tmp1.r = pv_real_MUL(inner_coeff##my.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2.v[0],indx3.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2.v[1],indx3.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2.v[2],indx3.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2.v[3],indx3.v[3]); \
    field##tmp2.r = pv_real_MUL(inner_coeff##Oy.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##pls1.v[0],indx3.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##pls1.v[1],indx3.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2##pls1.v[2],indx3.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2##pls1.v[3],indx3.v[3]); \
    field##tmp3.r = pv_real_MUL(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_real_MUL(outer_coeff##Oz.r, field##tmp1.r);	\
    field.r = pv_real_ADD(field.r, field##tmp1.r);			\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##mns1.v[0],indx3##pls1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##mns1.v[1],indx3##pls1.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2##mns1.v[2],indx3##pls1.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2##mns1.v[3],indx3##pls1.v[3]);	\
    field##tmp1.r = pv_real_MUL(inner_coeff##my.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2.v[0],indx3##pls1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2.v[1],indx3##pls1.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2.v[2],indx3##pls1.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2.v[3],indx3##pls1.v[3]); \
    field##tmp2.r = pv_real_MUL(inner_coeff##Oy.r, field##_in.r);		\
    									\
    field##_in.v[0] = CF3(F_ENUM, indx1.v[0], indx2##pls1.v[0],indx3##pls1.v[0]); \
    field##_in.v[1] = CF3(F_ENUM, indx1.v[1], indx2##pls1.v[1],indx3##pls1.v[1]); \
    field##_in.v[2] = CF3(F_ENUM, indx1.v[2], indx2##pls1.v[2],indx3##pls1.v[2]); \
    field##_in.v[3] = CF3(F_ENUM, indx1.v[3], indx2##pls1.v[3],indx3##pls1.v[3]);	\
    field##tmp3.r = pv_real_MUL(inner_coeff##ly.r, field##_in.r);		\
    									\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_real_ADD(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_real_MUL(outer_coeff##lz.r, field##tmp1.r);	\
    field.r = pv_real_ADD(field.r, field##tmp1.r);			\
  }

#endif

#endif
