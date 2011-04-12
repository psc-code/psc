
#ifndef SIMD_WRAP_H
#define SIMD_WRAP_H

#include <xmmintrin.h>
#include <emmintrin.h>

//////////////////////////////////////////////////////////////////////
/// Slow (serial) but correct version of 2D field interpolation
///
/// FIXME: gh2 and gh3 should switch positions in the argument list

#define F3YZ(m, j,k) F3_SSE2(pf, m, 0,j,k)

#define INTERP_FIELD_YZ_SLOW(m, jl2, jl3, gh3, gh2, var) do {		\
    for (int c = 0; c < VEC_SIZE; c++) {				\
      var.v[c]								\
	= (gh3##mz.v[c]*( gh2##my.v[c]*F3YZ(m,jl2.v[c]-1,jl3.v[c]-1) +	\
			  gh2##Oy.v[c]*F3YZ(m,jl2.v[c]  ,jl3.v[c]-1) +	\
			  gh2##ly.v[c]*F3YZ(m,jl2.v[c]+1,jl3.v[c]-1)) +	\
	   gh3##Oz.v[c]*( gh2##my.v[c]*F3YZ(m,jl2.v[c]-1,jl3.v[c]  ) +	\
			  gh2##Oy.v[c]*F3YZ(m,jl2.v[c]  ,jl3.v[c]  ) +	\
			  gh2##ly.v[c]*F3YZ(m,jl2.v[c]+1,jl3.v[c]  )) +	\
	   gh3##lz.v[c]*( gh2##my.v[c]*F3YZ(m,jl2.v[c]-1,jl3.v[c]+1) +	\
			  gh2##Oy.v[c]*F3YZ(m,jl2.v[c]  ,jl3.v[c]+1) +	\
			  gh2##ly.v[c]*F3YZ(m,jl2.v[c]+1,jl3.v[c]+1)));	\
    }									\
  } while (0)

#define F3XZ(m, j,k) F3_SSE2(pf, m, j,psc.ilo[1],k)
#define INTERP_FIELD_XZ_SLOW(m, jl1, jl3, gh3, gh1, var) do {		\
    for (int c = 0; c < VEC_SIZE; c++) {				\
      var.v[c]								\
	= (gh3##mz.v[c]*( gh1##mx.v[c]*F3XZ(m,jl1.v[c]-1,jl3.v[c]-1) +	\
			  gh1##Ox.v[c]*F3XZ(m,jl1.v[c]  ,jl3.v[c]-1) +	\
			  gh1##lx.v[c]*F3XZ(m,jl1.v[c]+1,jl3.v[c]-1)) +	\
	   gh3##Oz.v[c]*( gh1##mx.v[c]*F3XZ(m,jl1.v[c]-1,jl3.v[c]  ) +	\
			  gh1##Ox.v[c]*F3XZ(m,jl1.v[c]  ,jl3.v[c]  ) +	\
			  gh1##lx.v[c]*F3XZ(m,jl1.v[c]+1,jl3.v[c]  )) +	\
	   gh3##lz.v[c]*( gh1##mx.v[c]*F3XZ(m,jl1.v[c]-1,jl3.v[c]+1) +	\
			  gh1##Ox.v[c]*F3XZ(m,jl1.v[c]  ,jl3.v[c]+1) +	\
			  gh1##lx.v[c]*F3XZ(m,jl1.v[c]+1,jl3.v[c]+1)));	\
    }									\
  } while (0)


//---------------------------------------------
/// Calculates the offset for accessing pointers to 
/// flattened (2-D YZ only) fields using vector operations

#define FIELD_OFF_yz(jy, jz)			\
  ((((jz) - psc.ilg[2])				\
    *psc.img[1] + ((jy)-psc.ilg[1]))*psc.img[0])

#define FIELD_OFF_xz(jx,jz)			\
  (((jz) - psc.ilg[2])					\
    *psc.img[1]*psc.img[0] + ((jx) - psc.ilg[0]))

#if SSE2_DOUBLE

union packed_vector{
  __m128d r;
  sse2_real v[2] __attribute__ ((aligned (32))); //FIXME : Might break for any non gcc
} ;

union packed_int{
  __m128i r;
  int v[4] __attribute__ ((aligned (16))); 
                                           
} ; 

//-----------------------------
// real functions

//Arithmetic
#define pv_add_real(var1_r, var2_r) _mm_add_pd( var1_r, var2_r )
#define pv_sub_real(var1_r, var2_r) _mm_sub_pd( var1_r, var2_r )
#define pv_mul_real(var1_r, var2_r) _mm_mul_pd( var1_r, var2_r )
#define pv_div_real(var1_r, var2_r) _mm_div_pd( var1_r, var2_r )
#define pv_sqrt_real(var1_r) _mm_sqrt_pd( var1_r)

//Loads and sets
#define pv_set1_real(number) _mm_set1_pd( number )
#define pv_loads_real(pointer) _mm_load_sd( pointer )
#define pv_load_real(mem_loc) _mm_load_pd( mem_loc )
#define pv_loadu_real(mem_loc) _mm_loadu_pd( mem_loc )

//Stores
#define pv_store_real(mem_loc, var_r) _mm_store_pd( (mem_loc), var_r)
#define pv_storeu_real(mem_loc, var_r) _mm_storeu_pd( (mem_loc), var_r)

// Unpacks and shuffles
#define pv_unpacklo_real(var1_r, var2_r) _mm_unpacklo_pd( var1_r, var2_r )
#define pv_unpackhi_real(var1_r, var2_r) _mm_unpackhi_pd( var1_r, var2_r )

// int functions
#define pv_add_int(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r )
#define pv_sub_int(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r )
#define pv_mul_uint(var1_r, var2_r) func_mul_epu32( var1_r, var2_r )
#define pv_set1_int(number) _mm_set1_epi32( number )
#define pv_extract_int(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

//conversion functions (round or pad)
#define pv_cvt_real_to_int(real_vec) _mm_cvtpd_epi32( real_vec )
#define pv_cvt_int_to_real(int_vec) _mm_cvtepi32_pd( int_vec )

//----------------------------------
// Load particles from memory. This doesn't seem to hurt, and
// on a more modern processor (fishercat) it seems to give a factor of 
// two speed increase over serial loading
 
#define LOAD_PARTS(p_struct) {						\
    pvReal __A, __B, __C, __D, __E, __F, __G, __H;			\
    __A.r = pv_loadu_real(&pp->particles[n].xi);			\
    __E.r = pv_loadu_real(&pp->particles[n+1].xi);			\
									\
    __B.r = pv_loadu_real(&pp->particles[n].zi);			\
    __F.r = pv_loadu_real(&pp->particles[n+1].zi);			\
									\
    p_struct.xi.r = pv_unpacklo_real(__A.r, __E.r);			\
    p_struct.yi.r = pv_unpackhi_real(__A.r, __E.r);			\
									\
    p_struct.zi.r = pv_unpacklo_real(__B.r, __F.r);			\
    p_struct.pxi.r = pv_unpackhi_real(__B.r, __F.r);			\
									\
    __C.r = pv_loadu_real(&pp->particles[n].pyi);			\
    __G.r = pv_loadu_real(&pp->particles[n+1].pyi);			\
									\
    __D.r = pv_loadu_real(&pp->particles[n].qni);			\
    __H.r = pv_loadu_real(&pp->particles[n+1].qni);			\
									\
    p_struct.pyi.r = pv_unpacklo_real(__C.r, __G.r);			\
    p_struct.pzi.r = pv_unpackhi_real(__C.r, __G.r);			\
									\
    p_struct.qni.r = pv_unpacklo_real(__D.r, __H.r);			\
    p_struct.mni.r = pv_unpackhi_real(__D.r, __H.r);			\
									\
    __A.r = pv_loads_real(&pp->particles[n].wni);			\
    __E.r = pv_loads_real(&pp->particles[n+1].wni);			\
									\
    p_struct.wni.r = pv_unpacklo_real(__A.r, __E.r);			\
  }

//----------------------------------
// Store particles to memory. This either has no effect,
// or produces a performance boost. I really can't tell,
// the timing jitters too much at the moment.

#define STORE_PARTS(p_struct) {				\
    pvReal __A, __B, __C, __E, __F, __G;	\
    __A.r = pv_unpacklo_real(p_struct.xi.r, p_struct.yi.r);	\
						\
    pv_storeu_real(&pp->particles[n].xi, __A.r);	\
						\
    __B.r = pv_unpacklo_real(p_struct.zi.r, p_struct.pxi.r);	\
    __C.r = pv_unpacklo_real(p_struct.pyi.r, p_struct.pzi.r);	\
						\
    pv_storeu_real(&pp->particles[n].zi, __B.r);	\
						\
    __E.r = pv_unpackhi_real(p_struct.xi.r, p_struct.yi.r);		\
    __F.r = pv_unpackhi_real(p_struct.zi.r, p_struct.pxi.r);		\
						\
    pv_storeu_real(&pp->particles[n].pyi, __C.r);	\
						\
    __G.r = pv_unpackhi_real(p_struct.pyi.r, p_struct.pzi.r);	\
						\
    pv_storeu_real(&pp->particles[n+1].xi, __E.r);	\
    pv_storeu_real(&pp->particles[n+1].zi, __F.r);		\
    pv_storeu_real(&pp->particles[n+1].pyi, __G.r);		\
  }


//---------------------------------------------
// Field interpolation for true 2D yz pusher.
// This may be ugly, but it smokes the serial version.

#define INTERP_FIELD_2D_SIMD(inner_dir, outer_dir, F_ENUM, indx2, indx3, outer_coeff, inner_coeff, field) { \
    pvReal field##tmp1, field##tmp2, field##tmp3;			\
    pvReal field##_in_A1, field##_in_B1;				\
    pvReal field##_in_A2, field##_in_B2;				\
    pvReal field##_in_A3, field##_in_B3;				\
    int i2_0 = indx2.v[0],						\
      i2_1 = indx2.v[1],						\
      i3_0 = indx3.v[0],						\
      i3_1 = indx3.v[1];						\
									\
    field##_in_A1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 - 1,i3_0 - 1)); \
    field##_in_B1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 - 1,i3_1 - 1)); \
    field##tmp1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_B1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##m##inner_dir.r, field##tmp1.r);	\
									\
    field##_in_A2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0,i3_0 - 1)); \
    field##_in_B2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1,i3_1 - 1)); \
    field##tmp2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_B2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##O##inner_dir.r, field##tmp2.r); \
    									\
    field##_in_A3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 + 1,i3_0 - 1)); \
    field##_in_B3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 + 1,i3_1 - 1)); \
    field##tmp3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_B3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##l##inner_dir.r, field##tmp3.r);	\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_mul_real(outer_coeff##m##outer_dir.r, field##tmp1.r);		\
    									\
    field##_in_A1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 - 1,i3_0)); \
    field##_in_B1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 - 1,i3_1)); \
    field##tmp1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_B1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##m##inner_dir.r, field##tmp1.r); \
									\
    field##_in_A2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0,i3_0)); \
    field##_in_B2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1,i3_1)); \
    field##tmp2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_B2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##O##inner_dir.r, field##tmp2.r);		\
									\
    field##_in_A3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 + 1,i3_0)); \
    field##_in_B3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 + 1,i3_1)); \
    field##tmp3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_B3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##l##inner_dir.r, field##tmp3.r);	\
									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##O##outer_dir.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
    									\
    field##_in_A1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 - 1,i3_0 +1)); \
    field##_in_B1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 - 1,i3_1 +1)); \
    field##tmp1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_B1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##m##inner_dir.r, field##tmp1.r);		\
    									\
    field##_in_A2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0,i3_0 +1)); \
    field##_in_B2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1,i3_1 +1)); \
    field##tmp2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_B2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##O##inner_dir.r, field##tmp2.r);	\
    									\
    field##_in_A3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 +1,i3_0 +1)); \
    field##_in_B3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 +1,i3_1 +1)); \
    field##tmp3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_B3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##l##inner_dir.r, field##tmp3.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##l##outer_dir.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
  }


#else // Single Precision

/// Union of real vector and array
union packed_vector{
  __m128 r; ///< Vector
  sse2_real v[4] __attribute__ ((aligned (16))); ///< Array
} ;

/// \FIXME The __attribute__ ((aligned())) statement may not function on certain non-gcc compilers. I think microsoft/intel compilers use delspec(). Unfortunately, the standard says nothing about these keywords. We may need to macro-ize alignment keywords.

/// Union of integer vector and integer array.
union packed_int{
  __m128i r; ///< Vector
  int v[4] __attribute__ ((aligned (16))); ///< Array
}; 

//---------------------------------------
// real functions

//Arithmetic
/// Real addition
#define pv_add_real(var1_r, var2_r) _mm_add_ps( var1_r, var2_r ) 
/// Real subtraction
#define pv_sub_real(var1_r, var2_r) _mm_sub_ps( var1_r, var2_r )
/// Real multiplication
#define pv_mul_real(var1_r, var2_r) _mm_mul_ps( var1_r, var2_r )
/// Real division
#define pv_div_real(var1_r, var2_r) _mm_div_ps( var1_r, var2_r ) 
/// Real sqrt
#define pv_sqrt_real(var1_r) _mm_sqrt_ps( var1_r) 

// Sets and Loads
/// Sets all members to the given const. float
#define pv_set1_real(number) _mm_set1_ps( number )  
/// loads an aligned vector from mem_loc
#define pv_load_real(mem_loc) _mm_load_ps( mem_loc ) 
/// loads an unaligned vector from mem_loc
#define pv_loadu_real(mem_loc) _mm_loadu_ps( mem_loc )
/// loads a single word from mem_loc into the lowest element and zeros the rest of the vector 
#define pv_loads_real(mem_loc) _mm_load_ss( mem_loc ) 

//Stores
/// stores a vector to aligned memory at mem_loc
#define pv_store_real(mem_loc, var_r) _mm_store_ps(mem_loc, var_r) 
/// stores a vector to unaligned memory at mem_loc
#define pv_storeu_real(mem_loc, var_r) _mm_storeu_ps(mem_loc, var_r) 

//Unpacks and shuffles
/// Interleave lower elements of given vectors.
#define pv_unpacklo_real( var1_r, var2_r ) _mm_unpacklo_ps( var1_r, var2_r) 
/// Interleave upper elements of given vectors.
#define pv_unpackhi_real( var1_r, var2_r ) _mm_unpackhi_ps( var1_r, var2_r) 


//---------------------------------------
// int functions

//Arithmetic
/// Integer addition
#define pv_add_int(var1_r, var2_r) _mm_add_epi32( var1_r, var2_r ) 
/// Integer subtraction.
#define pv_sub_int(var1_r, var2_r) _mm_sub_epi32( var1_r, var2_r ) 
/// Unsigned Integer multiplication. 
///
/// Calls the function func_mul_epu32() in psc_sse2.c
#define pv_mul_uint(var1_r, var2_r) func_mul_epu32( var1_r, var2_r )

// Sets and loads
/// Set all elements of an integer vector to the same value
#define pv_set1_int(number) _mm_set1_epi32( number ) 

//------------------------------------
//conversion functions (round or pad)

/// Round real to integer
#define pv_cvt_real_to_int(real_vec) _mm_cvtps_epi32( real_vec ) 
/// Promote from integer to real
#define pv_cvt_int_to_real(int_vec) _mm_cvtepi32_ps( int_vec ) 

/////////////////////////
/// Load particles from memory. 
///
/// On a modern processor packing data into vectors
/// using this shuffle/unpack method is significantly faster than doing single element
/// loads. The reason for this is somewhat complicated, but important to understand.
/// Most processors can only load either the entire vector or the lowest, 'preferred', 
/// element. That is to say, given a vector [a0, a1, a2, a3], there is an instruction to load
/// all four elements a0-a3 from memory, and an instruction to load just a0, but no instructions
/// exist to load a1-a3 individually. Therefore if one attempts to load a0-a3 from seperate memory
/// locations the compiler has to translate the code in a sub-optimal way. Generally it will 
/// copy each element from main memory onto the stack. Once the vector is assembled on the stack, 
/// the processor will either perform a load into a vector register or simply use the value of the
/// stack everytime it is requested. This method is slow. Shuffles/unpacks, while complicated to 
/// write, are significantly faster, and should be used whenever possible. A similar method is used
/// for field loading.

 
#define LOAD_PARTS(p_struct) {					\
    pvReal __A, __B, __C, __D, __E, __F, __G, __H, __alpha, __beta, __gamma, __delta;	\
    __A.r = pv_loadu_real(&pp->particles[n].xi);			\
    __C.r = pv_loadu_real(&pp->particles[n+2].xi);			\
    __B.r = pv_loadu_real(&pp->particles[n+1].xi);			\
    __D.r = pv_loadu_real(&pp->particles[n+3].xi);			\
								\
    __alpha.r = pv_unpacklo_real(__A.r, __C.r);			\
    __gamma.r = pv_unpackhi_real(__A.r, __C.r);			\
    __beta.r = pv_unpacklo_real(__B.r, __D.r);			\
    __delta.r = pv_unpackhi_real(__B.r, __D.r);			\
								\
    __E.r = pv_loadu_real(&pp->particles[n].pyi);			\
    __G.r = pv_loadu_real(&pp->particles[n+2].pyi);			\
								\
    p_struct.xi.r = pv_unpacklo_real(__alpha.r, __beta.r);			\
    p_struct.yi.r = pv_unpackhi_real(__alpha.r, __beta.r);			\
    p_struct.zi.r = pv_unpacklo_real(__gamma.r, __delta.r);		\
    p_struct.pxi.r = pv_unpackhi_real(__gamma.r, __delta.r);		\
								\
    __F.r = pv_loadu_real(&pp->particles[n+1].pyi);			\
    __H.r = pv_loadu_real(&pp->particles[n+3].pyi);			\
								\
    __gamma.r = pv_unpackhi_real(__E.r, __G.r);			\
    __alpha.r = pv_unpacklo_real(__E.r, __G.r);			\
    __beta.r = pv_unpacklo_real(__F.r, __H.r);			\
    __delta.r = pv_unpackhi_real(__F.r, __H.r);			\
								\
    p_struct.pyi.r = pv_unpacklo_real(__alpha.r, __beta.r);		\
    p_struct.pzi.r = pv_unpackhi_real(__alpha.r, __beta.r);		\
    p_struct.qni.r = pv_unpacklo_real(__gamma.r, __delta.r);		\
    p_struct.mni.r = pv_unpackhi_real(__gamma.r, __delta.r);		\
								\
    __A.r = pv_loads_real(&pp->particles[n].wni);			\
    __C.r = pv_loads_real(&pp->particles[n+2].wni);			\
    __B.r = pv_loads_real(&pp->particles[n+1].wni);			\
    __D.r = pv_loads_real(&pp->particles[n+3].wni);			\
								\
    __alpha.r = pv_unpacklo_real(__A.r, __C.r);			\
    __beta.r = pv_unpacklo_real(__B.r, __D.r);			\
    p_struct.wni.r = pv_unpacklo_real(__alpha.r, __beta.r);		\
  }

//----------------------------------
/// Store particles to memory. 

/// For single precision it doesn't appear to be worthwhile
/// to shuffle and do a vector store. Some tests have even
/// indicated it could hurt performance by occupying fp
/// units and decreasing pipelining efficiency, but by and large
/// the tests are inconclusive
#define STORE_PARTS(p_struct) {				\
  for(int m=0; m < VEC_SIZE; m++){		\
      (pp->particles[n+m]).xi = p_struct.xi.v[m];	\
      (pp->particles[n+m]).yi = p_struct.yi.v[m];	\
      (pp->particles[n+m]).zi = p_struct.zi.v[m];	\
      (pp->particles[n+m]).pxi = p_struct.pxi.v[m];	\
      (pp->particles[n+m]).pyi = p_struct.pyi.v[m];	\
      (pp->particles[n+m]).pzi = p_struct.pzi.v[m];	\
  }						\
}

////////////////////////
/// Field interpolation for true 2D yz pusher.
///
/// On a modern processor packing data into vectors
/// using this shuffle/unpack method is significantly faster than doing single element
/// loads. The reason for this is somewhat complicated, but important to understand.
/// Most processors can only load either the entire vector or the lowest, 'preferred', 
/// element. That is to say, given a vector [a0, a1, a2, a3], there is an instruction to load
/// all four elements a0-a3 from memory, and an instruction to load just a0, but no instructions
/// exist to load a1-a3 individually. Therefore if one attempts to load a0-a3 from seperate memory
/// locations the compiler has to translate the code in a sub-optimal way. Generally it will 
/// copy each element from main memory onto the stack. Once the vector is assembled on the stack, 
/// the processor will either perform a load into a vector register or simply use the value of the
/// stack everytime it is requested. This method is slow. Shuffles/unpacks, while complicated to 
/// write, are significantly faster, and should be used whenever possible. A similar method is used
/// for particle loading.
#define INTERP_FIELD_2D_SIMD(inner_dir, outer_dir, F_ENUM, indx2, indx3, outer_coeff, inner_coeff, field) { \
									\
    pvReal field##tmp1, field##tmp2, field##tmp3;			\
    pvReal field##shuffAC_1, field##shuffBD_1, field##shuffAC_2, field##shuffBD_2, field##shuffAC_3, field##shuffBD_3; \
    pvReal field##_in_A1, field##_in_B1, field##_in_C1, field##_in_D1;	\
    pvReal field##_in_A2, field##_in_B2, field##_in_C2, field##_in_D2;	\
    pvReal field##_in_A3, field##_in_B3, field##_in_C3, field##_in_D3;	\
    int i2_0 = indx2.v[0],						\
      i2_1 = indx2.v[1],						\
      i2_2 = indx2.v[2],						\
      i2_3 = indx2.v[3],						\
      i3_0 = indx3.v[0],						\
      i3_1 = indx3.v[1],						\
      i3_2 = indx3.v[2],						\
      i3_3 = indx3.v[3];						\
									\
    field##_in_A1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 - 1,i3_0 - 1)); \
    field##_in_B1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 - 1,i3_1 - 1)); \
    field##_in_C1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2 - 1,i3_2 - 1)); \
    field##_in_D1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3 - 1,i3_3 - 1)); \
    field##shuffAC_1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_C1.r); \
    field##shuffBD_1.r = pv_unpacklo_real( field##_in_B1.r, field##_in_D1.r); \
    field##tmp1.r = pv_unpacklo_real( field##shuffAC_1.r, field##shuffBD_1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##m##inner_dir.r, field##tmp1.r);	\
									\
    field##_in_A2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0,i3_0 - 1)); \
    field##_in_B2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1,i3_1 - 1)); \
    field##_in_C2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2,i3_2 - 1)); \
    field##_in_D2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3,i3_3 - 1)); \
    field##shuffAC_2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_C2.r); \
    field##shuffBD_2.r = pv_unpacklo_real( field##_in_B2.r, field##_in_D2.r); \
    field##tmp2.r = pv_unpacklo_real( field##shuffAC_2.r, field##shuffBD_2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##O##inner_dir.r, field##tmp2.r);	\
									\
    field##_in_A3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 + 1,i3_0 - 1)); \
    field##_in_B3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 + 1,i3_1 - 1)); \
    field##_in_C3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2 + 1,i3_2 - 1)); \
    field##_in_D3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3 + 1,i3_3 - 1)); \
    field##shuffAC_3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_C3.r); \
    field##shuffBD_3.r = pv_unpacklo_real( field##_in_B3.r, field##_in_D3.r); \
    field##tmp3.r = pv_unpacklo_real( field##shuffAC_3.r, field##shuffBD_3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##l##inner_dir.r, field##tmp3.r);	\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field.r = pv_mul_real(outer_coeff##m##outer_dir.r, field##tmp1.r);		\
    									\
    field##_in_A1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 - 1,i3_0)); \
    field##_in_B1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 - 1,i3_1)); \
    field##_in_C1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2 - 1,i3_2)); \
    field##_in_D1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3 - 1,i3_3)); \
    field##shuffAC_1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_C1.r); \
    field##shuffBD_1.r = pv_unpacklo_real( field##_in_B1.r, field##_in_D1.r); \
    field##tmp1.r = pv_unpacklo_real( field##shuffAC_1.r, field##shuffBD_1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##m##inner_dir.r, field##tmp1.r);		\
									\
    field##_in_A2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0,i3_0)); \
    field##_in_B2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1,i3_1)); \
    field##_in_C2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2,i3_2)); \
    field##_in_D2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3,i3_3)); \
    field##shuffAC_2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_C2.r); \
    field##shuffBD_2.r = pv_unpacklo_real( field##_in_B2.r, field##_in_D2.r); \
    field##tmp2.r = pv_unpacklo_real( field##shuffAC_2.r, field##shuffBD_2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##O##inner_dir.r, field##tmp2.r);		\
									\
    field##_in_A3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 + 1,i3_0)); \
    field##_in_B3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 + 1,i3_1)); \
    field##_in_C3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2 + 1,i3_2)); \
    field##_in_D3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3 + 1,i3_3)); \
    field##shuffAC_3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_C3.r); \
    field##shuffBD_3.r = pv_unpacklo_real( field##_in_B3.r, field##_in_D3.r); \
    field##tmp3.r = pv_unpacklo_real( field##shuffAC_3.r, field##shuffBD_3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##l##inner_dir.r, field##tmp3.r);	\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##O##outer_dir.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
    									\
    field##_in_A1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 - 1,i3_0 +1)); \
    field##_in_B1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 - 1,i3_1 +1)); \
    field##_in_C1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2 - 1,i3_2 +1)); \
    field##_in_D1.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3 - 1,i3_3 +1)); \
    field##shuffAC_1.r = pv_unpacklo_real( field##_in_A1.r, field##_in_C1.r); \
    field##shuffBD_1.r = pv_unpacklo_real( field##_in_B1.r, field##_in_D1.r); \
    field##tmp1.r = pv_unpacklo_real( field##shuffAC_1.r, field##shuffBD_1.r); \
    field##tmp1.r = pv_mul_real(inner_coeff##m##inner_dir.r, field##tmp1.r);		\
    									\
    field##_in_A2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0,i3_0 +1)); \
    field##_in_B2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1,i3_1 +1)); \
    field##_in_C2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2,i3_2 +1)); \
    field##_in_D2.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3,i3_3 +1)); \
    field##shuffAC_2.r = pv_unpacklo_real( field##_in_A2.r, field##_in_C2.r); \
    field##shuffBD_2.r = pv_unpacklo_real( field##_in_B2.r, field##_in_D2.r); \
    field##tmp2.r = pv_unpacklo_real( field##shuffAC_2.r, field##shuffBD_2.r); \
    field##tmp2.r = pv_mul_real(inner_coeff##O##inner_dir.r, field##tmp2.r);	\
    									\
    field##_in_A3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_0 +1,i3_0 +1)); \
    field##_in_B3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_1 +1,i3_1 +1)); \
    field##_in_C3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_2 +1,i3_2 +1)); \
    field##_in_D3.r = pv_loads_real(F_ENUM##point +			\
				    FIELD_OFF_##inner_dir##outer_dir(i2_3 +1,i3_3 +1)); \
    field##shuffAC_3.r = pv_unpacklo_real( field##_in_A3.r, field##_in_C3.r); \
    field##shuffBD_3.r = pv_unpacklo_real( field##_in_B3.r, field##_in_D3.r); \
    field##tmp3.r = pv_unpacklo_real( field##shuffAC_3.r, field##shuffBD_3.r); \
    field##tmp3.r = pv_mul_real(inner_coeff##l##inner_dir.r, field##tmp3.r);		\
    									\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp2.r);		\
    field##tmp1.r = pv_add_real(field##tmp1.r, field##tmp3.r);		\
    field##tmp1.r = pv_mul_real(outer_coeff##l##outer_dir.r, field##tmp1.r);	\
    field.r = pv_add_real(field.r, field##tmp1.r);			\
  }

#endif

#define INTERP_FIELD_YZ(m, jl2, jl3, gh3, gh2, var) INTERP_FIELD_2D_SIMD(y, z, m, jl2, jl3, gh3, gh2, var)
//#define INTERP_FIELD_YZ INTERP_FIELD_YZ_SLOW
#define INTERP_FIELD_XZ(m, jl1, jl3, gh3, gh1, var) INTERP_FIELD_2D_SIMD(x, z, m, jl1, jl3, gh3, gh1, var)
//#define INTERP_FIELD_XZ INTERP_FIELD_XZ_SLOW
#endif

// Automatically added by psc_doxyfy.sh to conform to doxygen standards
// It would be a good idea to fill this in with a real description

/// \file simd_wrap.h Wrappers for SSE2 functions to allow precision switching.
