#ifndef PSC_SPU_H
#define PSC_SPU_H 
// #include <simdmath.h>
#include "../psc_cbe.h"





// Some global variables on the spu
extern psc_cell_ctx_t spu_ctx;
extern psc_cell_block_t psc_block;


/// \FIXME Ordering of particle load/store operations is 
/// based on the estimated latencies of 
/// the comparable sse2 instructions. It may not be the best layout for 
/// the spe. Likewise this method may be very slow

#if 0

#if CBE_DOUBLE
#define spu_sqrt( vec ) sqrtd2( vec )
#define spu_div( vec1, vec2 ) divd2( (vec1), (vec2) ) /// \FIXME: this could be slow
#define spu_round_real( vec ) roundd2( (vec) )
#define spu_round_int( vec ) llroundd2( (vec) )

typedef vector double v_real;
typedef vector signed long long v_int;

const vector unsigned char uplo_pat = 
  {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, // first word from first vec
   0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17}; // first word from second vec

const vector unsigned char uphi_pat = 
  {0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, // second word from first vec
   0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f}; // second word from second vec

const vector signed long long element_assign[2] = {{-1ll, 0}, {0, -1ll}};

#define LOAD_PARTICLES_SPU(n) {						\
    v_real __A, __B, __C, __D, __E, __F, __G, __H;			\
    __A = *((v_real * ) &(part[(n)].xi));				\
    __E = *((v_real * ) &(part[(n)+1].xi));				\
									\
    __B = *((v_real * ) &(part[(n)].zi));				\
    __F = *((v_real * ) &(part[(n)+1].zi));				\
									\
    xi = spu_shuffle(__A, __E, uplo_pat);				\
    yi = spu_shuffle(__A, __E, uphi_pat);				\
									\
    zi = spu_shuffle(__B, __F, uplo_pat);				\
    pxi = spu_shuffle(__B, __F, uphi_pat);				\
									\
    __C = *((v_real * ) &(part[(n)].pyi));				\
    __G = *((v_real * ) &(part[(n)+1].pyi));				\
									\
    __D = *((v_real * ) &(part[(n)].qni));				\
    __H = *((v_real * ) &(part[(n)+1].qni));				\
									\
    pyi = spu_shuffle(__C, __G, uplo_pat);				\
    pzi = spu_shuffle(__C, __G, uphi_pat);				\
									\
    qni = spu_shuffle(__D, __H, uplo_pat);				\
    mni = spu_shuffle(__D, __H, uphi_pat);				\
									\
    __A = *((v_real * ) &(part[(n)].wni));				\
    __E = *((v_real * ) &(part[(n)+1].wni));				\
									\
    wni = spu_shuffle(__A, __E, uplo_pat);				\
  }

#define STORE_PARTICLES_SPU(n) {		\
    v_real __A, __B, __C, __D, __E, __F;	\
						\
						\
    __A = spu_shuffle(xi, yi, uplo_pat);	\
						\
    *((v_real * ) &(part[(n)].xi)) = __A;	\
						\
    __B = spu_shuffle(zi, pxi, uplo_pat);	\
    __C = spu_shuffle(pyi, pzi, uplo_pat);	\
						\
    *((v_real * ) &(part[(n)].zi)) = __B;	\
						\
    __E = spu_shuffle(xi, yi, uphi_pat);	\
    __F = spu_shuffle(zi, pxi, uphi_pat);	\
						\
    *((v_real * ) &(part[(n)].pyi)) = __C ;	\
						\
    __D = spu_shuffle(pyi, pzi, uphi_pat);	\
						\
    *((v_real * ) &(part[(n)+1].xi)) = __E;	\
    *((v_real * ) &(part[(n)+1].zi)) = __F;	\
    *((v_real * ) &(part[(n)+1].pyi)) = __D;	\
  }


#else 
#define spu_sqrt( vec ) sqrtf4( vec )
#define spu_div( vec1, vec2 ) divf4( (vec1), (vec2) ) 
#define spu_round_real( vec ) roundf4( (vec) )
#define spu_round_int( vec ) llroundf4( (vec))

typedef vector float v_real;
typedef llroundf4_t v_int;

const vector unsigned char uplo_pat = 
    {0x00, 0x01, 0x02, 0x03,  // first word from first vec
     0x10, 0x11, 0x12, 0x13, // first word from second vec
     0x04, 0x05, 0x06, 0x07, // second word from first vec
     0x14, 0x15, 0x16, 0x17}; // second word from second vec
 
const vector unsigned char uphi_pat = 
    {0x08, 0x09, 0x0a, 0x0b, // third word from first vec 
     0x18, 0x19, 0x1a, 0x1b,  // third word from second vec
     0x0c, 0x0d, 0x0e, 0x0f, // fourth word from first vec
     0x1c, 0x1d, 0x1e, 0x1f}; // fourth word from second vec

#define LOAD_PARTICLES_SPU(n) {			\
  v_real __A, __B, __C, __D, __E, __F, __G, __H, __alpha, __beta, __gamma, __delta; \
  __A = *((v_real * ) &(part[(n)].xi));					\
  __B = *((v_real * ) &(part[(n)+1].xi));				\
  __C = *((v_real * ) &(part[(n)+2].xi));				\
  __D = *((v_real * ) &(part[(n)+3].xi));				\
									\
  __alpha = spu_shuffle(__A, __C, uplo_pat);				\
  __gamma = spu_shuffle(__A, __C, uphi_pat);				\
  __beta = spu_shuffle(__B, __D, uplo_pat);				\
  __delta = spu_shuffle(__B, __D, uphi_pat);				\
									\
  __E = *((v_real *) &(part[(n)].pyi));					\
  __G = *((v_real *) &(part[(n)+2].pyi));				\
									\
  xi = spu_shuffle(__alpha, __beta, uplo_pat);				\
  yi = spu_shuffle(__alpha, __beta, uphi_pat);				\
  zi = spu_shuffle(__gamma, __delta, uplo_pat);				\
  pxi = spu_shuffle(__gamma, __delta, uphi_pat);			\
									\
  __F = *((v_real *) &(part[(n)+1].pyi));				\
  __H = *((v_real *) &(part[(n)+3].pyi));				\
									\
  __gamma = spu_shuffle(__E, __G, uphi_pat);				\
  __alpha = spu_shuffle(__E, __G, uplo_pat);				\
  __beta = spu_shuffle(__F, __H, uplo_pat);				\
  __delta = spu_shuffle(__F, __H, uphi_pat);				\
									\
  pyi = spu_shuffle(__alpha, __beta, uplo_pat);				\
  pzi = spu_shuffle(__alpha, __beta, uphi_pat);				\
  qni = spu_shuffle(__gamma, __delta, uplo_pat);			\
  mni = spu_shuffle(__gamma, __delta, uphi_pat);			\
									\
  __A = *((v_real *) &(part[(n)].wni));					\
  __C = *((v_real *) &(part[(n)+2].wni));				\
  __B = *((v_real *) &(part[(n)+1].wni));				\
  __D = *((v_real *) &(part[(n)+3].wni));				\
									\
  __alpha = spu_shuffle(__A, __C, uplo_pat);				\
  __beta = spu_shuffle(__B, __D, uphi_pat);				\
  wni = spu_shuffle(__alpha, __beta, uplo_pat);				\
  }

#define STORE_PARTICLES_SPU(n) {			\
    v_real __A, __B, __C, __D;				\
    v_real __alpha, __beta, __gamma, __delta;		\
  __alpha = spu_shuffle(xi, zi, uplo_pat);		\
  __beta = spu_shuffle(yi, pxi, uplo_pat);		\
  __A = spu_shuffle(__alpha, __beta, uplo_pat);		\
  *((v_real *) &(part[(n)].xi)) = __A;			\
  __gamma = spu_shuffle(pyi, qni, uplo_pat);		\
  __delta = spu_shuffle(pzi, mni, uplo_pat);		\
  __C = spu_shuffle(__gamma, __delta, uplo_pat);	\
  *((v_real *) &(part[(n)].pyi)) = __C;			\
  __B = spu_shuffle(__alpha, __beta, uphi_pat);		\
  *((v_real *) &(part[(n)+1].xi)) = __B;		\
  __D = spu_shuffle(__gamma, __delta, uphi_pat);	\
  *((v_real *) &(part[(n)+1].pyi)) = __D;		\
  __alpha = spu_shuffle(xi, zi, uphi_pat);		\
  __beta = spu_shuffle(yi, pxi, uphi_pat);		\
  __A = spu_shuffle(__alpha, __beta, uplo_pat);		\
  *((v_real *) &(part[(n)+2].xi)) = __A;		\
  __gamma = spu_shuffle(pyi, qni, uphi_pat);		\
b   __delta = spu_shuffle(pzi, mni, uphi_pat);		\
  __C = spu_shuffle(__gamma, __delta, uplo_pat);	\
  *((v_real *) &(part[(n)+2].pyi)) = __C;		\

  __B = spu_shuffle(__alpha, __beta, uphi_pat);		\
  *((v_real *) &(part[(n)+3].xi)) = __B;		\
  __D = spu_shuffle(__gamma, __delta, uphi_pat);	\
  *((v_real *) &(part[(n)+3].pyi)) = __D;		\
  }

  
#endif

#define F3_OFF_SPU(fldnr, jx,jy,jz)					\
  ((((((jz)-ilg[2]))						\
     *img[1] + ((jy)-ilg[1]))					\
    *img[0] + ((jx)-ilg[0]))					\
   *NR_FIELDS + fldnr)




#endif 
void spu_dma_get(volatile void *ls, unsigned long long ea, unsigned long size);
#endif
