#ifndef SPU_PARTICLES_H
#define SPU_PARTICLES_H

#include "psc_spu.h"
#include <psc_particles_cbe.h>

#define LOAD_PARTICLES_SPU {						\
    particle_cbe_t *part1 = buff.lb1;					\
    particle_cbe_t *part2 = buff.lb2;					\
    v_real __A, __B, __C, __D, __E, __F, __G, __H;			\
    __A = *((v_real * ) &(part1->xi));				\
    __E = *((v_real * ) &(part2->xi));				\
									\
    __B = *((v_real * ) &(part1->zi));				\
    __F = *((v_real * ) &(part2->zi));				\
									\
    xi = spu_shuffle(__A, __E, uplo_pat);				\
    yi = spu_shuffle(__A, __E, uphi_pat);				\
									\
    zi = spu_shuffle(__B, __F, uplo_pat);				\
    pxi = spu_shuffle(__B, __F, uphi_pat);				\
									\
    __C = *((v_real * ) &(part1->pyi));				\
    __G = *((v_real * ) &(part2->pyi));				\
									\
    __D = *((v_real * ) &(part1->qni));				\
    __H = *((v_real * ) &(part2->qni));				\
									\
    pyi = spu_shuffle(__C, __G, uplo_pat);				\
    pzi = spu_shuffle(__C, __G, uphi_pat);				\
									\
    qni = spu_shuffle(__D, __H, uplo_pat);				\
    mni = spu_shuffle(__D, __H, uphi_pat);				\
									\
    __A = *((v_real * ) &(part1->wni));				\
    __E = *((v_real * ) &(part2->wni));				\
									\
    wni = spu_shuffle(__A, __E, uplo_pat);				\
  }

#define STORE_PARTICLES_SPU {		\
  v_real __A, __B, __C, __D, __E, __F;					\
  particle_cbe_t *part1 = buff.lb1;					\
  particle_cbe_t *part2 = buff.lb2;					\
  									\
						\
    __A = spu_shuffle(xi, yi, uplo_pat);	\
						\
    *((v_real * ) &(part1->xi)) = __A;	\
						\
    __B = spu_shuffle(zi, pxi, uplo_pat);	\
    __C = spu_shuffle(pyi, pzi, uplo_pat);	\
						\
    *((v_real * ) &(part1->zi)) = __B;	\
						\
    __E = spu_shuffle(xi, yi, uphi_pat);	\
    __F = spu_shuffle(zi, pxi, uphi_pat);	\
						\
    *((v_real * ) &(part1->pyi)) = __C ;	\
						\
    __D = spu_shuffle(pyi, pzi, uphi_pat);	\
						\
    *((v_real * ) &(part2->xi)) = __E;	\
    *((v_real * ) &(part2->zi)) = __F;	\
    *((v_real * ) &(part2->pyi)) = __D;	\
  }


struct part_buffs {
  particle_cbe_t *plb1, *plb2; // LS add of pre-load buffers
  particle_cbe_t *lb1, *lb2; // LS add of load buffers
  particle_cbe_t *sb1, *sb2; // LS add of store buffers 
};

extern struct part_buffs buff; 

extern const unsigned int tag_pget;
extern const unsigned int tag_pput;

#endif
