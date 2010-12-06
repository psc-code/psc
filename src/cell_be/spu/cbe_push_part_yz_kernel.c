#include <stdlib.h>
#include <simdmath.h>
#include <spu_intrinsics.h>
#include <assert.h>
#include <string.h>
#include "psc_spu.h"
#include "psc_alf.h"
#include "psc_particles_cbe.h"


#define print_vec( vec ) printf("%g %g \n", spu_extract((vec), 0), spu_extract((vec), 1))
#define print_ivec( vec ) printf("%lld %lld \n", spu_extract((vec), 0), spu_extract((vec), 1))

enum { 
  s_EX, s_EZ, s_HY, s_NRFLDS,
};

enum {
  s_JX, s_JZ, s_NRCURR,
};


v_real one, half, threefourths, onepfive;

static inline void
find_index(v_real *xi, v_real *dxi, v_int *j, v_real *h)
{
  v_real tmp;
  tmp = spu_mul(*xi, *dxi);
  *h = spu_round_real(tmp);
  *j = spu_round_int(tmp);
  *h = spu_sub(*h, tmp);
} 


static inline void
find_index_minus_shift(v_real *xi, v_real *dxi, v_int *j, v_real *h, v_real *shift)
{
  v_real tmp;
  tmp = spu_msub(*xi, *dxi, *shift);
  *h = spu_round_real(tmp);
  *j = spu_round_int(tmp);
  *h = spu_sub(*h, tmp);
} 

static inline void
ip_to_grid_m(v_real *h, v_real * xmx)
{
  *xmx = spu_add(half, *h);
  *xmx = spu_mul(*xmx, *xmx);
  *xmx = spu_mul(half, *xmx);
}

static inline void
ip_to_grid_O(v_real *h, v_real *xOx)
{
  *xOx = spu_mul( *h , *h);
  *xOx = spu_sub(threefourths, *xOx);

}

static inline void
ip_to_grid_l(v_real *h, v_real * xlx)
{
  *xlx = spu_sub(half, *h);
  *xlx = spu_mul( *xlx, *xlx);
  *xlx = spu_mul(half, *xlx);
}

static inline void
form_factor_m(v_real *h, v_real * restrict xmx)
{
  *xmx = spu_sub(*h,one);
  *xmx = spu_add(onepfive,*xmx); // h-1 always <0
  *xmx = spu_mul(*xmx, *xmx);
  *xmx = spu_mul(half, *xmx);
}

static inline void
form_factor_O(v_real *h, v_real * restrict xOx)
{
  *xOx = spu_mul(*h, *h);
  *xOx = spu_sub(threefourths, *xOx);
}

static inline void
form_factor_l(v_real *h, v_real * restrict xlx)
{
  *xlx = spu_add(one, *h); 
  *xlx = spu_sub(onepfive, *xlx);//h+1 always >0
  *xlx = spu_mul(*xlx, *xlx);
  *xlx = spu_mul(half, *xlx);
}


int push_yz_prep_in_dtl(void *p_task_context,
			void *p_parm_ctx_buffer,
			void *p_dtl,
			unsigned int current_count,
			unsigned int total_count __attribute__ ((unused)))
{
  //  printf("Bringing in\n");
  push_wb_context_t *p_parm = (push_wb_context_t*)p_parm_ctx_buffer;
  unsigned long long p_start = p_parm->p_start; 
  unsigned long long p_pad = p_parm->p_pad;
  unsigned int lo_npart = p_parm->lo_npart;
  unsigned int maxpart = p_parm->maxpart;

  push_task_context_t *p_ctx = (push_task_context_t*)p_task_context;
  int *wb_lg = p_parm->wb_lg;
  int *wb_hg = p_parm->wb_hg;
  unsigned long long p_fields = p_ctx->p_fields;
  int *ilg = p_ctx->ilg;
  int *img = p_ctx->img;
  ALF_ACCEL_DTL_BEGIN(p_dtl, ALF_BUF_IN, 0);
  for( int jz = wb_lg[2]; jz <= wb_hg[2]; jz++){
    for( int jy = wb_lg[1]; jy <= wb_hg[1]; jy++){
      ALF_ACCEL_DTL_ENTRY_ADD(p_dtl, VEC_SIZE*s_NRFLDS*sizeof(fields_cbe_real_t), 
			      ALF_DATA_BYTE, 
			      (p_fields + F3_OFF_SPU(EX,wb_lg[0],jy,jz)*sizeof(fields_cbe_real_t)));
    }
  }
  
  ALF_ACCEL_DTL_END(p_dtl);
  
  // Skip particles we already processed
  p_start += (current_count * maxpart) * sizeof(particle_cbe_t ); 

  unsigned int size;
  if((current_count + 1) * maxpart > lo_npart){
    size = lo_npart - current_count * maxpart;
  } else {
    size = maxpart;
  }


  ALF_ACCEL_DTL_BEGIN(p_dtl, ALF_BUF_OVL_INOUT, 0);

  // The Cell has a maximum limit of 16KB per DTL entry. Host 
  // partitioning automatically breaks a too big entry into smaller
  // chunks. Accelerator side doesn't, so we have to do it ourselves. 
  // Credit where credit is due, this method of breaking things into
  // chunks is very similar to that used in the dot_prod_multi_spu.c 
  // IBM alf example, but it's such a straight forward approach
  // that I don't feel bad using it.

  unsigned int max_chunk = 16000 / ((int) sizeof(particle_cbe_t ));

  for(int k = 0; k<size; k+=max_chunk){
    unsigned int chunk;
    if(k+max_chunk <= size)
      chunk = max_chunk;
    else 
      chunk = size - k;
    ALF_ACCEL_DTL_ENTRY_ADD(p_dtl, chunk*sizeof(particle_cbe_t ), ALF_DATA_BYTE, 
			    p_start + k*sizeof(particle_cbe_t ));
  }
  // Add the padding onto the end of the buffer to avoid reading bad memory
  // when grabbing the last wni
  ALF_ACCEL_DTL_ENTRY_ADD(p_dtl, sizeof(particle_cbe_t ), ALF_DATA_BYTE, 
  			  p_pad);
  p_pad += sizeof(particle_cbe_t );
  int num_pad; 
  if( size < maxpart ){
    if( (maxpart - size) < (VEC_SIZE - 1)){
      num_pad = maxpart-size;
    } else {
      num_pad = (VEC_SIZE - 1);
    }
    
    ALF_ACCEL_DTL_ENTRY_ADD(p_dtl, sizeof(particle_cbe_t ), ALF_DATA_BYTE,
			  p_pad);
  }
  ALF_ACCEL_DTL_END(p_dtl);
  
  return 0;
}

int
push_yz_prep_out_dtl(void *p_task_context __attribute__ ((unused)),
			void *p_parm_ctx_buffer,
			void *p_dtl,
			unsigned int current_count,
			unsigned int total_count)
{
  if(current_count == (total_count - 1)){
    //    printf("Dumping out\n");
    push_wb_context_t *p_parm = (push_wb_context_t*)p_parm_ctx_buffer;
    unsigned long long p_cache = p_parm->p_cache;
    int fld_size = (p_parm->wb_hg[0] - p_parm->wb_lg[0] + 1)
      * (p_parm->wb_hg[1] - p_parm->wb_lg[1] + 1)
      * (p_parm->wb_hg[2] - p_parm->wb_lg[2] + 1);

    // FIXME It is possible (though not likely) that the current data might
    // exceed the 16KB single dtl entry maximum. To be safe, a chunk system 
    // similar to the particles should be put in here.

    ALF_ACCEL_DTL_BEGIN(p_dtl, ALF_BUF_OUT, 0);
    ALF_ACCEL_DTL_ENTRY_ADD(p_dtl, 4*fld_size,
			    PSC_ALF_DATATYPE, p_cache);
    ALF_ACCEL_DTL_END(p_dtl);
  }
  return 0;
}

#define F2_SPU_OFF(jy,jz)			\
  ((((jz) - wb_lg[2]))				\
   *img[1] + ((jy) - wb_lg[1]))*s_NRFLDS

#define JSX(indx2, indx3) s_cur[(((indx3) - wb_lg[2])*img[1] + ((indx2) - wb_lg[1]))*s_NRCURR + s_JX]
#define JSZ(indx2, indx3) s_cur[(((indx3) - wb_lg[2])*img[1] + ((indx2) - wb_lg[1]))*s_NRCURR + s_JZ]


#define IP_FIELD_SPU(f_offset, word, indx2, indx3, outer_coeff, inner_coeff, field) { \
    v_real field##tmp1, field##tmp2, field##tmp3;			\
    v_real field##_in_A1, field##_in_B1;				\
    v_real field##_in_A2, field##_in_B2;				\
    v_real field##_in_A3, field##_in_B3;				\
    int i2_0 = (int)spu_extract(indx2,0);				\
    int i2_1 = (int)spu_extract(indx2,1);				\
    int i3_0 = (int)spu_extract(indx3,0);				\
    int i3_1 = (int)spu_extract(indx3,1);				\
									\
    field##_in_A1 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0 - 1,i3_0 - 1));			\
    field##_in_B1 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1 - 1,i3_1 - 1));			\
    field##tmp1 = spu_shuffle( field##_in_A1, field##_in_B1, up##word##_pat); \
    field##tmp1 = spu_mul(inner_coeff##my, field##tmp1);		\
    									\
    field##_in_A2 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0,i3_0 - 1));			\
    field##_in_B2 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1,i3_1 - 1));			\
    field##tmp2 = spu_shuffle( field##_in_A2, field##_in_B2, up##word##_pat); \
    field##tmp2 = spu_mul(inner_coeff##Oy, field##tmp2);		\
    									\
    field##_in_A3 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0 + 1,i3_0 - 1));			\
    field##_in_B3 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1 + 1,i3_1 - 1));			\
    field##tmp3 = spu_shuffle( field##_in_A3, field##_in_B3, up##word##_pat); \
    field##tmp3 = spu_mul(inner_coeff##ly, field##tmp3);		\
    									\
    field##tmp1 = spu_add(field##tmp1, field##tmp2);			\
    field##tmp1 = spu_add(field##tmp1, field##tmp3);			\
    field = spu_mul(outer_coeff##mz, field##tmp1);			\
    									\
    field##_in_A1 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0 - 1,i3_0));			\
    field##_in_B1 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1 - 1,i3_1));			\
    field##tmp1 = spu_shuffle( field##_in_A1, field##_in_B1, up##word##_pat); \
    field##tmp1 = spu_mul(inner_coeff##my, field##tmp1);		\
    									\
    field##_in_A2 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0,i3_0));				\
    field##_in_B2 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1,i3_1));				\
    field##tmp2 = spu_shuffle( field##_in_A2, field##_in_B2, up##word##_pat); \
    field##tmp2 = spu_mul(inner_coeff##Oy, field##tmp2);		\
    									\
    field##_in_A3 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0 + 1,i3_0));			\
    field##_in_B3 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1 + 1,i3_1));			\
    field##tmp3 = spu_shuffle( field##_in_A3, field##_in_B3, up##word##_pat); \
    field##tmp3 = spu_mul(inner_coeff##ly, field##tmp3);		\
    									\
    field##tmp1 = spu_add(field##tmp1, field##tmp2);			\
    field##tmp1 = spu_add(field##tmp1, field##tmp3);			\
    field##tmp1 = spu_mul(outer_coeff##Oz, field##tmp1);		\
    field = spu_add(field, field##tmp1);				\
    									\
    field##_in_A1 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0 - 1,i3_0 +1));			\
    field##_in_B1 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1 - 1,i3_1 +1));			\
    field##tmp1 = spu_shuffle( field##_in_A1, field##_in_B1, up##word##_pat); \
    field##tmp1 = spu_mul(inner_coeff##my, field##tmp1);		\
    									\
    field##_in_A2 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0,i3_0 +1));			\
    field##_in_B2 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1,i3_1 +1));			\
    field##tmp2 = spu_shuffle( field##_in_A2, field##_in_B2, up##word##_pat); \
    field##tmp2 = spu_mul(inner_coeff##Oy, field##tmp2);		\
    									\
    field##_in_A3 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_0 +1,i3_0 +1));			\
    field##_in_B3 = *(flds + f_offset +					\
		      F2_SPU_OFF(i2_1 +1,i3_1 +1));			\
    field##tmp3 = spu_shuffle( field##_in_A3, field##_in_B3, up##word##_pat); \
    field##tmp3 = spu_mul(inner_coeff##ly, field##tmp3);		\
    									\
    field##tmp1 = spu_add(field##tmp1, field##tmp2);			\
    field##tmp1 = spu_add(field##tmp1, field##tmp3);			\
    field##tmp1 = spu_mul(outer_coeff##lz, field##tmp1);		\
    field = spu_add(field, field##tmp1);				\
  }


int cbe_push_part_yz_kernel(void *p_task_context, // const accross all pars here
			    void *p_parm_context, // WB specific: domain, num part, etc
			    void *p_input_buffer, // pointer to fields comes in here
			    void *p_output_buffer, // will be currents
			    void *p_inout_buffer, // particles here
			    unsigned int current_count,
			    unsigned int total_count )
{
  // Happy initialization here
  push_task_context_t *p_ctx = (push_task_context_t *) p_task_context;
  push_wb_context_t *p_parm = (push_wb_context_t *) p_parm_context;
  
  v_real * flds = (v_real *) p_input_buffer; 
  v_real * s_cur = (v_real *)(p_parm_context + sizeof(push_wb_context_t)+ (16-sizeof(push_wb_context_t)%16));

  int *wb_lg = p_parm->wb_lg;
  int *wb_hg = p_parm->wb_hg;
  int img[3] = {1, wb_hg[1] - wb_lg[1] + 1, wb_hg[2] - wb_lg[2] + 1};
  v_real yl = spu_splats(p_ctx->yl);
  v_real zl = spu_splats(p_ctx->zl);
  v_real dyi = spu_splats(p_ctx->dyi);
  v_real dzi = spu_splats(p_ctx->dzi);
  v_real dqs = spu_splats(p_ctx->dqs);
  v_real fnqs = spu_splats(p_ctx->fnqs);
  v_real fnqxs = spu_splats(p_ctx->fnqxs);
  v_real fnqys = spu_splats(p_ctx->fnqys);
  v_real fnqzs = spu_splats(p_ctx->fnqzs);
  
  int fld_size = (p_parm->wb_hg[0] - p_parm->wb_lg[0] + 1)
    * (p_parm->wb_hg[1] - p_parm->wb_lg[1] + 1)
    * (p_parm->wb_hg[2] - p_parm->wb_lg[2] + 1);
  
  if(current_count == 0){
    memset((void *)s_cur, 0, 4*fld_size*sizeof(fields_cbe_real_t));
  }
  
  particle_cbe_t  *part = (particle_cbe_t  *) p_inout_buffer;
  
  // Hacky to get it to work. Will fix later.
#if CBE_DOUBLE
  one = spu_splats(1.0);
  half = spu_splats(0.5);
  threefourths = spu_splats(0.75);
  onepfive = spu_splats(1.5);
  v_real third = spu_splats(1./3.);
  v_real zero = spu_splats(0.0);
#else 
  one = spu_splats(1.0f);
  half = spu_splats(0.5f);
  onepfive = spu_splats(1.5f);
  threefourths = spu_splats(0.75f);
#endif

  /// \FIXME There is actually a faster way to do this using immediate loads of a single value.
  /// I have left it this slower way because it's only done once per WB (so the performance hit
  /// isn't too bad) and it's a bit more transparent what's happening, which I need at the moment.
  
  unsigned int lo_npart = p_parm->lo_npart;
  unsigned int maxpart = p_parm->maxpart;
  unsigned int size;
  if((current_count + 1) * maxpart > lo_npart){
    size = lo_npart - current_count * maxpart;
  } else {
    size = maxpart;
  }
  

  for(int n = 0; n < size; n += VEC_SIZE)
    {

      v_real xi, yi, zi, pxi, pyi, pzi, qni, mni, wni;
	
      LOAD_PARTICLES_SPU(n);

      v_real vxi,vyi,vzi,root,tmpx,tmpy,tmpz; 
      
      tmpx = spu_mul(pxi,pxi);
      tmpy = spu_mul(pyi,pyi);
      tmpz = spu_mul(pzi, pzi);
    
      tmpx = spu_add(tmpx, tmpy);
      tmpz = spu_add(one, tmpz);
      tmpx = spu_add(tmpx, tmpz);
      root = spu_sqrt(tmpx);
      
      vxi = spu_div(pxi, root);
      vyi = spu_div(pyi, root);
      vzi = spu_div(pzi, root);
      
      tmpy = spu_mul(vyi, yl);
      tmpz = spu_mul(vzi, zl);
      
      yi = spu_add(yi, tmpy);
      zi = spu_add(zi, tmpz);

      v_real gmy, gmz, gOy, gOz, gly, glz, H2, H3, h2, h3;
      v_int j2, j3, l2, l3;
      
      find_index(&yi, &dyi, &j2, &H2);
      find_index(&zi, &dzi, &j3, &H3);
      
      ip_to_grid_m(&H2, &gmy);
      ip_to_grid_m(&H3, &gmz);
      
      ip_to_grid_O(&H2, &gOy);
      ip_to_grid_O(&H3, &gOz);

      ip_to_grid_l(&H2, &gly);
      ip_to_grid_l(&H3, &glz);

      find_index_minus_shift(&yi, &dyi, &l2, &h2, &half);
      find_index_minus_shift(&zi, &dzi, &l3, &h3, &half);

      v_real hmy, hmz, hOy, hOz, hly, hlz;

      ip_to_grid_m(&h2, &hmy);
      ip_to_grid_m(&h3, &hmz);
      
      ip_to_grid_O(&h2, &hOy);
      ip_to_grid_O(&h3, &hOz);

      ip_to_grid_l(&h2, &hly);
      ip_to_grid_l(&h3, &hlz);

      // Field interpolation here. urg...

      v_real exq, eyq, ezq, hxq, hyq, hzq;

      IP_FIELD_SPU(s_EX,lo,j2,j3,g,g,exq);
      IP_FIELD_SPU(s_EX,hi,l2,j3,g,h,eyq);
      IP_FIELD_SPU(s_EZ,lo,j2,l3,h,g,ezq);
      IP_FIELD_SPU(s_EZ,hi,l2,l3,h,h,hxq);
      IP_FIELD_SPU(s_HY,lo,j2,l3,h,g,hyq);
      IP_FIELD_SPU(s_HY,hi,l2,j3,g,h,hzq);

      v_real s0y[5], s0z[5], s1y[5], s1z[5];
      for(int mp=0; mp<5; mp++){
	s0y[mp] = spu_splats(0.0);
	s1y[mp] = spu_splats(0.0);
	s0z[mp] = spu_splats(0.0);
	s1z[mp] = spu_splats(0.0);
      }
      form_factor_m(&H2, &s0y[1]);
      form_factor_O(&H2, &s0y[2]);
      form_factor_l(&H2, &s0y[3]);
      
      form_factor_m(&H3, &s0z[1]);
      form_factor_O(&H3, &s0z[2]);
      form_factor_l(&H3, &s0z[3]);


      // It would be nice to offload the momentum to a function, but it may not
      // be necessary

      v_real dq;

      dq = spu_mul(qni, dqs);
      dq = spu_div(dq, mni);

      /*
      pxi = spu_madd(dq, exq, pxi);
      pyi = spu_madd(dq, eyq, pyi);
      pzi = spu_madd(dq, ezq, pzi);
      
      */
      v_real dqex, dqey, dqez;
      dqex = spu_mul(dq, exq);
      dqey = spu_mul(dq, eyq);
      dqez = spu_mul(dq, ezq);

      pxi = spu_add(pxi, dqex);
      pyi = spu_add(pyi, dqey);
      pzi = spu_add(pzi, dqez);
      
      v_real taux, tauy, tauz, txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;

      tmpx = spu_mul(pxi, pxi);
      tmpy = spu_mul(pyi, pyi);
      tmpz = spu_mul(pzi, pzi);
      root = spu_add(one, tmpx);
      tmpy = spu_add(tmpy, tmpz);
      root = spu_add(root, tmpy);
      root = spu_sqrt(root);
      root = spu_div(dq, root);

      taux = spu_mul(hxq, root);
      tauy = spu_mul(hyq, root);
      tauz = spu_mul(hzq, root);


      txx = spu_mul(taux, taux);
      tyy = spu_mul(tauy, tauy);
      tzz = spu_mul(tauz, tauz);
      t2xy = spu_mul(taux, tauy);
      t2xz = spu_mul(taux, tauz);
      t2yz = spu_mul(tauy, tauz);
      t2xy = spu_add(t2xy, t2xy);
      t2xz = spu_add(t2xz, t2xz);
      t2yz = spu_add(t2yz, t2yz);

      v_real tau;
      tau = spu_add(one, txx);
      tmpx = spu_add(tyy, tzz);
      tau = spu_add(tau, tmpx);
      tau = spu_div(one, tau);

      taux = spu_add(taux, taux);
      tauy = spu_add(tauy, tauy);
      tauz = spu_add(tauz, tauz);

      //pxp
      tmpx = spu_add(one, txx);
      tmpx = spu_sub(tmpx, tyy);
      tmpx = spu_sub(tmpx, tzz);
      tmpx = spu_mul(tmpx, pxi);

      tmpy = spu_add(t2xy, tauz);
      tmpy = spu_mul(tmpy, pyi);
      
      tmpz = spu_sub(t2xz, tauy);
      tmpz = spu_mul(tmpz, pzi);
      
      pxp = spu_add(tmpx, tmpy);
      pxp = spu_add(pxp, tmpz);
      pxp = spu_mul(pxp, tau);

      //pyp
      tmpx = spu_sub(t2xy, tauz);
      tmpx = spu_mul(tmpx, pxi);
      
      tmpy = spu_sub(one, txx);
      tmpy = spu_add(tmpy, tyy);
      tmpy = spu_sub(tmpy, tzz);
      tmpy = spu_mul(tmpy, pyi);

      tmpz = spu_add(t2yz, taux);
      tmpz = spu_mul(tmpz, pzi);

      pyp = spu_add(tmpx, tmpy);
      pyp = spu_add(pyp, tmpz);
      pyp = spu_mul(pyp, tau);

      // pzp
      tmpx = spu_add(t2xz, tauy);
      tmpx = spu_mul(tmpx, pxi);

      tmpy = spu_sub(t2yz, taux);
      tmpy = spu_mul(tmpy, pyi);

      tmpz = spu_sub(one, txx);
      tmpz = spu_sub(tmpz, tyy);
      tmpz = spu_add(tmpz, tzz);
      tmpz = spu_mul(tmpz, pzi);

      pzp = spu_add(tmpx, tmpy);
      pzp = spu_add(pzp, tmpz);
      pzp = spu_mul(pzp, tau);

      pxi = spu_add(pxp, dqex);
      pyi = spu_add(pyp, dqey);
      pzi = spu_add(pzp, dqez);

      /*
      pxi = spu_madd(dq, exq, pxp);
      pyi = spu_madd(dq, eyq, pyp);
      pzi = spu_madd(dq, ezq, pzp);
      */
      // finish advancing particles
      tmpx = spu_mul(pxi,pxi);
      tmpy = spu_mul(pyi,pyi);
      tmpz = spu_mul(pzi,pzi);
    
      tmpx = spu_add(tmpx, tmpy);
      tmpz = spu_add(one, tmpz);
      tmpx = spu_add(tmpx, tmpz);
      root = spu_sqrt(tmpx);
      
      vxi = spu_div(pxi, root);
      vyi = spu_div(pyi, root);
      vzi = spu_div(pzi, root);
      
      tmpy = spu_mul(vyi, yl);
      tmpz = spu_mul(vzi, zl);
      
      yi = spu_add(yi, tmpy);
      zi = spu_add(zi, tmpz);
     

      STORE_PARTICLES_SPU(n);


      yi = spu_add(yi, tmpy);
      zi = spu_add(zi, tmpz);

      v_int k2,k3;

      find_index(&yi, &dyi, &k2, &h2);
      find_index(&zi, &dzi, &k3, &h3);
      
      // God help me, there's some things I just can't fgure out how to parallelize
      // The g-- here are just temporary variables. I can't do the assignments in parallel, 
      // but I'll be damned if I can't do the FLOPS in parallel

      form_factor_m(&h2, &gmy);
      form_factor_O(&h2, &gOy);
      form_factor_l(&h2, &gly);
      
      form_factor_m(&h3, &gmz);
      form_factor_O(&h3, &gOz);
      form_factor_l(&h3, &glz);

      signed long long j2_scal[2], j3_scal[2];
      j2_scal[0] = spu_extract(j2, 0);
      j2_scal[1] = spu_extract(j2, 1);
      j3_scal[0] = spu_extract(j3, 0);
      j3_scal[1] = spu_extract(j3, 1);

      signed long long dfy[2];
      signed long long  dfz[2];
      dfy[0] = spu_extract(k2,0) - j2_scal[0];
      dfy[1] = spu_extract(k2,1) - j2_scal[1];
      dfz[0] = spu_extract(k3,0) - j3_scal[0];
      dfz[1] = spu_extract(k3,1) - j3_scal[1];
      
  
      for(int p=0; p < VEC_SIZE; p++){
	s1y[(int)dfy[p] + 1] = spu_sel(s1y[(int)dfy[p] + 1], gmy, (vector unsigned long long) element_assign[p]);
	s1y[(int)dfy[p] + 2] = spu_sel(s1y[(int)dfy[p] + 2], gOy, (vector unsigned long long) element_assign[p]);
	s1y[(int)dfy[p] + 3] = spu_sel(s1y[(int)dfy[p] + 3], gly, (vector unsigned long long) element_assign[p]);
	s1z[(int)dfz[p] + 1] = spu_sel(s1z[(int)dfz[p] + 1], gmz, (vector unsigned long long) element_assign[p]);
	s1z[(int)dfz[p] + 2] = spu_sel(s1z[(int)dfz[p] + 2], gOz, (vector unsigned long long) element_assign[p]);
	s1z[(int)dfz[p] + 3] = spu_sel(s1z[(int)dfz[p] + 3], glz, (vector unsigned long long) element_assign[p]);
      }

      for(int m=0; m<5; m++){
	s1y[m] = spu_sub(s1y[m], s0y[m]);
	s1z[m] = spu_sub(s1z[m], s0z[m]);
      }
      
      // This section doesn't actually
      // do anything. I'm leaving it here
      // because I want a record that I did
      // it, if only in a branch which will be
      // deprecated.
      /*      
      v_int lowb2, lowb3, hib2, hib3;
      v_int itmp; 
      v_int diff2 = spu_promote(dfy[0], 0);
      diff2 = spu_insert(dfy[1],diff2,1);
      v_int diff3 = spu_promote(dfz[0], 0);
      diff3 = spu_insert(dfz[1],diff3,1);
      
      itmp = (vector signed long long)spu_rl((vector signed int)diff2, 1);
      lowb2 = spu_and(diff2, itmp);
      hib2 = negatell2(diff2);
      itmp = (vector signed long long)spu_rl((vector signed int)hib2, 1);
      hib2 = spu_and(hib2,itmp);

      itmp = (vector signed long long)spu_rl((vector signed int)diff3, 1);
      lowb3 = spu_and(diff3, itmp);
      hib3 = negatell2(diff3);
      itmp = (vector signed long long)spu_rl((vector signed int)hib3, 1);
      hib3 = spu_and(hib3,itmp);

      s1y[0] = spu_and(s1y[0], (v_real)lowb2);
      s0y[0] = spu_and(s0y[0], (v_real)lowb2);
      s1z[0] = spu_and(s1z[0], (v_real)lowb3);
      s0z[0] = spu_and(s0z[0], (v_real)lowb3);
     
      s1y[4] = spu_and(s1y[4], (v_real)hib2);
      s0y[4] = spu_and(s0y[4], (v_real)hib2);
      s1z[4] = spu_and(s1z[4], (v_real)hib3);
      s0z[4] = spu_and(s0z[4], (v_real)hib3);
      */

      // I'm pretty sure this is exactly the same as the
      // branching in the serial code. Will test in 
      // sse2 code.

      int l2min = 1 + ((dfy[0] | dfy[1]) >> 1),
	l2max = 3 + (((dfy[0]>>1)^dfy[0]) | ((dfy[1]>>1)^dfy[1])),
	l3min = 1 + ((dfz[0] | dfz[1]) >> 1),
	l3max = 3 + (((dfz[0]>>1)^dfz[0]) | ((dfz[1]>>1)^dfz[1]));
      
      v_real fnqx, fnqy, fnqz;
	  
      fnqx = spu_mul(wni, fnqs);
      fnqx = spu_mul(qni, fnqx);
      fnqx = spu_mul(vxi, fnqx);

      fnqy = spu_mul(wni, fnqys);
      fnqy = spu_mul(qni, fnqy);
    
      fnqz = spu_mul(wni, fnqzs);
      fnqz = spu_mul(qni, fnqz);
        
      v_real jzh[5]; // As per Will's suggestion, using the minimal 
                     // variable here, and for jyh below, gives a nice
                     // performance boost
      memset(jzh,0,5*sizeof(v_real));
      v_real jx_store[VEC_SIZE], jz_store[VEC_SIZE];
      for(int mp = 0; mp<VEC_SIZE;mp++){
	jx_store[mp] = spu_splats(0.0);
	jz_store[mp] = spu_splats(0.0);
      }
      v_real jyh;
      for(int l3i=l3min; l3i<=l3max; l3i++){
	jyh = spu_splats(0.0);
	for(int l2i=l2min; l2i<=l2max; l2i++){
	  v_real wx, wy, wz;
	  wx = spu_mul(half,s1y[l2i]);
	  wx = spu_add(s0y[l2i], wx);
	  wx = spu_mul(s0z[l3i], wx);
	  tmpx = spu_mul(half, s0y[l2i]);
	  tmpy = spu_mul(third, s1y[l2i]);
	  tmpx = spu_add(tmpx, tmpy);
	  tmpx = spu_mul(tmpx, s1z[l3i]);
	  wx = spu_add(wx, tmpx);
	  
	  wy = spu_mul(half, s1z[l3i]);
	  wy = spu_add(s0z[l3i], wy);
	  wy = spu_mul(s1y[l2i], wy);
	  
	  wz = spu_mul(half, s1y[l2i]);
	  wz = spu_add(s0y[l2i], wz);
	  wz = spu_mul(s1z[l3i], wz);
	  
	  wx = spu_mul(fnqx, wx);
	  wy = spu_mul(fnqy, wy);
	  jyh = spu_sub(jyh, wy);
	  wz = spu_mul(fnqz, wz);
	  jzh[l2i] = spu_sub(jzh[l2i], wz);
	  
	  jx_store[0] = spu_shuffle(wx, jyh, uplo_pat);
	  jx_store[1] = spu_shuffle(wx, jyh, uphi_pat);
	  
	  jz_store[0] = spu_shuffle(jzh[l2i],one,uplo_pat);
	  jz_store[1] = spu_shuffle(jzh[l2i],one,uphi_pat);
	  
	  JSX(j2_scal[0] + l2i - 2, j3_scal[0] + l3i - 2 ) += jx_store[0];
	  JSZ(j2_scal[0] + l2i - 2, j3_scal[0] + l3i - 2 ) += jz_store[0];
	  JSX(j2_scal[1] + l2i - 2, j3_scal[1] + l3i - 2 ) += jx_store[1];
	  JSZ(j2_scal[1] + l2i - 2 ,j3_scal[1] + l3i - 2 ) += jz_store[1];
	}
      }    

    }
#define JOX(indx2, indx3) out[(((indx3) - wb_lg[2])*psc.img[1] + ((indx2) - wb_lg[1]))*s_NRCURR + s_JX]
  if(current_count == (total_count - 1)){
    memcpy(p_output_buffer, (void *)s_cur, 4*fld_size*sizeof(fields_cbe_real_t));       
    v_real * out = (v_real *)p_output_buffer;
    //p_output_buffer = s_cur;
  }
  
  return 0;
}

ALF_ACCEL_EXPORT_API_LIST_BEGIN
   ALF_ACCEL_EXPORT_API("cbe_push_part_yz_kernel", cbe_push_part_yz_kernel);
   ALF_ACCEL_EXPORT_API("push_yz_prep_in_dtl", push_yz_prep_in_dtl);
   ALF_ACCEL_EXPORT_API("push_yz_prep_out_dtl", push_yz_prep_out_dtl);
ALF_ACCEL_EXPORT_API_LIST_END

#undef print_vec
