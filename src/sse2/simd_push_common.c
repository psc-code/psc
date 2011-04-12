#include "psc_sse2.h"
#include "simd_wrap.h"
#include "simd_push_common.h"
#include <math.h>

////////////////////
/// For SIMD code, calculates the velocities from momenta.
///
/// 'root' is used to update p2X, so need to have it
/// stored to a variable in the push_part_yx_ function.
inline void 
calc_vi(struct particle_vec *p, pvReal * restrict vxi, pvReal * restrict vyi,
	pvReal * restrict vzi, pvReal * restrict root){
  pvReal tmpx, tmpy, tmpz;
  
  tmpx.r = pv_mul_real(p->pxi.r, p->pxi.r);
  tmpy.r = pv_mul_real(p->pyi.r, p->pyi.r);
  tmpz.r = pv_mul_real(p->pzi.r, p->pzi.r);
  
  tmpx.r = pv_add_real(tmpx.r, tmpy.r);
  tmpz.r = pv_add_real(ones.r, tmpz.r);
  tmpx.r = pv_add_real(tmpx.r, tmpz.r);
  root->r = pv_sqrt_real(tmpx.r);
                                             
  vxi->r = pv_div_real(p->pxi.r, root->r);
  vyi->r = pv_div_real(p->pyi.r, root->r);
  vzi->r = pv_div_real(p->pzi.r, root->r);
}

///////////////////////////
///  For SIMD code, finds the index corresponding to position xi using the 
/// C99 built-in round.
///
/// This is slower, and should really only be used
/// when there is no weighting to suppress strange behavior
/// in the intrinsic round find_index_Iround.
inline void
find_index_Cround(pvReal *xi, ///< Particle coordinate
		  pvReal *dxi, ///< 1/dx
		  pvInt * restrict j, ///< Resultant coordinate on the grid
		  pvReal * restrict h) ///< Difference between j - xi*dxi
{ 
  pvReal tmp;
  tmp.r = pv_mul_real(xi->r, dxi->r);
  for(int m = 0; m < VEC_SIZE; m++){
      j->v[m] = round(tmp.v[m]);
  }
  h->r = pv_cvt_int_to_real(j->r);
  h->r = pv_sub_real(h->r, tmp.r);
}

///////////////
/// For SIMD code, finds index using intrinsic round.
///
/// Sometimes .5 is treated as .5-epsilon (not sure why)
/// and is rounded the wrong direction. This isn't a problem
/// most of the time due to the weighting coeffecients, but
/// can be a problem for any unweighted dimension. If such a 
/// problem occurs, using find_index_Cround should be used.
inline void
find_index_Iround(pvReal *xi, ///< Particle coordinate
		  pvReal *dxi, ///< 1/dx
		  pvInt * restrict j,///< Resultant coordinate on the grid
		  pvReal * restrict h)///< Difference between j - xi*dxi
{
  pvReal tmp;
  tmp.r = pv_mul_real(xi->r, dxi->r);
  j->r = pv_cvt_real_to_int(tmp.r);
  h->r = pv_cvt_int_to_real(j->r);
  h->r = pv_sub_real(h->r, tmp.r);
}

/////////////////
/// For SIMD code, same as find_index_Iround
/// except it subtracts the vector 'shift' before rounding.
///
/// Always uses the intrinsic round. Currently no C round version
/// exists, and shouldn't be necessary. Calculates:
/// j = round(xi*dxi - shift).
/// I suppose there isn't a good reason why this has it's own function...
inline void
find_index_minus_shift(pvReal *xi, ///< Particle coordinate
		       pvReal *dxi, ///< 1/dx
		       pvInt * restrict j, ///< Resultant coordinate on the grid
		       pvReal * restrict h, ///< Difference between j - xi*dxi
		       pvReal *shift){ ///< Amount to be subtracted from xi*dxi before finding index
  pvReal tmp;
  tmp.r = pv_mul_real(xi->r, dxi->r);
  tmp.r = pv_sub_real(tmp.r, shift->r);
  j->r = pv_cvt_real_to_int(tmp.r);
  h->r = pv_cvt_int_to_real(j->r);
  h->r = pv_sub_real(h->r, tmp.r);
}

///////
/// For SIMD code, finds the first (m) weight factor
/// for field interpolation on the grid
inline void
ip_to_grid_m(pvReal *h, pvReal * restrict xmx)
{
  xmx->r = pv_add_real(half.r, h->r);
  xmx->r = pv_mul_real(xmx->r, xmx->r);
  xmx->r = pv_mul_real(half.r, xmx->r);
}

///////
/// For SIMD code, finds the second (O here, 0 in FORTRAN) weight factor
/// for field interpolation on the grid
inline void
ip_to_grid_O(pvReal *h, pvReal * restrict xOx)
{
  xOx->r = pv_mul_real(h->r, h->r);
  xOx->r = pv_sub_real(threefourths.r, xOx->r);
}

/////////////
/// For SIMD code, finds the third (l here, 1 in FORTRAN) weight factor
/// for field interpolation on the grid

inline void
ip_to_grid_l(pvReal *h, pvReal * restrict xlx)
{
  xlx->r = pv_sub_real(half.r, h->r);
  xlx->r = pv_mul_real(xlx->r, xlx->r);
  xlx->r = pv_mul_real(half.r, xlx->r);
}

//---------------------------------------------
// form_factor_... - I'm not sure my understanding 
// of this is correct. I believe these calculate 
// the extent of the quasi-particle on the grid
// before updating the charge and current densities

inline void
form_factor_m(pvReal *h, pvReal * restrict xmx)
{
  xmx->r = pv_sub_real(h->r,ones.r);
  xmx->r = pv_add_real(onepfive.r, xmx->r); // h-1 always <0
  xmx->r = pv_mul_real(xmx->r, xmx->r);
  xmx->r = pv_mul_real(half.r, xmx->r);
}

inline void
form_factor_O(pvReal *h, pvReal * restrict xOx)
{
  xOx->r = pv_mul_real(h->r, h->r);
  xOx->r = pv_sub_real(threefourths.r, xOx->r);
}

inline void
form_factor_l(pvReal *h, pvReal * restrict xlx)
{
  xlx->r = pv_add_real(ones.r, h->r); 
  xlx->r = pv_sub_real(onepfive.r, xlx->r);//h+1 always >0
  xlx->r = pv_mul_real(xlx->r, xlx->r);
  xlx->r = pv_mul_real(half.r, xlx->r);
}

////////////////////////
/// For SIMD code, dvances the momentum a full time
/// step using the interpolated fields.
///
/// As far as I can tell, this should work the same
/// regardless of the dimensionality of the system.
/// If you'd rather use h{x,y,z} instead of b{x,y,z}, 
/// you should just be able to pass them instead.
/// This piece of the code is very complicated in
/// SIMD, but from what I can tell it's bug free 
/// and runs pretty quick.
void
push_pi_dt(struct particle_vec * p,
	   pvReal * exq, pvReal * eyq,  pvReal * ezq,
	   pvReal * bxq, pvReal * byq,  pvReal * bzq,
	   pvReal * dqs)
{

  pvReal tmpx, tmpy, tmpz, root;
  // Half step momentum  with E-field

  pvReal dq, dqex, dqey, dqez;
  dq.r = pv_div_real(dqs->r, p->mni.r);
  dq.r = pv_mul_real(p->qni.r, dq.r);
  
  dqex.r = pv_mul_real(dq.r, exq->r);
  dqey.r = pv_mul_real(dq.r, eyq->r);
  dqez.r = pv_mul_real(dq.r, ezq->r);
  
  p->pxi.r = pv_add_real(p->pxi.r, dqex.r);
  p->pyi.r = pv_add_real(p->pyi.r, dqey.r);
  p->pzi.r = pv_add_real(p->pzi.r, dqez.r);


// Rotate with B-field
  pvReal taux, tauy, tauz, txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;
  
  tmpx.r = pv_mul_real(p->pxi.r, p->pxi.r);
  tmpy.r = pv_mul_real(p->pyi.r, p->pyi.r);
  tmpz.r = pv_mul_real(p->pzi.r, p->pzi.r);
  root.r = pv_add_real(ones.r, tmpx.r);
  tmpy.r = pv_add_real(tmpy.r, tmpz.r);
  root.r = pv_add_real(root.r, tmpy.r);
  root.r = pv_sqrt_real(root.r);
  root.r = pv_div_real(dq.r, root.r);

  taux.r = pv_mul_real(bxq->r, root.r);
  tauy.r = pv_mul_real(byq->r, root.r);
  tauz.r = pv_mul_real(bzq->r, root.r);

  txx.r = pv_mul_real(taux.r, taux.r);
  tyy.r = pv_mul_real(tauy.r, tauy.r);
  tzz.r = pv_mul_real(tauz.r, tauz.r);
  t2xy.r = pv_mul_real(taux.r, tauy.r);
  t2xz.r = pv_mul_real(taux.r, tauz.r);
  t2yz.r = pv_mul_real(tauy.r, tauz.r);
  t2xy.r = pv_add_real(t2xy.r, t2xy.r);
  t2xz.r = pv_add_real(t2xz.r, t2xz.r);
  t2yz.r = pv_add_real(t2yz.r, t2yz.r);

  pvReal tau;

  tau.r = pv_add_real(ones.r, txx.r);
  tmpx.r = pv_add_real(tyy.r, tzz.r);
  tau.r = pv_add_real(tau.r, tmpx.r);
  tau.r = pv_div_real(ones.r, tau.r); //recp is evil! evil evil evil evil!!!!
    
  //Never use tau_ without a two in front
  taux.r = pv_add_real(taux.r, taux.r);
  tauy.r = pv_add_real(tauy.r, tauy.r);
  tauz.r = pv_add_real(tauz.r, tauz.r);
    
  //pxp
  tmpx.r = pv_add_real(ones.r, txx.r);
  tmpx.r = pv_sub_real(tmpx.r, tyy.r);
  tmpx.r = pv_sub_real(tmpx.r, tzz.r);
  tmpx.r = pv_mul_real(tmpx.r, p->pxi.r);
    
  tmpy.r = pv_add_real(t2xy.r, tauz.r);
  tmpy.r = pv_mul_real(tmpy.r, p->pyi.r);

  tmpz.r = pv_sub_real(t2xz.r, tauy.r);
  tmpz.r = pv_mul_real(tmpz.r, p->pzi.r);

  pxp.r = pv_add_real(tmpx.r, tmpy.r);
  pxp.r = pv_add_real(pxp.r, tmpz.r);
  pxp.r = pv_mul_real(pxp.r, tau.r);
    
  //pyp
  tmpx.r = pv_sub_real(t2xy.r, tauz.r);
  tmpx.r = pv_mul_real(tmpx.r, p->pxi.r);

  tmpy.r = pv_sub_real(ones.r, txx.r);
  tmpy.r = pv_add_real(tmpy.r, tyy.r);
  tmpy.r = pv_sub_real(tmpy.r, tzz.r);
  tmpy.r = pv_mul_real(tmpy.r, p->pyi.r);
  
  tmpz.r = pv_add_real(t2yz.r, taux.r);
  tmpz.r = pv_mul_real(tmpz.r, p->pzi.r);

  pyp.r = pv_add_real(tmpx.r, tmpy.r);
  pyp.r = pv_add_real(pyp.r, tmpz.r);
  pyp.r = pv_mul_real(pyp.r, tau.r);
   
  //pzp
  tmpx.r = pv_add_real(t2xz.r, tauy.r);
  tmpx.r = pv_mul_real(tmpx.r, p->pxi.r);

  tmpy.r = pv_sub_real(t2yz.r, taux.r);
  tmpy.r = pv_mul_real(tmpy.r, p->pyi.r);

  tmpz.r = pv_sub_real(ones.r, txx.r);
  tmpz.r = pv_sub_real(tmpz.r, tyy.r);
  tmpz.r = pv_add_real(tmpz.r, tzz.r);
  tmpz.r = pv_mul_real(tmpz.r, p->pzi.r);

  pzp.r = pv_add_real(tmpx.r, tmpy.r);
  pzp.r = pv_add_real(pzp.r, tmpz.r);
  pzp.r = pv_mul_real(pzp.r, tau.r);

// Half step momentum  with E-field

  p->pxi.r = pv_add_real(pxp.r, dqex.r);
  p->pyi.r = pv_add_real(pyp.r, dqey.r);
  p->pzi.r = pv_add_real(pzp.r, dqez.r);
}

////////////////////
/// \file simd_push_common.c Certain functions common to all SIMD particle pushers.
///
/// The functions contained in this file used to be declared static at the 
/// the top of sse2_push_part_yz.c. As I started writing the xz pusher, I relized 
/// these seem to be common among all the differenent pushers. As such, there's 
/// not really any harm in giving them a seperate file. This way any bug fixes
/// don't need to be replicated seven times. 
