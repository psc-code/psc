
#include "psc_debug.h"

#include "psc_particles_as_cuda2.h"

#include "psc_fields_cuda2.h"

#include "../cuda/psc_cuda.h"

#define F3_CACHE F3_CUDA2

#include "../psc_push_particles/inc_params.c"
#include "../psc_push_particles/inc_interpolate.c"

static inline void
calc_vxi(particle_cuda2_real_t vxi[3], particle_cuda2_t *part)
{
  particle_cuda2_real_t root = 1.f 
    / particle_cuda2_real_sqrt(1.f + sqr(part->pxi4.x) + sqr(part->pxi4.y) + sqr(part->pxi4.z));
  vxi[0] = part->pxi4.x * root;
  vxi[1] = part->pxi4.y * root;
  vxi[2] = part->pxi4.z * root;
}

static inline void
push_xi(particle_cuda2_t *part, particle_cuda2_real_t vxi[3], particle_cuda2_real_t dt)
{
#if DIM == DIM_YZ
  part->xi4.y += vxi[1] * dt;
  part->xi4.z += vxi[2] * dt;
#elif DIM == DIM_XYZ
  part->xi4.x += vxi[0] * dt;
  part->xi4.y += vxi[1] * dt;
  part->xi4.z += vxi[2] * dt;
#endif
}

static inline void
push_pxi(particle_cuda2_t *part, particle_cuda2_real_t exq, particle_cuda2_real_t eyq, particle_cuda2_real_t ezq,
	 particle_cuda2_real_t hxq, particle_cuda2_real_t hyq, particle_cuda2_real_t hzq, particle_cuda2_real_t dq)
{
  particle_cuda2_real_t pxm = part->pxi4.x + dq*exq;
  particle_cuda2_real_t pym = part->pxi4.y + dq*eyq;
  particle_cuda2_real_t pzm = part->pxi4.z + dq*ezq;
  
  particle_cuda2_real_t root = dq / particle_cuda2_real_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
  particle_cuda2_real_t taux = hxq*root;
  particle_cuda2_real_t tauy = hyq*root;
  particle_cuda2_real_t tauz = hzq*root;
  
  particle_cuda2_real_t tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
  particle_cuda2_real_t pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
	       (2.f*taux*tauy+2.f*tauz)*pym + 
	       (2.f*taux*tauz-2.f*tauy)*pzm)*tau;
  particle_cuda2_real_t pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
	       (1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
	       (2.f*tauy*tauz+2.f*taux)*pzm)*tau;
  particle_cuda2_real_t pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
	       (2.f*tauy*tauz-2.f*taux)*pym +
	       (1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
  
  part->pxi4.x = pxp + dq * exq;
  part->pxi4.y = pyp + dq * eyq;
  part->pxi4.z = pzp + dq * ezq;
}

static inline void
find_idx_off_1st_rel(particle_cuda2_real_t xi[3], int lg[3], particle_cuda2_real_t og[3], particle_cuda2_real_t shift,
		     particle_cuda2_real_t dxi[3])
{
  for (int d = 0; d < 3; d++) {
    particle_cuda2_real_t pos = xi[d] * dxi[d] + shift;
    lg[d] = particle_cuda2_real_fint(pos);
    og[d] = pos - lg[d];
  }
}

static inline void
find_idx_off_pos_1st_rel(particle_cuda2_real_t xi[3], int lg[3], particle_cuda2_real_t og[3],
			 particle_cuda2_real_t pos[3], particle_cuda2_real_t shift,
			 particle_cuda2_real_t dxi[3])
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * dxi[d] + shift;
    lg[d] = particle_cuda2_real_fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

#if CALC_J == CALC_J_1VB_SPLIT
#include "cuda2_calc_j_split_inc.c"
#elif CALC_J == CALC_J_1VB_VAR1
#include "cuda2_calc_j_var1_inc.c"
#endif

// ----------------------------------------------------------------------
// push_one

static void
push_one(struct psc_mparticles *mprts, struct psc_mfields *mflds, int n, int p)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  struct psc_fields *flds = psc_mfields_get_patch(mflds, p);

  particle_cuda2_t prt;
  _LOAD_PARTICLE_POS(prt, mprts_sub->h_xi4, n);
  _LOAD_PARTICLE_MOM(prt, mprts_sub->h_pxi4, n);
  
  // field interpolation
  
  int lg[3], lh[3];
  particle_cuda2_real_t og[3], oh[3], xm[3];
  find_idx_off_pos_1st_rel(&prt.xi4.x, lg, og, xm, 0.f, prm.dxi); // FIXME passing xi hack
  find_idx_off_1st_rel(&prt.xi4.x, lh, oh, -.5f, prm.dxi);
  
  particle_cuda2_real_t exq, eyq, ezq, hxq, hyq, hzq;
  INTERPOLATE_1ST_EC(flds, exq, eyq, ezq, hxq, hyq, hzq);
  
  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
  int kind = cuda_float_as_int(prt.xi4.w);
  particle_cuda2_real_t dq = prm.dq_kind[kind];
  push_pxi(&prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  
  particle_cuda2_real_t vxi[3];
  calc_vxi(vxi, &prt);
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  push_xi(&prt, vxi, prm.dt);
  
  int lf[3];
  particle_cuda2_real_t of[3], xp[3];
  find_idx_off_pos_1st_rel(&prt.xi4.x, lf, of, xp, 0.f, prm.dxi);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
  calc_j(flds, xm, xp, lf, lg, &prt, vxi);

  _STORE_PARTICLE_POS(prt, mprts_sub->h_xi4, n);
  _STORE_PARTICLE_MOM(prt, mprts_sub->h_pxi4, n);
}

