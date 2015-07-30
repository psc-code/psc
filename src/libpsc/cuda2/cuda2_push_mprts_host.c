
#include "psc_cuda2.h"
#include "psc_debug.h"

#include "psc_particles_cuda2.h"

#include "psc_fields_cuda2.h"

#include "../cuda/psc_cuda.h"

#define MAX_NR_KINDS (10)

struct params_1vb {
  particle_cuda2_real_t dt;
  particle_cuda2_real_t fnqs, fnqxs, fnqys, fnqzs;
  particle_cuda2_real_t dxi[3];
  particle_cuda2_real_t dq_kind[MAX_NR_KINDS];
  particle_cuda2_real_t fnqx_kind[MAX_NR_KINDS];
  particle_cuda2_real_t fnqy_kind[MAX_NR_KINDS];
  particle_cuda2_real_t fnqz_kind[MAX_NR_KINDS];
};

static struct params_1vb prm;

static void
params_1vb_set(struct psc *psc, int p)
{
  prm.dt = ppsc->dt;
  prm.fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  prm.fnqxs = ppsc->patch[p].dx[0] * prm.fnqs / prm.dt;
  prm.fnqys = ppsc->patch[p].dx[1] * prm.fnqs / prm.dt;
  prm.fnqzs = ppsc->patch[p].dx[2] * prm.fnqs / prm.dt;
  for (int d = 0; d < 3; d++) {
    prm.dxi[d] = 1.f / ppsc->patch[p].dx[d];
  }

  assert(ppsc->nr_kinds <= MAX_NR_KINDS);
  for (int k = 0; k < ppsc->nr_kinds; k++) {
    prm.dq_kind[k] = .5f * ppsc->coeff.eta * prm.dt * ppsc->kinds[k].q / ppsc->kinds[k].m;
    prm.fnqx_kind[k] = prm.fnqxs * ppsc->kinds[k].q;
    prm.fnqy_kind[k] = prm.fnqys * ppsc->kinds[k].q;
    prm.fnqz_kind[k] = prm.fnqzs * ppsc->kinds[k].q;
  }
}

#define INTERPOLATE_1ST_EC_YZ(pf, exq, eyq, ezq, hxq, hyq, hzq)        	\
  do {									\
    particle_cuda2_real_t g0y = 1.f - og[1];				\
    particle_cuda2_real_t g0z = 1.f - og[2];				\
    particle_cuda2_real_t g1y = og[1];					\
    particle_cuda2_real_t g1z = og[2];					\
    									\
    exq = (g0z*(g0y*F3_CUDA2(pf, EX, 0,lg[1]  ,lg[2]  ) +		\
		g1y*F3_CUDA2(pf, EX, 0,lg[1]+1,lg[2]  )) +		\
	   g1z*(g0y*F3_CUDA2(pf, EX, 0,lg[1]  ,lg[2]+1) +		\
		g1y*F3_CUDA2(pf, EX, 0,lg[1]+1,lg[2]+1)));		\
    eyq = (g0z*F3_CUDA2(pf, EY, 0,lg[1]  ,lg[2]  ) +			\
	   g1z*F3_CUDA2(pf, EY, 0,lg[1]  ,lg[2]+1));			\
    ezq = (g0y*F3_CUDA2(pf, EZ, 0,lg[1]  ,lg[2]  ) +			\
	   g1y*F3_CUDA2(pf, EZ, 0,lg[1]+1,lg[2]  ));			\
									\
    hxq = F3_CUDA2(pf, HX, 0,lg[1]  ,lg[2]  );				\
    hyq = (g0y*F3_CUDA2(pf, HY, 0,lg[1]  ,lg[2]  ) +			\
	   g1y*F3_CUDA2(pf, HY, 0,lg[1]+1,lg[2]  ));			\
    hzq = (g0z*F3_CUDA2(pf, HZ, 0,lg[1]  ,lg[2]  ) +			\
	   g1z*F3_CUDA2(pf, HZ, 0,lg[1]  ,lg[2]+1));			\
    									\
    assert_finite(exq); assert_finite(eyq); assert_finite(ezq);		\
    assert_finite(hxq); assert_finite(hyq); assert_finite(hzq);		\
  } while (0)

#define INTERPOLATE_1ST_EC_XYZ(pf, exq, eyq, ezq, hxq, hyq, hzq)	\
  do {									\
    particle_cuda2_real_t g0x = 1.f - og[0];				\
    particle_cuda2_real_t g0y = 1.f - og[1];				\
    particle_cuda2_real_t g0z = 1.f - og[2];				\
    particle_cuda2_real_t g1x = og[0];					\
    particle_cuda2_real_t g1y = og[1];					\
    particle_cuda2_real_t g1z = og[2];					\
    									\
    exq = (g0z*(g0y*F3_CUDA2(pf, EX, lg[0]  ,lg[1]  ,lg[2]  ) +		\
		g1y*F3_CUDA2(pf, EX, lg[0]  ,lg[1]+1,lg[2]  )) +	\
	   g1z*(g0y*F3_CUDA2(pf, EX, lg[0]  ,lg[1]  ,lg[2]+1) +		\
		g1y*F3_CUDA2(pf, EX, lg[0]  ,lg[1]+1,lg[2]+1)));	\
									\
    eyq = (g0x*(g0z*F3_CUDA2(pf, EY, lg[0]  ,lg[1]  ,lg[2]  ) +		\
		g1z*F3_CUDA2(pf, EY, lg[0]  ,lg[1]  ,lg[2]+1)) +	\
	   g1x*(g0z*F3_CUDA2(pf, EY, lg[0]+1,lg[1]  ,lg[2]  ) +		\
		g1z*F3_CUDA2(pf, EY, lg[0]+1,lg[1]  ,lg[2]+1)));	\
									\
    ezq = (g0y*(g0x*F3_CUDA2(pf, EZ, lg[0]  ,lg[1]  ,lg[2]  ) +		\
		g1x*F3_CUDA2(pf, EZ, lg[0]+1,lg[1]  ,lg[2]  )) +	\
	   g1y*(g0x*F3_CUDA2(pf, EZ, lg[0]  ,lg[1]+1,lg[2]  ) +		\
		g1x*F3_CUDA2(pf, EZ, lg[0]+1,lg[1]+1,lg[2]  )));	\
									\
    hxq = (g0x*F3_CUDA2(pf, HX, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	   g1x*F3_CUDA2(pf, HX, lg[0]+1,lg[1]  ,lg[2]  ));		\
									\
    hyq = (g0y*F3_CUDA2(pf, HY, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	   g1y*F3_CUDA2(pf, HY, lg[0]  ,lg[1]+1,lg[2]  ));		\
									\
    hzq = (g0z*F3_CUDA2(pf, HZ, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	   g1z*F3_CUDA2(pf, HZ, lg[0]  ,lg[1]  ,lg[2]+1));		\
									\
    assert_finite(exq); assert_finite(eyq); assert_finite(ezq);		\
    assert_finite(hxq); assert_finite(hyq); assert_finite(hzq);		\
  } while (0)


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
push_xi_yz(particle_cuda2_t *part, particle_cuda2_real_t vxi[3], particle_cuda2_real_t dt)
{
  part->xi4.y += vxi[1] * dt;
  part->xi4.z += vxi[2] * dt;
}

static inline void
push_xi_xyz(particle_cuda2_t *part, particle_cuda2_real_t vxi[3], particle_cuda2_real_t dt)
{
  part->xi4.x += vxi[0] * dt;
  part->xi4.y += vxi[1] * dt;
  part->xi4.z += vxi[2] * dt;
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

// ======================================================================
// current 1vb (yz)

static inline void
calc_3d_dx1_yz(particle_cuda2_real_t dx1[3], particle_cuda2_real_t x[3], particle_cuda2_real_t dx[3], int off[3])
{
  if (off[2] == 0) {
    dx1[1] = .5f * off[1] - x[1];
   if (dx[1] == 0.f) {
     dx1[0] = 0.f;
     dx1[2] = 0.f;
   } else {
     dx1[0] = dx[0] / dx[1] * dx1[1];
     dx1[2] = dx[2] / dx[1] * dx1[1];
   }
  } else {
    dx1[2] = .5f * off[2] - x[2];
   if (dx[2] == 0.f) {
     dx1[0] = 0.f;
     dx1[1] = 0.f;
   } else {
     dx1[0] = dx[0] / dx[2] * dx1[2];
     dx1[1] = dx[1] / dx[2] * dx1[2];
   }
  }
}

static inline void
curr_3d_vb_cell_yz(struct psc_fields *pf, int i[3], particle_cuda2_real_t x[3], particle_cuda2_real_t dx[3],
		   particle_cuda2_real_t fnq[3], particle_cuda2_real_t dxt[3], int off[3])
{
  particle_cuda2_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];
  particle_cuda2_real_t xa[3] = { 0.,
			    x[1] + .5f * dx[1],
			    x[2] + .5f * dx[2], };
  F3_CUDA2(pf, JXI, 0,i[1]  ,i[2]  ) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h);
  F3_CUDA2(pf, JXI, 0,i[1]+1,i[2]  ) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h);
  F3_CUDA2(pf, JXI, 0,i[1]  ,i[2]+1) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h);
  F3_CUDA2(pf, JXI, 0,i[1]+1,i[2]+1) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h);

  F3_CUDA2(pf, JYI, 0,i[1]  ,i[2]  ) += fnq[1] * dx[1] * (.5f - xa[2]);
  F3_CUDA2(pf, JYI, 0,i[1]  ,i[2]+1) += fnq[1] * dx[1] * (.5f + xa[2]);
  F3_CUDA2(pf, JZI, 0,i[1]  ,i[2]  ) += fnq[2] * dx[2] * (.5f - xa[1]);
  F3_CUDA2(pf, JZI, 0,i[1]+1,i[2]  ) += fnq[2] * dx[2] * (.5f + xa[1]);
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    dxt[2] -= dx[2];
    x[1] += dx[1] - off[1];
    x[2] += dx[2] - off[2];
    i[1] += off[1];
    i[2] += off[2];
  }
}

static inline void
calc_j_3d_yz(struct psc_fields *flds, particle_cuda2_real_t *xm, particle_cuda2_real_t *xp,
		int *lf, int *lg, particle_cuda2_t *prt, particle_cuda2_real_t *vxi)
{
  int i[3] = { 0, lg[1], lg[2] };					
  int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			
  particle_cuda2_real_t dx[3] = { 0., xp[1] - xm[1], xp[2] - xm[2] };		
  particle_cuda2_real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; 
  									
  particle_cuda2_real_t dx1[3];
  int off[3];
  int first_dir, second_dir = -1;
  /* FIXME, make sure we never div-by-zero? */
  if (idiff[1] == 0 && idiff[2] == 0) {
    first_dir = -1;
  } else if (idiff[1] == 0) {
    first_dir = 2;
  } else if (idiff[2] == 0) {
    first_dir = 1;
  } else {
    dx1[1] = .5f * idiff[1] - x[1];
    if (dx[1] == 0.f) {
      dx1[2] = 0.f;
    } else {
      dx1[2] = dx[2] / dx[1] * dx1[1];
    }
    if (particle_cuda2_real_abs(x[2] + dx1[2]) > .5f) {
      first_dir = 2;
    } else {
      first_dir = 1;
    }
    second_dir = 3 - first_dir;
  }

  int kind = cuda_float_as_int(prt->xi4.w);
  particle_cuda2_real_t fnq[3] = { particle_cuda2_wni(prt) * prm.fnqx_kind[kind],
				   particle_cuda2_wni(prt) * prm.fnqy_kind[kind],
				   particle_cuda2_wni(prt) * prm.fnqz_kind[kind] };
  dx[0] = vxi[0] * prm.dt * prm.dxi[0];

  if (first_dir >= 0) {
    off[3 - first_dir] = 0;
    off[first_dir] = idiff[first_dir];
    calc_3d_dx1_yz(dx1, x, dx, off);
    curr_3d_vb_cell_yz(flds, i, x, dx1, fnq, dx, off);
  }

  if (second_dir >= 0) {
    off[first_dir] = 0;
    off[second_dir] = idiff[second_dir];
    calc_3d_dx1_yz(dx1, x, dx, off);
    curr_3d_vb_cell_yz(flds, i, x, dx1, fnq, dx, off);
  }

  curr_3d_vb_cell_yz(flds, i, x, dx, fnq, NULL, NULL);
}

// ======================================================================

static inline void
curr_3d_vb_one_yz(struct psc_fields *pf, int i[3], particle_cuda2_real_t x[3], particle_cuda2_real_t dx[3],
		  particle_cuda2_real_t fnq[3], particle_cuda2_real_t dxt[3], int off[3])
{
  particle_cuda2_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];
  particle_cuda2_real_t xa[3] = { 0.,
			    x[1] + .5f * dx[1],
			    x[2] + .5f * dx[2], };
  F3_CUDA2(pf, JXI, 0,i[1]  ,i[2]  ) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h);
  F3_CUDA2(pf, JXI, 0,i[1]+1,i[2]  ) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h);
  F3_CUDA2(pf, JXI, 0,i[1]  ,i[2]+1) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h);
  F3_CUDA2(pf, JXI, 0,i[1]+1,i[2]+1) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h);

  F3_CUDA2(pf, JYI, 0,i[1]  ,i[2]  ) += fnq[1] * dx[1] * (.5f - xa[2]);
  F3_CUDA2(pf, JYI, 0,i[1]  ,i[2]+1) += fnq[1] * dx[1] * (.5f + xa[2]);
  F3_CUDA2(pf, JZI, 0,i[1]  ,i[2]  ) += fnq[2] * dx[2] * (.5f - xa[1]);
  F3_CUDA2(pf, JZI, 0,i[1]+1,i[2]  ) += fnq[2] * dx[2] * (.5f + xa[1]);
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    dxt[2] -= dx[2];
    x[1] += dx[1] - off[1];
    x[2] += dx[2] - off[2];
    i[1] += off[1];
    i[2] += off[2];
  }
}

static inline void
calc_j2_3d_yz(struct psc_fields *flds, particle_cuda2_real_t *xm, particle_cuda2_real_t *xp,
	      int *lf, int *lg, particle_cuda2_t *prt, particle_cuda2_real_t *vxi)
{
  int i[3] = { 0, lg[1], lg[2] };					
  int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			
  particle_cuda2_real_t dxt[3] = { 0., xp[1] - xm[1], xp[2] - xm[2] };		
  particle_cuda2_real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; 
  									
  int kind = cuda_float_as_int(prt->xi4.w);
  particle_cuda2_real_t fnq[3] = { particle_cuda2_wni(prt) * prm.fnqx_kind[kind],
				   particle_cuda2_wni(prt) * prm.fnqy_kind[kind],
				   particle_cuda2_wni(prt) * prm.fnqz_kind[kind] };
  dxt[0] = vxi[0] * prm.dt * prm.dxi[0];

  particle_cuda2_real_t dx[3];
  int off[3];
  int first_dir, second_dir = -1;

  if (idiff[2] == 0) { // not intersecting z
    if (idiff[1] == 0) { // not intersecting y
      curr_3d_vb_one_yz(flds, i, x, dxt, fnq, NULL, NULL);
    } else { // intersecting y
      //first_dir = 1;
      off[2] = 0;
      off[1] = idiff[1];
      calc_3d_dx1_yz(dx, x, dxt, off);
      curr_3d_vb_one_yz(flds, i, x, dx, fnq, dxt, off);
      curr_3d_vb_one_yz(flds, i, x, dxt, fnq, NULL, NULL);
    }
  } else { // intersecting z
    if (idiff[1] == 0) { // not intersecting y
      off[1] = 0;
      off[2] = idiff[2];
      calc_3d_dx1_yz(dx, x, dxt, off);
      curr_3d_vb_one_yz(flds, i, x, dx, fnq, dxt, off);
      curr_3d_vb_one_yz(flds, i, x, dxt, fnq, NULL, NULL);
    } else { // intersecting y
      dx[1] = .5f * idiff[1] - x[1];
      if (dxt[1] == 0.f) {
	dx[2] = 0.f;
      } else {
	dx[2] = dxt[2] / dxt[1] * dx[1];
      }
      if (particle_cuda2_real_abs(x[2] + dx[2]) > .5f) {
	first_dir = 2;
	off[1] = 0;
	off[2] = idiff[2];
	calc_3d_dx1_yz(dx, x, dxt, off);
	curr_3d_vb_one_yz(flds, i, x, dx, fnq, dxt, off);
      } else {
	first_dir = 1;
	off[2] = 0;
	off[1] = idiff[1];
	calc_3d_dx1_yz(dx, x, dxt, off);
	curr_3d_vb_one_yz(flds, i, x, dx, fnq, dxt, off);
      }
      second_dir = 3 - first_dir;
      off[first_dir] = 0;
      off[second_dir] = idiff[second_dir];
      calc_3d_dx1_yz(dx, x, dxt, off);
      curr_3d_vb_one_yz(flds, i, x, dx, fnq, dxt, off);
      curr_3d_vb_one_yz(flds, i, x, dxt, fnq, NULL, NULL);
    }
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// push_one_yz

static void
push_one_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds, int n, int p)
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
  
  // FIELD INTERPOLATION
  particle_cuda2_real_t exq, eyq, ezq, hxq, hyq, hzq;
  INTERPOLATE_1ST_EC_YZ(flds, exq, eyq, ezq, hxq, hyq, hzq);
  
  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
  int kind = cuda_float_as_int(prt.xi4.w);
  particle_cuda2_real_t dq = prm.dq_kind[kind];
  push_pxi(&prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  
  particle_cuda2_real_t vxi[3];
  calc_vxi(vxi, &prt);
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  push_xi_yz(&prt, vxi, prm.dt);
  
  int lf[3];
  particle_cuda2_real_t of[3], xp[3];
  find_idx_off_pos_1st_rel(&prt.xi4.x, lf, of, xp, 0.f, prm.dxi);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
  calc_j2_3d_yz(flds, xm, xp, lf, lg, &prt, vxi);

  _STORE_PARTICLE_POS(prt, mprts_sub->h_xi4, n);
  _STORE_PARTICLE_MOM(prt, mprts_sub->h_pxi4, n);
}

// ----------------------------------------------------------------------
// push_one_xyz

static void
push_one_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds, int n, int p)
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
  
  // FIELD INTERPOLATION
  particle_cuda2_real_t exq, eyq, ezq, hxq, hyq, hzq;
  INTERPOLATE_1ST_EC_XYZ(flds, exq, eyq, ezq, hxq, hyq, hzq);
  
  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
  int kind = cuda_float_as_int(prt.xi4.w);
  particle_cuda2_real_t dq = prm.dq_kind[kind];
  push_pxi(&prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  
  particle_cuda2_real_t vxi[3];
  calc_vxi(vxi, &prt);
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  push_xi_xyz(&prt, vxi, prm.dt);
  
  int lf[3];
  particle_cuda2_real_t of[3], xp[3];
  find_idx_off_pos_1st_rel(&prt.xi4.x, lf, of, xp, 0.f, prm.dxi);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
#warning should be xyz
  //  calc_j_3d_yz(flds, xm, xp, lf, lg, &prt, vxi);

  _STORE_PARTICLE_POS(prt, mprts_sub->h_xi4, n);
  _STORE_PARTICLE_MOM(prt, mprts_sub->h_pxi4, n);
}

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts_yz_gold

void
cuda2_1vbec_push_mprts_yz_gold(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  params_1vb_set(ppsc, 0);

  psc_mfields_zero_range(mflds, JXI, JXI + 3);

  for (int b = 0; b < mprts_sub->nr_blocks_total; b++) {
    int p = b / mprts_sub->nr_blocks;
    for (int n = mprts_sub->h_b_off[b]; n < mprts_sub->h_b_off[b+1]; n++) {
      push_one_yz(mprts, mflds, n, p);
    }
  }
}

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts_xyz_gold

void
cuda2_1vbec_push_mprts_xyz_gold(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  params_1vb_set(ppsc, 0);

  psc_mfields_zero_range(mflds, JXI, JXI + 3);

  for (int b = 0; b < mprts_sub->nr_blocks_total; b++) {
    int p = b / mprts_sub->nr_blocks;
    for (int n = mprts_sub->h_b_off[b]; n < mprts_sub->h_b_off[b+1]; n++) {
      push_one_xyz(mprts, mflds, n, p);
    }
  }
}

