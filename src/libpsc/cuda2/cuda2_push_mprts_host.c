
#include "psc_cuda2.h"
#include "psc_debug.h"

#include "psc_particles_single.h"

#include "psc_fields_cuda2.h"

static inline int
fint(particle_single_real_t x)
{
  return floorf(x);
}

#define MAX_NR_KINDS (10)

struct params_1vb {
  float dt;
  float fnqs, fnqxs, fnqys, fnqzs;
  float dxi[3];
  float dq_kind[MAX_NR_KINDS];
  float fnqx_kind[MAX_NR_KINDS];
  float fnqy_kind[MAX_NR_KINDS];
  float fnqz_kind[MAX_NR_KINDS];
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

#define INTERPOLATE_1ST_EC(pf, exq, eyq, ezq, hxq, hyq, hzq)        	\
  do {									\
    float g0y = 1.f - og[1];					\
    float g0z = 1.f - og[2];					\
    float g1y = og[1];					\
    float g1z = og[2];					\
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


static inline void
calc_vxi(float vxi[3], particle_single_t *part)
{
  float root = 1.f 
    / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline void
push_xi(particle_single_t *part, float vxi[3], float dt)
{
  part->yi += vxi[1] * dt;
  part->zi += vxi[2] * dt;
}

static inline void
push_pxi(particle_single_t *part, float exq, float eyq, float ezq,
	 float hxq, float hyq, float hzq, float dq)
{
  float pxm = part->pxi + dq*exq;
  float pym = part->pyi + dq*eyq;
  float pzm = part->pzi + dq*ezq;
  
  float root = dq / sqrtf(1.f + pxm*pxm + pym*pym + pzm*pzm);
  float taux = hxq*root;
  float tauy = hyq*root;
  float tauz = hzq*root;
  
  float tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
  float pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
	       (2.f*taux*tauy+2.f*tauz)*pym + 
	       (2.f*taux*tauz-2.f*tauy)*pzm)*tau;
  float pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
	       (1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
	       (2.f*tauy*tauz+2.f*taux)*pzm)*tau;
  float pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
	       (2.f*tauy*tauz-2.f*taux)*pym +
	       (1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
  
  part->pxi = pxp + dq * exq;
  part->pyi = pyp + dq * eyq;
  part->pzi = pzp + dq * ezq;
}

static inline void
find_idx_off_1st_rel(float xi[3], int lg[3], float og[3], float shift,
		     float dxi[3])
{
  for (int d = 0; d < 3; d++) {
    float pos = xi[d] * dxi[d] + shift;
    lg[d] = fint(pos);
    og[d] = pos - lg[d];
  }
}

static inline void
find_idx_off_pos_1st_rel(float xi[3], int lg[3], float og[3],
			 float pos[3], float shift,
			 float dxi[3])
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * dxi[d] + shift;
    lg[d] = fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

// ======================================================================
// current 1vb (yz)

static inline void
calc_3d_dx1(float dx1[3], float x[3], float dx[3], int off[3])
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
curr_3d_vb_cell(struct psc_fields *pf, int i[3], float x[3], float dx[3],
		float fnq[3], float dxt[3], int off[3])
{
  float h = (1.f/12.f) * dx[0] * dx[1] * dx[2];
  float xa[3] = { 0.,
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
calc_jxyz_3d(struct psc_fields *flds, float *xm, float *xp,
	     int *lf, int *lg, particle_single_t *prt, float *vxi)
{
  int i[3] = { 0, lg[1], lg[2] };					
  int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			
  float dx[3] = { 0., xp[1] - xm[1], xp[2] - xm[2] };		
  float x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; 
  									
  float dx1[3];
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
    if (fabsf(x[2] + dx1[2]) > .5f) {
      first_dir = 2;
    } else {
      first_dir = 1;
    }
    second_dir = 3 - first_dir;
  }

  float fnq[3] = { particle_single_wni(prt) * prm.fnqx_kind[prt->kind],
		   particle_single_wni(prt) * prm.fnqy_kind[prt->kind],
		   particle_single_wni(prt) * prm.fnqz_kind[prt->kind] };
  dx[0] = vxi[0] * prm.dt * prm.dxi[0];

  if (first_dir >= 0) {
    off[3 - first_dir] = 0;
    off[first_dir] = idiff[first_dir];
    calc_3d_dx1(dx1, x, dx, off);
    curr_3d_vb_cell(flds, i, x, dx1, fnq, dx, off);
  }

  if (second_dir >= 0) {
    off[first_dir] = 0;
    off[second_dir] = idiff[second_dir];
    calc_3d_dx1(dx1, x, dx, off);
    curr_3d_vb_cell(flds, i, x, dx1, fnq, dx, off);
  }

  curr_3d_vb_cell(flds, i, x, dx, fnq, NULL, NULL);
}

// ======================================================================

static void
push_one(struct psc_fields *flds, struct psc_particles *prts, int n)
{
  particle_single_t *prt = particles_single_get_one(prts, n);
  
  // field interpolation
  
  int lg[3], lh[3];
  float og[3], oh[3], xm[3];
  find_idx_off_pos_1st_rel(&prt->xi, lg, og, xm, 0.f, prm.dxi); // FIXME passing xi hack
  find_idx_off_1st_rel(&prt->xi, lh, oh, -.5f, prm.dxi);
  
  // FIELD INTERPOLATION
  float exq, eyq, ezq, hxq, hyq, hzq;
  INTERPOLATE_1ST_EC(flds, exq, eyq, ezq, hxq, hyq, hzq);
  
  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
  float dq = prm.dq_kind[prt->kind];
  push_pxi(prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  
  float vxi[3];
  calc_vxi(vxi, prt);
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  push_xi(prt, vxi, prm.dt);
  
  int lf[3];
  float of[3], xp[3];
  find_idx_off_pos_1st_rel(&prt->xi, lf, of, xp, 0.f, prm.dxi);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
  calc_jxyz_3d(flds, xm, xp, lf, lg, prt, vxi);
}

static void
psc_push_particles_push_a_yz(struct psc_particles *prts,
			     struct psc_fields *flds)
{
  params_1vb_set(ppsc, flds->p);

  psc_fields_zero_range(flds, JXI, JXI + 3);
  for (int n = 0; n < prts->n_part; n++) {
    push_one(flds, prts, n);
  }
}

void
cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    psc_push_particles_push_a_yz(psc_mparticles_get_patch(mprts, p),
				 psc_mfields_get_patch(mflds, p));
  }
}

