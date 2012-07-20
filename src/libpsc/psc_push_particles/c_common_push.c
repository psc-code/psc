
static inline void
calc_vxi(particle_real_t vxi[3], particle_t *part)
{
  particle_real_t root = 1.f 
    / particle_real_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline void
push_xi(particle_t *part, particle_real_t vxi[3], particle_real_t dt)
{
  part->yi += vxi[1] * dt;
  part->zi += vxi[2] * dt;
}

static inline void
push_pxi(particle_t *part, particle_real_t exq, particle_real_t eyq, particle_real_t ezq,
	 particle_real_t hxq, particle_real_t hyq, particle_real_t hzq, particle_real_t dq)
{
  particle_real_t pxm = part->pxi + dq*exq;
  particle_real_t pym = part->pyi + dq*eyq;
  particle_real_t pzm = part->pzi + dq*ezq;
  
  particle_real_t root = dq / particle_real_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
  particle_real_t taux = hxq*root;
  particle_real_t tauy = hyq*root;
  particle_real_t tauz = hzq*root;
  
  particle_real_t tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
  particle_real_t pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
	       (2.f*taux*tauy+2.f*tauz)*pym + 
	       (2.f*taux*tauz-2.f*tauy)*pzm)*tau;
  particle_real_t pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
	       (1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
	       (2.f*tauy*tauz+2.f*taux)*pzm)*tau;
  particle_real_t pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
	       (2.f*tauy*tauz-2.f*taux)*pym +
	       (1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
  
  part->pxi = pxp + dq * exq;
  part->pyi = pyp + dq * eyq;
  part->pzi = pzp + dq * ezq;
}

static inline void
find_idx_off_1st_rel(particle_real_t xi[3], int lg[3], particle_real_t og[3], particle_real_t shift,
		     particle_real_t dxi[3])
{
  for (int d = 0; d < 3; d++) {
    particle_real_t pos = xi[d] * dxi[d] + shift;
    lg[d] = particle_real_fint(pos);
    og[d] = pos - lg[d];
  }
}

static inline void
find_idx_off_pos_1st_rel(particle_real_t xi[3], int lg[3], particle_real_t og[3],
			 particle_real_t pos[3], particle_real_t shift,
			 particle_real_t dxi[3])
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * dxi[d] + shift;
    lg[d] = particle_real_fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

#define INTERPOLATE_SETUP_1ST				\
  particle_real_t g0y = 1.f - og[1];			\
  particle_real_t g0z = 1.f - og[2];			\
  particle_real_t g1y = og[1];				\
  particle_real_t g1z = og[2];				\
  							\
  particle_real_t h0y = 1.f - oh[1];			\
  particle_real_t h0z = 1.f - oh[2];			\
  particle_real_t h1y = oh[1];				\
  particle_real_t h1z = oh[2]


#define INTERPOLATE_FIELD_1ST(m, gy, gz)				\
  (gz##0z*(gy##0y*F3(pf, m, 0,l##gy[1]  ,l##gz[2]  ) +			\
	   gy##1y*F3(pf, m, 0,l##gy[1]+1,l##gz[2]  )) +			\
   gz##1z*(gy##0y*F3(pf, m, 0,l##gy[1]  ,l##gz[2]+1) +			\
	   gy##1y*F3(pf, m, 0,l##gy[1]+1,l##gz[2]+1)))

#define INTERPOLATE_FIELD_1ST_CACHE(m, gy, gz)				\
    (gz##0z*(gy##0y*F3_CACHE(pf, m, 0,l##gy[1]  ,l##gz[2]  ) +		\
	     gy##1y*F3_CACHE(pf, m, 0,l##gy[1]+1,l##gz[2]  )) +		\
     gz##1z*(gy##0y*F3_CACHE(pf, m, 0,l##gy[1]  ,l##gz[2]+1) +		\
	     gy##1y*F3_CACHE(pf, m, 0,l##gy[1]+1,l##gz[2]+1)))

#ifdef F3_CURR

// ======================================================================
// current 1vb (yz)

static inline void
calc_dx1(particle_real_t dx1[2], particle_real_t x[2], particle_real_t dx[2], int off[2])
{
  if (off[1] == 0) {
    dx1[0] = .5f * off[0] - x[0];
   if (dx[0] == 0.f) {
    dx1[1] = 0.f;
   } else {
    dx1[1] = dx[1] / dx[0] * dx1[0];
   }
  } else {
    dx1[1] = .5f * off[1] - x[1];
   if (dx[1] == 0.f) {
    dx1[0] = 0.f;
   } else {
    dx1[0] = dx[0] / dx[1] * dx1[1];
   }
  }
}

static inline void
curr_2d_vb_cell(struct psc_fields *pf, int i[2], particle_real_t x[2], particle_real_t dx[2],
		particle_real_t fnq[2], particle_real_t dxt[2], int off[2])
{
  F3_CURR(pf, JYI, 0,i[0],i[1]  ) += fnq[0] * dx[0] * (.5f - x[1] - .5f * dx[1]);
  F3_CURR(pf, JYI, 0,i[0],i[1]+1) += fnq[0] * dx[0] * (.5f + x[1] + .5f * dx[1]);
  F3_CURR(pf, JZI, 0,i[0],i[1]  ) += fnq[1] * dx[1] * (.5f - x[0] - .5f * dx[0]);
  F3_CURR(pf, JZI, 0,i[0]+1,i[1]) += fnq[1] * dx[1] * (.5f + x[0] + .5f * dx[0]);
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    x[0] += dx[0] - off[0];
    x[1] += dx[1] - off[1];
    i[0] += off[0];
    i[1] += off[1];
  }
}

#endif

#ifdef F3_CACHE

#include "psc_fields_single.h"

static struct psc_fields *
cache_fields_from_em(fields_t *pf)
{
  struct psc_fields *fld = psc_fields_create(psc_fields_comm(pf));
  // FIXME, can do -1 .. 1?
  psc_fields_set_param_int3(fld, "ib", (int[3]) { 0, -2, -2 });
  psc_fields_set_param_int3(fld, "im", (int[3]) { 1,
	pf->im[1] + 2 * pf->ib[1] + 4, pf->im[2] + 2 * pf->ib[2] + 4});
  psc_fields_set_param_int(fld, "nr_comp", 9); // JX .. HZ
  psc_fields_setup(fld);
  for (int iz = fld->ib[2]; iz < fld->ib[2] + fld->im[2]; iz++) {
    for (int iy = fld->ib[1]; iy < fld->ib[1] + fld->im[1]; iy++) {
      F3_CACHE(fld, EX, 0,iy,iz) = F3(pf, EX, 0,iy,iz);
      F3_CACHE(fld, EY, 0,iy,iz) = F3(pf, EY, 0,iy,iz);
      F3_CACHE(fld, EZ, 0,iy,iz) = F3(pf, EZ, 0,iy,iz);
      F3_CACHE(fld, HX, 0,iy,iz) = F3(pf, HX, 0,iy,iz);
      F3_CACHE(fld, HY, 0,iy,iz) = F3(pf, HY, 0,iy,iz);
      F3_CACHE(fld, HZ, 0,iy,iz) = F3(pf, HZ, 0,iy,iz);
    }
  }
  return fld;
}

static void __unused
cache_fields_to_j(struct psc_fields *fld, fields_t *pf)
{
  for (int iz = fld->ib[2]; iz < fld->ib[2] + fld->im[2]; iz++) {
    for (int iy = fld->ib[1]; iy < fld->ib[1] + fld->im[1]; iy++) {
      F3(pf, JXI, 0,iy,iz) += F3_CACHE(fld, JXI, 0,iy,iz);
      F3(pf, JYI, 0,iy,iz) += F3_CACHE(fld, JYI, 0,iy,iz);
      F3(pf, JZI, 0,iy,iz) += F3_CACHE(fld, JZI, 0,iy,iz);
    }
  }
}

#endif

