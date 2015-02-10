
#include "psc_debug.h"

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
calc_3d_dx1(particle_real_t dx1[3], particle_real_t x[3], particle_real_t dx[3], int off[3])
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

static inline void
curr_3d_vb_cell(struct psc_fields *pf, int i[3], particle_real_t x[3], particle_real_t dx[3],
		particle_real_t fnq[3], particle_real_t dxt[3], int off[3])
{
  particle_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];
  particle_real_t xa[3] = { 0.,
			    x[1] + .5f * dx[1],
			    x[2] + .5f * dx[2], };
  F3_CURR(pf, JXI, 0,i[1]  ,i[2]  ) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h);
  F3_CURR(pf, JXI, 0,i[1]+1,i[2]  ) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h);
  F3_CURR(pf, JXI, 0,i[1]  ,i[2]+1) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h);
  F3_CURR(pf, JXI, 0,i[1]+1,i[2]+1) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h);

  F3_CURR(pf, JYI, 0,i[1]  ,i[2]  ) += fnq[1] * dx[1] * (.5f - xa[2]);
  F3_CURR(pf, JYI, 0,i[1]  ,i[2]+1) += fnq[1] * dx[1] * (.5f + xa[2]);
  F3_CURR(pf, JZI, 0,i[1]  ,i[2]  ) += fnq[2] * dx[2] * (.5f - xa[1]);
  F3_CURR(pf, JZI, 0,i[1]+1,i[2]  ) += fnq[2] * dx[2] * (.5f + xa[1]);
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

#define CALC_JX_2D(pf, part, vxi)					\
  do {									\
    int lf[3];								\
    particle_real_t of[3];						\
    find_idx_off_1st_rel(&part->xi, lf, of, 0.f, prm.dxi);		\
    									\
    particle_real_t fnqx = vxi[0] * particle_wni(part) * prm.fnqx_kind[part->kind]; \
    F3_CURR(pf, JXI, 0,lf[1]  ,lf[2]  ) += (1.f - of[1]) * (1.f - of[2]) * fnqx; \
    F3_CURR(pf, JXI, 0,lf[1]+1,lf[2]  ) += (      of[1]) * (1.f - of[2]) * fnqx; \
    F3_CURR(pf, JXI, 0,lf[1]  ,lf[2]+1) += (1.f - of[1]) * (      of[2]) * fnqx; \
    F3_CURR(pf, JXI, 0,lf[1]+1,lf[2]+1) += (      of[1]) * (      of[2]) * fnqx; \
  } while (0)

#define CALC_JYZ_2D(pf, xm, xp)						\
  do {									\
    int i[2] = { lg[1], lg[2] };					\
    int idiff[2] = { lf[1] - lg[1], lf[2] - lg[2] };			\
    particle_real_t dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };		\
    particle_real_t x[2] = { xm[1] - (i[0] + .5f), xm[2] - (i[1] + .5f) }; \
    									\
    particle_real_t dx1[2];						\
    int off[2];								\
    int first_dir, second_dir = -1;					\
    /* FIXME, make sure we never div-by-zero? */			\
    if (idiff[0] == 0 && idiff[1] == 0) {				\
      first_dir = -1;							\
    } else if (idiff[0] == 0) {						\
      first_dir = 1;							\
    } else if (idiff[1] == 0) {						\
      first_dir = 0;							\
    } else {								\
      dx1[0] = .5f * idiff[0] - x[0];					\
      if (dx[0] == 0.f) {						\
	dx1[1] = 0.f;							\
    } else {								\
	dx1[1] = dx[1] / dx[0] * dx1[0];				\
      }									\
      if (particle_real_abs(x[1] + dx1[1]) > .5f) {			\
	first_dir = 1;							\
      } else {								\
	first_dir = 0;							\
      }									\
      second_dir = 1 - first_dir;					\
    }									\
    									\
    particle_real_t fnq[2] = { particle_wni(prt) * prm.fnqy_kind[prt->kind], \
			       particle_wni(prt) * prm.fnqz_kind[prt->kind] }; \
    									\
    if (first_dir >= 0) {						\
      off[1-first_dir] = 0;						\
      off[first_dir] = idiff[first_dir];				\
      calc_dx1(dx1, x, dx, off);					\
      curr_2d_vb_cell(pf, i, x, dx1, fnq, dx, off);			\
    }									\
    									\
    if (second_dir >= 0) {						\
      off[first_dir] = 0;						\
      off[second_dir] = idiff[second_dir];				\
      calc_dx1(dx1, x, dx, off);					\
      curr_2d_vb_cell(pf, i, x, dx1, fnq, dx, off);			\
    }									\
    									\
    curr_2d_vb_cell(pf, i, x, dx, fnq, NULL, NULL);			\
  } while (0)

#define CALC_JXYZ_3D(pf, xm, xp)					\
  do {									\
    int i[3] = { 0, lg[1], lg[2] };					\
    int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };			\
    particle_real_t dx[3] = { 0., xp[1] - xm[1], xp[2] - xm[2] };	\
    particle_real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) }; \
    									\
    particle_real_t dx1[3];						\
    int off[3];								\
    int first_dir, second_dir = -1;					\
    /* FIXME, make sure we never div-by-zero? */			\
    if (idiff[1] == 0 && idiff[2] == 0) {				\
      first_dir = -1;							\
    } else if (idiff[1] == 0) {						\
      first_dir = 2;							\
    } else if (idiff[2] == 0) {						\
      first_dir = 1;							\
    } else {								\
      dx1[1] = .5f * idiff[1] - x[1];					\
      if (dx[1] == 0.f) {						\
	dx1[2] = 0.f;							\
      } else {								\
	dx1[2] = dx[2] / dx[1] * dx1[1];				\
      }									\
      if (particle_real_abs(x[2] + dx1[2]) > .5f) {			\
	first_dir = 2;							\
      } else {								\
	first_dir = 1;							\
      }									\
      second_dir = 3 - first_dir;					\
    }									\
									\
    particle_real_t fnq[3] = { particle_wni(prt) * prm.fnqx_kind[prt->kind], \
			       particle_wni(prt) * prm.fnqy_kind[prt->kind], \
			       particle_wni(prt) * prm.fnqz_kind[prt->kind] }; \
    dx[0] = vxi[0] * prm.dt * prm.dxi[0];				\
    									\
    if (first_dir >= 0) {						\
      off[3 - first_dir] = 0;						\
      off[first_dir] = idiff[first_dir];				\
      calc_3d_dx1(dx1, x, dx, off);					\
      curr_3d_vb_cell(pf, i, x, dx1, fnq, dx, off);			\
    }									\
									\
    if (second_dir >= 0) {						\
      off[first_dir] = 0;						\
      off[second_dir] = idiff[second_dir];				\
      calc_3d_dx1(dx1, x, dx, off);					\
      curr_3d_vb_cell(pf, i, x, dx1, fnq, dx, off);			\
    }									\
									\
    curr_3d_vb_cell(pf, i, x, dx, fnq, NULL, NULL);			\
  } while (0)

#endif

#ifdef F3_CACHE

#include "psc_fields_single.h"

static struct psc_fields *
cache_fields_from_em(fields_t *pf)
{
  struct psc_fields *fld = psc_fields_create(psc_fields_comm(pf));
  psc_fields_set_type(fld, F3_CACHE_TYPE);
  // FIXME, can do -1 .. 1?
  psc_fields_set_param_int3(fld, "ib", (int[3]) { 0, -2, -2 });
  psc_fields_set_param_int3(fld, "im", (int[3]) { 1,
	pf->im[1] + 2 * pf->ib[1] + 4, pf->im[2] + 2 * pf->ib[2] + 4});
  psc_fields_set_param_int(fld, "nr_comp", 9); // JX .. HZ
  psc_fields_set_param_int(fld, "p", pf->p);
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

static void _mrc_unused
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

