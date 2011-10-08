
static inline void
calc_vxi(creal vxi[3], particle_t *part)
{
  creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline void
push_xi(particle_t *part, creal vxi[3], creal dt)
{
  part->yi += vxi[1] * dt;
  part->zi += vxi[2] * dt;
}

static inline void
push_pxi(particle_t *part, creal exq, creal eyq, creal ezq,
	 creal hxq, creal hyq, creal hzq, creal dqs)
{
  creal dq = dqs * part->qni / part->mni;
  creal pxm = part->pxi + dq*exq;
  creal pym = part->pyi + dq*eyq;
  creal pzm = part->pzi + dq*ezq;
  
  creal root = dq / creal_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
  creal taux = hxq*root;
  creal tauy = hyq*root;
  creal tauz = hzq*root;
  
  creal tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
  creal pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
	       (2.f*taux*tauy+2.f*tauz)*pym + 
	       (2.f*taux*tauz-2.f*tauy)*pzm)*tau;
  creal pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
	       (1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
	       (2.f*tauy*tauz+2.f*taux)*pzm)*tau;
  creal pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
	       (2.f*tauy*tauz-2.f*taux)*pym +
	       (1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
  
  part->pxi = pxp + dq * exq;
  part->pyi = pyp + dq * eyq;
  part->pzi = pzp + dq * ezq;
}

static inline void
find_idx_off_1st(creal xi[3], int lg[3], creal og[3], creal shift,
		 double xb[3], double dxi[3])
{
  for (int d = 0; d < 3; d++) {
    creal pos = (xi[d] - xb[d]) * dxi[d] + shift;
    lg[d] = fint(pos);
    og[d] = pos - lg[d];
  }
}

static inline void
find_idx_off_pos_1st(creal xi[3], int lg[3], creal og[3], creal pos[3], creal shift,
		     double xb[3], double dxi[3])
{
  for (int d = 0; d < 3; d++) {
    pos[d] = (xi[d] - xb[d]) * dxi[d] + shift;
    lg[d] = fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

#define INTERPOLATE_SETUP_1ST			\
  creal g0y = 1.f - og[1];			\
  creal g0z = 1.f - og[2];			\
  creal g1y = og[1];				\
  creal g1z = og[2];				\
						\
  creal h0y = 1.f - oh[1];			\
  creal h0z = 1.f - oh[2];			\
  creal h1y = oh[1];				\
  creal h1z = oh[2]


#define INTERPOLATE_FIELD_1ST(m, gy, gz)				\
    (gz##0z*(gy##0y*F3(pf, m, 0,l##gy[1]  ,l##gz[2]  ) +			\
	     gy##1y*F3(pf, m, 0,l##gy[1]+1,l##gz[2]  )) +			\
     gz##1z*(gy##0y*F3(pf, m, 0,l##gy[1]  ,l##gz[2]+1) +			\
	     gy##1y*F3(pf, m, 0,l##gy[1]+1,l##gz[2]+1)))
