
// ----------------------------------------------------------------------
// calc_v

CUDA_DEVICE static inline void
calc_v(particle_real_t v[3], const particle_real_t p[3])
{
#ifdef __CUDACC__
  particle_real_t root = rsqrt(1.f + sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
#else
  particle_real_t root = 1.f 
    / particle_real_sqrt(1.f + sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
#endif
  for (int d = 0; d < 3; d++) {
    v[d] = p[d] * root;
  }
}

// ----------------------------------------------------------------------
// push_xi

CUDA_DEVICE static inline void
push_xi(particle_t *prt, particle_real_t vxi[3], particle_real_t dt)
{
#if DIM == DIM_YZ

#if PSC_PARTICLES_AS_CUDA2 || PSC_PARTICLES_AS_ACC
  prt->xi[1] += vxi[1] * dt;
  prt->xi[2] += vxi[2] * dt;
#else
  prt->yi += vxi[1] * dt;
  prt->zi += vxi[2] * dt;
#endif

#elif DIM == DIM_XYZ

#if PSC_PARTICLES_AS_CUDA2 || PSC_PARTICLES_AS_ACC
  prt->xi[0] += vxi[0] * dt;
  prt->xi[1] += vxi[1] * dt;
  prt->xi[2] += vxi[2] * dt;
#else
  prt->xi += vxi[0] * dt;
  prt->yi += vxi[1] * dt;
  prt->zi += vxi[2] * dt;
#endif

#endif
}

// ----------------------------------------------------------------------
// push_pxi

CUDA_DEVICE static inline void
push_pxi(particle_t *prt, particle_real_t exq, particle_real_t eyq, particle_real_t ezq,
	 particle_real_t hxq, particle_real_t hyq, particle_real_t hzq, particle_real_t dq)
{
#if PSC_PARTICLES_AS_CUDA2 || PSC_PARTICLES_AS_ACC
  particle_real_t pxm = prt->pxi[0] + dq*exq;
  particle_real_t pym = prt->pxi[1] + dq*eyq;
  particle_real_t pzm = prt->pxi[2] + dq*ezq;
#else
  particle_real_t pxm = prt->pxi + dq*exq;
  particle_real_t pym = prt->pyi + dq*eyq;
  particle_real_t pzm = prt->pzi + dq*ezq;
#endif
  
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
  
#if PSC_PARTICLES_AS_CUDA2 || PSC_PARTICLES_AS_ACC
  prt->pxi[0] = pxp + dq * exq;
  prt->pxi[1] = pyp + dq * eyq;
  prt->pxi[2] = pzp + dq * ezq;
#else
  prt->pxi = pxp + dq * exq;
  prt->pyi = pyp + dq * eyq;
  prt->pzi = pzp + dq * ezq;
#endif
}

// ----------------------------------------------------------------------
// find_idx_off_1st_rel

CUDA_DEVICE static inline void
find_idx_off_1st_rel(particle_real_t xi[3], int lg[3], particle_real_t og[3], particle_real_t shift)
{
  for (int d = 0; d < 3; d++) {
    particle_real_t pos = xi[d] * c_prm.dxi[d] + shift;
    lg[d] = particle_real_fint(pos);
    og[d] = pos - lg[d];
  }
}

// ----------------------------------------------------------------------
// find_idx_off_pos_1st_rel

CUDA_DEVICE static inline void
find_idx_off_pos_1st_rel(particle_real_t xi[3], int lg[3], particle_real_t og[3],
			 particle_real_t pos[3], particle_real_t shift)
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * c_prm.dxi[d] + shift;
    lg[d] = particle_real_fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

