
// ----------------------------------------------------------------------
// calc_vxi

static inline void
calc_vxi(particle_real_t vxi[3], particle_t *prt)
{
#ifdef PSC_PARTICLES_AS_CUDA2
  particle_real_t root = 1.f 
    / particle_real_sqrt(1.f + sqr(prt->pxi4.x) + sqr(prt->pxi4.y) + sqr(prt->pxi4.z));
  vxi[0] = prt->pxi4.x * root;
  vxi[1] = prt->pxi4.y * root;
  vxi[2] = prt->pxi4.z * root;
#else
  particle_real_t root = 1.f 
    / particle_real_sqrt(1.f + sqr(prt->pxi) + sqr(prt->pyi) + sqr(prt->pzi));
  vxi[0] = prt->pxi * root;
  vxi[1] = prt->pyi * root;
  vxi[2] = prt->pzi * root;
#endif
}

// ----------------------------------------------------------------------
// push_xi

#ifdef PSC_PARTICLES_AS_CUDA2

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

#else

#if DIM == DIM_YZ

static inline void
push_xi(particle_t *part, particle_real_t vxi[3], particle_real_t dt)
{
  part->yi += vxi[1] * dt;
  part->zi += vxi[2] * dt;
}

#endif

#endif

// ----------------------------------------------------------------------
// push_pxi

static inline void
push_pxi(particle_t *prt, particle_real_t exq, particle_real_t eyq, particle_real_t ezq,
	 particle_real_t hxq, particle_real_t hyq, particle_real_t hzq, particle_real_t dq)
{
#ifdef PSC_PARTICLES_AS_CUDA2
  particle_real_t pxm = prt->pxi4.x + dq*exq;
  particle_real_t pym = prt->pxi4.y + dq*eyq;
  particle_real_t pzm = prt->pxi4.z + dq*ezq;
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
  
#ifdef PSC_PARTICLES_AS_CUDA2
  prt->pxi4.x = pxp + dq * exq;
  prt->pxi4.y = pyp + dq * eyq;
  prt->pxi4.z = pzp + dq * ezq;
#else
  prt->pxi = pxp + dq * exq;
  prt->pyi = pyp + dq * eyq;
  prt->pzi = pzp + dq * ezq;
#endif
}

// ----------------------------------------------------------------------
// find_idx_off_1st_rel

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

// ----------------------------------------------------------------------
// find_idx_off_pos_1st_rel

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

