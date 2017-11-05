
#ifndef PUSH_DIM
#define PUSH_DIM DIM
#endif

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
// push_x

CUDA_DEVICE static inline void
push_x(particle_real_t x[3], const particle_real_t v[3], particle_real_t dt)
{
  IF_DIM_X( x[0] += v[0] * dt; );
  IF_DIM_Y( x[1] += v[1] * dt; );
  IF_DIM_Z( x[2] += v[2] * dt; );
}

// ----------------------------------------------------------------------
// push_p

CUDA_DEVICE static inline void
push_p(particle_real_t p[3], const particle_real_t E[3], const particle_real_t H[3],
       particle_real_t dq)
{
  particle_real_t pxm = p[0] + dq * E[0];
  particle_real_t pym = p[1] + dq * E[1];
  particle_real_t pzm = p[2] + dq * E[2];

  particle_real_t root = dq / particle_real_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
  particle_real_t taux = H[0] * root;
  particle_real_t tauy = H[1] * root;
  particle_real_t tauz = H[2] * root;
  
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
  
  p[0] = pxp + dq * E[0];
  p[1] = pyp + dq * E[1];
  p[2] = pzp + dq * E[2];
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

