
  particle_real_t dt = ppsc->dt;
  particle_real_t dqs = .5f * ppsc->coeff.eta * dt;
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;

#if (DIM & DIM_X)
  particle_real_t s0x[N_RHO] = {}, s1x[N_RHO];
  xl = .5f * dt;
  particle_real_t fnqxs = ppsc->patch[p].dx[0] * fnqs / dt;
  particle_real_t dxi = 1.f / ppsc->patch[p].dx[0];
#endif
#if (DIM & DIM_Y)
  particle_real_t s0y[N_RHO] = {}, s1y[N_RHO];
  yl = .5f * dt;
  particle_real_t fnqys = ppsc->patch[p].dx[1] * fnqs / dt;
  particle_real_t dyi = 1.f / ppsc->patch[p].dx[1];
#endif
#if (DIM & DIM_Z)
  particle_real_t s0z[N_RHO] = {}, s1z[N_RHO];
  zl = .5f * dt;
  particle_real_t fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;
  particle_real_t dzi = 1.f / ppsc->patch[p].dx[2];
#endif
