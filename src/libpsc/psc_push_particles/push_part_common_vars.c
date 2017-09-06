
#if (DIM & DIM_X)
  particle_real_t s0x[N_RHO] = {}, s1x[N_RHO];
#endif
#if (DIM & DIM_Y)
  particle_real_t s0y[N_RHO] = {}, s1y[N_RHO];
#endif
#if (DIM & DIM_Z)
  particle_real_t s0z[N_RHO] = {}, s1z[N_RHO];
#endif

  c_prm_set(ppsc);
