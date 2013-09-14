
#define DEPOSIT_TO_GRID_1ST_CC(part, pf, m, val) do {			\
    particle_real_t xi[3];						\
    particle_get_relative_pos(part, patch->xb, xi);			\
    particle_real_t u = xi[0] * dxi - .5;				\
    particle_real_t v = xi[1] * dyi - .5;				\
    particle_real_t w = xi[2] * dzi - .5;				\
    int jx = particle_real_fint(u);					\
    int jy = particle_real_fint(v);					\
    int jz = particle_real_fint(w);					\
    particle_real_t h1 = u - jx;					\
    particle_real_t h2 = v - jy;					\
    particle_real_t h3 = w - jz;					\
    									\
    particle_real_t g0x = 1.f - h1;					\
    particle_real_t g0y = 1.f - h2;					\
    particle_real_t g0z = 1.f - h3;					\
    particle_real_t g1x = h1;						\
    particle_real_t g1y = h2;						\
    particle_real_t g1z = h3;						\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (ppsc->domain.gdims[0] == 1) {					\
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;				\
    }									\
    if (ppsc->domain.gdims[1] == 1) {					\
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;				\
    }									\
    if (ppsc->domain.gdims[2] == 1) {					\
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;				\
    }									\
    									\
    assert(jx >= -1 && jx < patch->ldims[0]);				\
    assert(jy >= -1 && jy < patch->ldims[1]);				\
    assert(jz >= -1 && jz < patch->ldims[2]);				\
    									\
    particle_real_t fnq = particle_wni(part) * fnqs;			\
									\
    F3(pf, m, jx    ,jy  ,jz  ) += fnq*g0x*g0y*g0z * (val);		\
    F3(pf, m, jx+jxd,jy  ,jz  ) += fnq*g1x*g0y*g0z * (val);		\
    F3(pf, m, jx    ,jy+1,jz  ) += fnq*g0x*g1y*g0z * (val);		\
    F3(pf, m, jx+jxd,jy+1,jz  ) += fnq*g1x*g1y*g0z * (val);		\
    F3(pf, m, jx    ,jy  ,jz+1) += fnq*g0x*g0y*g1z * (val);		\
    F3(pf, m, jx+jxd,jy  ,jz+1) += fnq*g1x*g0y*g1z * (val);		\
    F3(pf, m, jx    ,jy+1,jz+1) += fnq*g0x*g1y*g1z * (val);		\
    F3(pf, m, jx+jxd,jy+1,jz+1) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

#define DEPOSIT_TO_GRID_1ST_NC(part, pf, m, val) do {			\
    particle_real_t *xi = &part->xi; /* don't shift back in time */	\
    particle_real_t u = xi[0] * dxi;					\
    particle_real_t v = xi[1] * dyi;					\
    particle_real_t w = xi[2] * dzi;					\
    int jx = particle_real_fint(u);					\
    int jy = particle_real_fint(v);					\
    int jz = particle_real_fint(w);					\
    particle_real_t h1 = u - jx;					\
    particle_real_t h2 = v - jy;					\
    particle_real_t h3 = w - jz;					\
    									\
    particle_real_t g0x = 1.f - h1;					\
    particle_real_t g0y = 1.f - h2;					\
    particle_real_t g0z = 1.f - h3;					\
    particle_real_t g1x = h1;						\
    particle_real_t g1y = h2;						\
    particle_real_t g1z = h3;						\
    									\
    if (ppsc->domain.gdims[0] == 1) {					\
      jx = 0; g0x = 1.; g1x = 0.;					\
    }									\
    if (ppsc->domain.gdims[1] == 1) {					\
      jy = 0; g0y = 1.; g1y = 0.;					\
    }									\
    if (ppsc->domain.gdims[2] == 1) {					\
      jz = 0; g0z = 1.; g1z = 0.;					\
    }									\
    									\
    assert(jx >= -1 && jx < patch->ldims[0]);				\
    assert(jy >= -1 && jy < patch->ldims[1]);				\
    assert(jz >= -1 && jz < patch->ldims[2]);				\
    									\
    particle_real_t fnq = particle_wni(part) * fnqs;			\
									\
    F3(pf, m, jx  ,jy  ,jz  ) += fnq*g0x*g0y*g0z * (val);		\
    F3(pf, m, jx+1,jy  ,jz  ) += fnq*g1x*g0y*g0z * (val);		\
    F3(pf, m, jx  ,jy+1,jz  ) += fnq*g0x*g1y*g0z * (val);		\
    F3(pf, m, jx+1,jy+1,jz  ) += fnq*g1x*g1y*g0z * (val);		\
    F3(pf, m, jx  ,jy  ,jz+1) += fnq*g0x*g0y*g1z * (val);		\
    F3(pf, m, jx+1,jy  ,jz+1) += fnq*g1x*g0y*g1z * (val);		\
    F3(pf, m, jx  ,jy+1,jz+1) += fnq*g0x*g1y*g1z * (val);		\
    F3(pf, m, jx+1,jy+1,jz+1) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

// FIXME, this function exists about 100x all over the place, should
// be consolidated

static inline void
particle_calc_vxi(particle_t *part, particle_real_t vxi[3])
{
  particle_real_t root =
    1.f / particle_real_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline int
particle_kind_ei(particle_t *part)
{
  if (particle_qni_div_mni(part) < 0.) {
    return 0;
  } else if (particle_qni_div_mni(part) > 0.) {
    return 1;
  } else {
    assert(0);
  }
}

