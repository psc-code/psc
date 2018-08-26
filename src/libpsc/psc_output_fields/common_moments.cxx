
#pragma once

#include "fields_item.hxx"

#include <string>

#define DEPOSIT_TO_GRID_1ST_CC(part, flds, m, val) do {			\
    using Fields = Fields3d<fields_t>;					\
    Fields F(flds);							\
    real_t *xi = part->x; /* don't shift back in time */		\
    real_t u = xi[0] * dxi - .5;					\
    real_t v = xi[1] * dyi - .5;					\
    real_t w = xi[2] * dzi - .5;					\
    int jx = fint(u);							\
    int jy = fint(v);							\
    int jz = fint(w);							\
    real_t h1 = u - jx;							\
    real_t h2 = v - jy;							\
    real_t h3 = w - jz;							\
    									\
    real_t g0x = 1.f - h1;						\
    real_t g0y = 1.f - h2;						\
    real_t g0z = 1.f - h3;						\
    real_t g1x = h1;							\
    real_t g1y = h2;							\
    real_t g1z = h3;							\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (flds.grid().isInvar(0)) {					\
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;				\
    }									\
    if (flds.grid().isInvar(1)) {					\
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;				\
    }									\
    if (flds.grid().isInvar(2)) {					\
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;				\
    }									\
    									\
    assert(jx >= -1 && jx < grid.ldims[0]);				\
    assert(jy >= -1 && jy < grid.ldims[1]);				\
    assert(jz >= -1 && jz < grid.ldims[2]);				\
    									\
    real_t fnq = prts.prt_wni(*part) * fnqs;				\
    									\
    F(m, jx    ,jy    ,jz    ) += fnq*g0x*g0y*g0z * (val);		\
    F(m, jx+jxd,jy    ,jz    ) += fnq*g1x*g0y*g0z * (val);		\
    F(m, jx    ,jy+jyd,jz    ) += fnq*g0x*g1y*g0z * (val);		\
    F(m, jx+jxd,jy+jyd,jz    ) += fnq*g1x*g1y*g0z * (val);		\
    F(m, jx    ,jy    ,jz+jzd) += fnq*g0x*g0y*g1z * (val);		\
    F(m, jx+jxd,jy    ,jz+jzd) += fnq*g1x*g0y*g1z * (val);		\
    F(m, jx    ,jy+jyd,jz+jzd) += fnq*g0x*g1y*g1z * (val);		\
    F(m, jx+jxd,jy+jyd,jz+jzd) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

#define DEPOSIT_TO_GRID_1ST_NC(part, flds, m, val) do {			\
    using Fields = Fields3d<fields_t>;					\
    Fields F(flds);							\
    real_t *xi = part->x; /* don't shift back in time */		\
    real_t u = xi[0] * dxi;						\
    real_t v = xi[1] * dyi;						\
    real_t w = xi[2] * dzi;						\
    int jx = fint(u);							\
    int jy = fint(v);							\
    int jz = fint(w);							\
    real_t h1 = u - jx;							\
    real_t h2 = v - jy;							\
    real_t h3 = w - jz;							\
    									\
    real_t g0x = 1.f - h1;						\
    real_t g0y = 1.f - h2;						\
    real_t g0z = 1.f - h3;						\
    real_t g1x = h1;							\
    real_t g1y = h2;							\
    real_t g1z = h3;							\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (grid.isInvar(0)) {						\
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;				\
    }									\
    if (grid.isInvar(1)) {					\
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;				\
    }									\
    if (grid.isInvar(2)) {					\
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;				\
    }									\
    									\
    assert(jx >= -1 && jx < grid.ldims[0]);				\
    assert(jy >= -1 && jy < grid.ldims[1]);				\
    assert(jz >= -1 && jz < grid.ldims[2]);				\
    									\
    real_t fnq = prts.prt_wni(*part) * fnqs;				\
									\
    F(m, jx    ,jy    ,jz    ) += fnq*g0x*g0y*g0z * (val);		\
    F(m, jx+jxd,jy    ,jz    ) += fnq*g1x*g0y*g0z * (val);		\
    F(m, jx    ,jy+jyd,jz    ) += fnq*g0x*g1y*g0z * (val);		\
    F(m, jx+jxd,jy+jyd,jz    ) += fnq*g1x*g1y*g0z * (val);		\
    F(m, jx    ,jy    ,jz+jzd) += fnq*g0x*g0y*g1z * (val);		\
    F(m, jx+jxd,jy    ,jz+jzd) += fnq*g1x*g0y*g1z * (val);		\
    F(m, jx    ,jy+jyd,jz+jzd) += fnq*g0x*g1y*g1z * (val);		\
    F(m, jx+jxd,jy+jyd,jz+jzd) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

#define DEPOSIT_TO_GRID_2ND_NC(part, flds, m, val) do {			\
    using Fields = Fields3d<fields_t>;					\
    Fields F(flds);							\
    real_t *xi = part->x; /* don't shift back in time */		\
    real_t u = xi[0] * dxi;						\
    real_t v = xi[1] * dyi;						\
    real_t w = xi[2] * dzi;						\
    int jx = nint(u);							\
    int jy = nint(v);							\
    int jz = nint(w);							\
    real_t h1 = jx - u;							\
    real_t h2 = jy - v;							\
    real_t h3 = jz - w;							\
    									\
    real_t gmx = .5f*(.5f+h1)*(.5f+h1);					\
    real_t gmy = .5f*(.5f+h2)*(.5f+h2);					\
    real_t gmz = .5f*(.5f+h3)*(.5f+h3);					\
    real_t g0x = .75f-h1*h1;						\
    real_t g0y = .75f-h2*h2;						\
    real_t g0z = .75f-h3*h3;						\
    real_t g1x = .5f*(.5f-h1)*(.5f-h1);					\
    real_t g1y = .5f*(.5f-h2)*(.5f-h2);					\
    real_t g1z = .5f*(.5f-h3)*(.5f-h3);					\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (grid.isInvar(0)) {						\
      jx = 0; g0x = 1.; g1x = 0.; gmx = 0.; jxd = 0;			\
    }									\
    if (grid.isInvar(1)) {						\
      jy = 0; g0y = 1.; g1y = 0.; gmy = 0.; jyd = 0;			\
    }									\
    if (grid.isInvar(2)) {						\
      jz = 0; g0z = 1.; g1z = 0.; gmz = 0.; jzd = 0;			\
    }									\
    									\
    assert(jx >= 0 && jx <= grid.ldims[0]);				\
    assert(jy >= 0 && jy <= grid.ldims[1]);				\
    assert(jz >= 0 && jz <= grid.ldims[2]);				\
    									\
    real_t fnq = prts.prt_wni(*part) * fnqs;				\
					     				\
    F(m, jx-jxd,jy-jyd,jz-jzd) += fnq*gmx*gmy*gmz * (val);		\
    F(m, jx    ,jy-jyd,jz-jzd) += fnq*g0x*gmy*gmz * (val);		\
    F(m, jx+jxd,jy-jyd,jz-jzd) += fnq*g1x*gmy*gmz * (val);		\
    F(m, jx-jxd,jy    ,jz-jzd) += fnq*gmx*g0y*gmz * (val);		\
    F(m, jx    ,jy    ,jz-jzd) += fnq*g0x*g0y*gmz * (val);		\
    F(m, jx+jxd,jy    ,jz-jzd) += fnq*g1x*g0y*gmz * (val);		\
    F(m, jx-jxd,jy+jyd,jz-jzd) += fnq*gmx*g1y*gmz * (val);		\
    F(m, jx    ,jy+jyd,jz-jzd) += fnq*g0x*g1y*gmz * (val);		\
    F(m, jx+jxd,jy+jyd,jz-jzd) += fnq*g1x*g1y*gmz * (val);		\
    F(m, jx-jxd,jy-jyd,jz    ) += fnq*gmx*gmy*g0z * (val);		\
    F(m, jx    ,jy-jyd,jz    ) += fnq*g0x*gmy*g0z * (val);		\
    F(m, jx+jxd,jy-jyd,jz    ) += fnq*g1x*gmy*g0z * (val);		\
    F(m, jx-jxd,jy    ,jz    ) += fnq*gmx*g0y*g0z * (val);		\
    F(m, jx    ,jy    ,jz    ) += fnq*g0x*g0y*g0z * (val);		\
    F(m, jx+jxd,jy    ,jz    ) += fnq*g1x*g0y*g0z * (val);		\
    F(m, jx-jxd,jy+jyd,jz    ) += fnq*gmx*g1y*g0z * (val);		\
    F(m, jx    ,jy+jyd,jz    ) += fnq*g0x*g1y*g0z * (val);		\
    F(m, jx+jxd,jy+jyd,jz    ) += fnq*g1x*g1y*g0z * (val);		\
    F(m, jx-jxd,jy-jyd,jz+jzd) += fnq*gmx*gmy*g1z * (val);		\
    F(m, jx    ,jy-jyd,jz+jzd) += fnq*g0x*gmy*g1z * (val);		\
    F(m, jx+jxd,jy-jyd,jz+jzd) += fnq*g1x*gmy*g1z * (val);		\
    F(m, jx-jxd,jy    ,jz+jzd) += fnq*gmx*g0y*g1z * (val);		\
    F(m, jx    ,jy    ,jz+jzd) += fnq*g0x*g0y*g1z * (val);		\
    F(m, jx+jxd,jy    ,jz+jzd) += fnq*g1x*g0y*g1z * (val);		\
    F(m, jx-jxd,jy+jyd,jz+jzd) += fnq*gmx*g1y*g1z * (val);		\
    F(m, jx    ,jy+jyd,jz+jzd) += fnq*g0x*g1y*g1z * (val);		\
    F(m, jx+jxd,jy+jyd,jz+jzd) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

// FIXME, this function exists about 100x all over the place, should
// be consolidated

template<typename particle_t>
static inline void
particle_calc_vxi(particle_t *part, typename particle_t::real_t vxi[3])
{
  typename particle_t::real_t root =
    1.f / std::sqrt(1.f + sqr(part->p[0]) + sqr(part->p[1]) + sqr(part->p[2]));
  vxi[0] = part->p[0] * root;
  vxi[1] = part->p[1] * root;
  vxi[2] = part->p[2] * root;
}

