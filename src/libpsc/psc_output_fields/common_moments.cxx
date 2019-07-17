
#pragma once

#include "fields_item.hxx"

#include <string>

template <typename FE, typename PRT, typename R>
void deposit_1st_cc(const Grid_t& grid, FE& flds, const PRT& prt, int m, R fnqs,
                    R dxi, R dyi, R dzi, R val)
{
  auto xi = prt.x(); /* don't shift back in time */
  R u = xi[0] * dxi - .5f;
  R v = xi[1] * dyi - .5f;
  R w = xi[2] * dzi - .5f;

  int jx = fint(u);
  int jy = fint(v);
  int jz = fint(w);
  R h1 = u - jx;
  R h2 = v - jy;
  R h3 = w - jz;

  R g0x = 1.f - h1;
  R g0y = 1.f - h2;
  R g0z = 1.f - h3;
  R g1x = h1;
  R g1y = h2;
  R g1z = h3;

  int jxd = 1, jyd = 1, jzd = 1;
  if (grid.isInvar(0)) {
    jx = 0;
    g0x = 1.;
    g1x = 0.;
    jxd = 0;
  }
  if (grid.isInvar(1)) {
    jy = 0;
    g0y = 1.;
    g1y = 0.;
    jyd = 0;
  }
  if (grid.isInvar(2)) {
    jz = 0;
    g0z = 1.;
    g1z = 0.;
    jzd = 0;
  }

  assert(jx >= -1 && jx < grid.ldims[0]);
  assert(jy >= -1 && jy < grid.ldims[1]);
  assert(jz >= -1 && jz < grid.ldims[2]);

  R fnq = prt.w() * fnqs;

  auto F = makeFields3d(flds);						\
  F(m, jx, jy, jz) += fnq * g0x * g0y * g0z * (val);
  F(m, jx + jxd, jy, jz) += fnq * g1x * g0y * g0z * (val);
  F(m, jx, jy + jyd, jz) += fnq * g0x * g1y * g0z * (val);
  F(m, jx + jxd, jy + jyd, jz) += fnq * g1x * g1y * g0z * (val);
  F(m, jx, jy, jz + jzd) += fnq * g0x * g0y * g1z * (val);
  F(m, jx + jxd, jy, jz + jzd) += fnq * g1x * g0y * g1z * (val);
  F(m, jx, jy + jyd, jz + jzd) += fnq * g0x * g1y * g1z * (val);
  F(m, jx + jxd, jy + jyd, jz + jzd) += fnq * g1x * g1y * g1z * (val);
}

#define DEPOSIT_TO_GRID_1ST_CC(prt, flds, m, val)                              \
  do {                                                                         \
    deposit_1st_cc(grid, flds, prt, m, fnqs, dxi, dyi, dzi, val);	\
  } while (0)

#define DEPOSIT_TO_GRID_1ST_NC(prt, flds, m, val)                              \
  do {                                                                         \
    auto F = makeFields3d(flds);                                               \
    auto xi = prt.x(); /* don't shift back in time */                          \
    real_t u = xi[0] * dxi;                                                    \
    real_t v = xi[1] * dyi;                                                    \
    real_t w = xi[2] * dzi;                                                    \
    int jx = fint(u);                                                          \
    int jy = fint(v);                                                          \
    int jz = fint(w);                                                          \
    real_t h1 = u - jx;                                                        \
    real_t h2 = v - jy;                                                        \
    real_t h3 = w - jz;                                                        \
                                                                               \
    real_t g0x = 1.f - h1;                                                     \
    real_t g0y = 1.f - h2;                                                     \
    real_t g0z = 1.f - h3;                                                     \
    real_t g1x = h1;                                                           \
    real_t g1y = h2;                                                           \
    real_t g1z = h3;                                                           \
                                                                               \
    int jxd = 1, jyd = 1, jzd = 1;                                             \
    if (grid.isInvar(0)) {                                                     \
      jx = 0;                                                                  \
      g0x = 1.;                                                                \
      g1x = 0.;                                                                \
      jxd = 0;                                                                 \
    }                                                                          \
    if (grid.isInvar(1)) {                                                     \
      jy = 0;                                                                  \
      g0y = 1.;                                                                \
      g1y = 0.;                                                                \
      jyd = 0;                                                                 \
    }                                                                          \
    if (grid.isInvar(2)) {                                                     \
      jz = 0;                                                                  \
      g0z = 1.;                                                                \
      g1z = 0.;                                                                \
      jzd = 0;                                                                 \
    }                                                                          \
                                                                               \
    assert(jx >= -1 && jx < grid.ldims[0]);                                    \
    assert(jy >= -1 && jy < grid.ldims[1]);                                    \
    assert(jz >= -1 && jz < grid.ldims[2]);                                    \
                                                                               \
    real_t fnq = prt.w() * fnqs;                                               \
                                                                               \
    F(m, jx, jy, jz) += fnq * g0x * g0y * g0z * (val);                         \
    F(m, jx + jxd, jy, jz) += fnq * g1x * g0y * g0z * (val);                   \
    F(m, jx, jy + jyd, jz) += fnq * g0x * g1y * g0z * (val);                   \
    F(m, jx + jxd, jy + jyd, jz) += fnq * g1x * g1y * g0z * (val);             \
    F(m, jx, jy, jz + jzd) += fnq * g0x * g0y * g1z * (val);                   \
    F(m, jx + jxd, jy, jz + jzd) += fnq * g1x * g0y * g1z * (val);             \
    F(m, jx, jy + jyd, jz + jzd) += fnq * g0x * g1y * g1z * (val);             \
    F(m, jx + jxd, jy + jyd, jz + jzd) += fnq * g1x * g1y * g1z * (val);       \
  } while (0)

#define DEPOSIT_TO_GRID_2ND_NC(prt, flds, m, val)                              \
  do {                                                                         \
    auto F = makeFields3d(flds);                                               \
    auto xi = prt.x(); /* don't shift back in time */                          \
    real_t u = xi[0] * dxi;                                                    \
    real_t v = xi[1] * dyi;                                                    \
    real_t w = xi[2] * dzi;                                                    \
    int jx = nint(u);                                                          \
    int jy = nint(v);                                                          \
    int jz = nint(w);                                                          \
    real_t h1 = jx - u;                                                        \
    real_t h2 = jy - v;                                                        \
    real_t h3 = jz - w;                                                        \
                                                                               \
    real_t gmx = .5f * (.5f + h1) * (.5f + h1);                                \
    real_t gmy = .5f * (.5f + h2) * (.5f + h2);                                \
    real_t gmz = .5f * (.5f + h3) * (.5f + h3);                                \
    real_t g0x = .75f - h1 * h1;                                               \
    real_t g0y = .75f - h2 * h2;                                               \
    real_t g0z = .75f - h3 * h3;                                               \
    real_t g1x = .5f * (.5f - h1) * (.5f - h1);                                \
    real_t g1y = .5f * (.5f - h2) * (.5f - h2);                                \
    real_t g1z = .5f * (.5f - h3) * (.5f - h3);                                \
                                                                               \
    int jxd = 1, jyd = 1, jzd = 1;                                             \
    if (grid.isInvar(0)) {                                                     \
      jx = 0;                                                                  \
      g0x = 1.;                                                                \
      g1x = 0.;                                                                \
      gmx = 0.;                                                                \
      jxd = 0;                                                                 \
    }                                                                          \
    if (grid.isInvar(1)) {                                                     \
      jy = 0;                                                                  \
      g0y = 1.;                                                                \
      g1y = 0.;                                                                \
      gmy = 0.;                                                                \
      jyd = 0;                                                                 \
    }                                                                          \
    if (grid.isInvar(2)) {                                                     \
      jz = 0;                                                                  \
      g0z = 1.;                                                                \
      g1z = 0.;                                                                \
      gmz = 0.;                                                                \
      jzd = 0;                                                                 \
    }                                                                          \
                                                                               \
    assert(jx >= 0 && jx <= grid.ldims[0]);                                    \
    assert(jy >= 0 && jy <= grid.ldims[1]);                                    \
    assert(jz >= 0 && jz <= grid.ldims[2]);                                    \
                                                                               \
    real_t fnq = prt.w() * fnqs;                                               \
                                                                               \
    F(m, jx - jxd, jy - jyd, jz - jzd) += fnq * gmx * gmy * gmz * (val);       \
    F(m, jx, jy - jyd, jz - jzd) += fnq * g0x * gmy * gmz * (val);             \
    F(m, jx + jxd, jy - jyd, jz - jzd) += fnq * g1x * gmy * gmz * (val);       \
    F(m, jx - jxd, jy, jz - jzd) += fnq * gmx * g0y * gmz * (val);             \
    F(m, jx, jy, jz - jzd) += fnq * g0x * g0y * gmz * (val);                   \
    F(m, jx + jxd, jy, jz - jzd) += fnq * g1x * g0y * gmz * (val);             \
    F(m, jx - jxd, jy + jyd, jz - jzd) += fnq * gmx * g1y * gmz * (val);       \
    F(m, jx, jy + jyd, jz - jzd) += fnq * g0x * g1y * gmz * (val);             \
    F(m, jx + jxd, jy + jyd, jz - jzd) += fnq * g1x * g1y * gmz * (val);       \
    F(m, jx - jxd, jy - jyd, jz) += fnq * gmx * gmy * g0z * (val);             \
    F(m, jx, jy - jyd, jz) += fnq * g0x * gmy * g0z * (val);                   \
    F(m, jx + jxd, jy - jyd, jz) += fnq * g1x * gmy * g0z * (val);             \
    F(m, jx - jxd, jy, jz) += fnq * gmx * g0y * g0z * (val);                   \
    F(m, jx, jy, jz) += fnq * g0x * g0y * g0z * (val);                         \
    F(m, jx + jxd, jy, jz) += fnq * g1x * g0y * g0z * (val);                   \
    F(m, jx - jxd, jy + jyd, jz) += fnq * gmx * g1y * g0z * (val);             \
    F(m, jx, jy + jyd, jz) += fnq * g0x * g1y * g0z * (val);                   \
    F(m, jx + jxd, jy + jyd, jz) += fnq * g1x * g1y * g0z * (val);             \
    F(m, jx - jxd, jy - jyd, jz + jzd) += fnq * gmx * gmy * g1z * (val);       \
    F(m, jx, jy - jyd, jz + jzd) += fnq * g0x * gmy * g1z * (val);             \
    F(m, jx + jxd, jy - jyd, jz + jzd) += fnq * g1x * gmy * g1z * (val);       \
    F(m, jx - jxd, jy, jz + jzd) += fnq * gmx * g0y * g1z * (val);             \
    F(m, jx, jy, jz + jzd) += fnq * g0x * g0y * g1z * (val);                   \
    F(m, jx + jxd, jy, jz + jzd) += fnq * g1x * g0y * g1z * (val);             \
    F(m, jx - jxd, jy + jyd, jz + jzd) += fnq * gmx * g1y * g1z * (val);       \
    F(m, jx, jy + jyd, jz + jzd) += fnq * g0x * g1y * g1z * (val);             \
    F(m, jx + jxd, jy + jyd, jz + jzd) += fnq * g1x * g1y * g1z * (val);       \
  } while (0)

// FIXME, this function exists about 100x all over the place, should
// be consolidated

template <typename Particle>
static inline void particle_calc_vxi(const Particle& prt,
                                     typename Particle::real_t vxi[3])
{
  typename Particle::real_t root =
    1.f / std::sqrt(1.f + sqr(prt.u()[0]) + sqr(prt.u()[1]) + sqr(prt.u()[2]));
  vxi[0] = prt.u()[0] * root;
  vxi[1] = prt.u()[1] * root;
  vxi[2] = prt.u()[2] * root;
}
