
#include "gtest/gtest.h"

#include "binary_collision.hxx"
#include "vec3.hxx"

struct TestParticle
{
  using real_t = double;

  Vec3<real_t> u;
  real_t q;
  real_t m;
};

struct TestParticleRef
{
  using real_t = double;

  TestParticleRef(TestParticle& prt)
    : prt_{prt}
  {}

  real_t q() const { return prt_.q; }
  real_t m() const { return prt_.m; }
  real_t  u(int d) const { return prt_.u[d]; }
  real_t& u(int d)       { return prt_.u[d]; }

private:
  TestParticle& prt_;
};

TEST(BinaryCollision, Test1)
{
  double eps = 1e-14;
  
  TestParticle prt1{{ 1., 0., 0.}, 1., 1. };
  TestParticle prt2{{ 0., 0., 0.}, 1., 1. };
  
  BinaryCollision<TestParticleRef> bc;

  double nudt = bc(prt1, prt2, .1);

  EXPECT_NEAR(prt1.u[0] + prt2.u[0], 1., eps);
  EXPECT_NEAR(prt1.u[1] + prt2.u[1], 0., eps);
  EXPECT_NEAR(prt1.u[2] + prt2.u[2], 0., eps);

#if 0
  printf("prt1: %g %g %g\n", prt1.u[0], prt1.u[1], prt1.u[2]);
  printf("prt2: %g %g %g\n", prt2.u[0], prt2.u[1], prt2.u[2]);
#endif
}

