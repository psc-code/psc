
#pragma once

namespace psc
{

template <typename C>
class CurrentDeposition1vb
{
public:
  using Curr = C;
  using real_t = typename Curr::real_t;
  using Real3 = Vec3<real_t>;

  GT_INLINE CurrentDeposition1vb(Real3 fnqs) : fnqs_{fnqs} {}

  GT_INLINE void operator()(Curr& curr, const Int3& i, real_t qni_wni,
                            const Real3& dx, const Real3& xa,
                            dim_xyz tag_dim) const
  {
    real_t h = (1.f / real_t(12.f)) * dx.prod();
    Real3 fnq = qni_wni * fnqs_;

    // clang-format off
    curr.add(0, i[0],     i[1],     i[2],     fnq[0] * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
    curr.add(0, i[0],     i[1] + 1, i[2],     fnq[0] * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
    curr.add(0, i[0],     i[1],     i[2] + 1, fnq[0] * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
    curr.add(0, i[0],     i[1] + 1, i[2] + 1, fnq[0] * (dx[0] * (      xa[1]) * (      xa[2]) + h));

    curr.add(1, i[0],     i[1],     i[2],     fnq[1] * (dx[1] * (1.f - xa[2]) * (1.f - xa[0]) + h));
    curr.add(1, i[0],     i[1],     i[2] + 1, fnq[1] * (dx[1] * (      xa[2]) * (1.f - xa[0]) - h));
    curr.add(1, i[0] + 1, i[1],     i[2],     fnq[1] * (dx[1] * (1.f - xa[2]) * (      xa[0]) - h));
    curr.add(1, i[0] + 1, i[1],     i[2] + 1, fnq[1] * (dx[1] * (      xa[2]) * (      xa[0]) + h));

    curr.add(2, i[0],     i[1],     i[2],     fnq[2] * (dx[2] * (1.f - xa[0]) * (1.f - xa[1]) + h));
    curr.add(2, i[0] + 1, i[1],     i[2],     fnq[2] * (dx[2] * (      xa[0]) * (1.f - xa[1]) - h));
    curr.add(2, i[0],     i[1] + 1, i[2],     fnq[2] * (dx[2] * (1.f - xa[0]) * (      xa[1]) - h));
    curr.add(2, i[0] + 1, i[1] + 1, i[2],     fnq[2] * (dx[2] * (      xa[0]) * (      xa[1]) + h));
    // clang-format on
  };

  GT_INLINE void operator()(Curr& curr, const Int3& i, real_t qni_wni,
                            const Real3& dx, const Real3& xa,
                            dim_xz tag_dim) const
  {
    real_t h = (1.f / real_t(12.f)) * dx.prod();
    Real3 fnq = qni_wni * fnqs_;

    // clang-format off
    curr.add(0, i[0],     i[1],     i[2],     fnq[0] * (dx[0] *                 (1.f - xa[2])    ));
    curr.add(0, i[0],     i[1],     i[2] + 1, fnq[0] * (dx[0] *                 (      xa[2])    ));

    curr.add(1, i[0],     i[1],     i[2],     fnq[1] * (dx[1] * (1.f - xa[2]) * (1.f - xa[0]) + h));
    curr.add(1, i[0],     i[1],     i[2] + 1, fnq[1] * (dx[1] * (      xa[2]) * (1.f - xa[0]) - h));
    curr.add(1, i[0] + 1, i[1],     i[2],     fnq[1] * (dx[1] * (1.f - xa[2]) * (      xa[0]) - h));
    curr.add(1, i[0] + 1, i[1],     i[2] + 1, fnq[1] * (dx[1] * (      xa[2]) * (      xa[0]) + h));

    curr.add(2, i[0],     i[1],     i[2],     fnq[2] * (dx[2] * (1.f - xa[0])                    ));
    curr.add(2, i[0] + 1, i[1],     i[2],     fnq[2] * (dx[2] * (      xa[0])                    ));
    // clang-format on
  }

  GT_INLINE void operator()(Curr& curr, const Int3& i, real_t qni_wni,
                            const Real3& dx, const Real3& xa,
                            dim_yz tag_dim) const
  {
    real_t h = (1.f / real_t(12.f)) * dx.prod();
    Real3 fnq = qni_wni * fnqs_;

    // clang-format off
    curr.add(0, i[0],     i[1],     i[2],     fnq[0] * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
    curr.add(0, i[0],     i[1] + 1, i[2],     fnq[0] * (dx[0] * (      xa[1]) * (1.f - xa[2]) - h));
    curr.add(0, i[0],     i[1],     i[2] + 1, fnq[0] * (dx[0] * (1.f - xa[1]) * (      xa[2]) - h));
    curr.add(0, i[0],     i[1] + 1, i[2] + 1, fnq[0] * (dx[0] * (      xa[1]) * (      xa[2]) + h));

    curr.add(1, i[0],     i[1],     i[2],     fnq[1] * (dx[1] * (1.f - xa[2])                    ));
    curr.add(1, i[0],     i[1],     i[2] + 1, fnq[1] * (dx[1] * (      xa[2])                    ));

    curr.add(2, i[0],     i[1],     i[2],     fnq[2] * (dx[2] *                 (1.f - xa[1])    ));
    curr.add(2, i[0],     i[1] + 1, i[2],     fnq[2] * (dx[2] *                 (      xa[1])    ));
    // clang-format on
  }

private:
  Real3 fnqs_;
};

}; // namespace psc
