
#pragma once

namespace psc
{

template <typename C>
class CurrentDeposition1vb
{
public:
  using Curr = C;
  using real_t = typename Curr::real_t;
  using real3_t = Vec3<real_t>;

  GT_INLINE CurrentDeposition1vb(real3_t fnqs) : fnqs_{fnqs} {}

  GT_INLINE void operator()(Curr& curr, const int i[3], real_t qni_wni,
                            const real_t dx[3], const real_t xa[3],
                            dim_xyz tag_dim) const
  {
    real_t h = (1.f / real_t(12.f)) * dx[0] * dx[1] * dx[2];
    real3_t fnq = qni_wni * fnqs_;

    curr.add(0, i[0], i[1], i[2],
             fnq[0] * (dx[0] * (1.f - xa[1]) * (1.f - xa[2]) + h));
    curr.add(0, i[0], i[1] + 1, i[2],
             fnq[0] * (dx[0] * (xa[1]) * (1.f - xa[2]) - h));
    curr.add(0, i[0], i[1], i[2] + 1,
             fnq[0] * (dx[0] * (1.f - xa[1]) * (xa[2]) - h));
    curr.add(0, i[0], i[1] + 1, i[2] + 1,
             fnq[0] * (dx[0] * (xa[1]) * (xa[2]) + h));

    curr.add(1, i[0], i[1], i[2],
             fnq[1] * (dx[1] * (1.f - xa[0]) * (1.f - xa[2]) + h));
    curr.add(1, i[0] + 1, i[1], i[2],
             fnq[1] * (dx[1] * (xa[0]) * (1.f - xa[2]) - h));
    curr.add(1, i[0], i[1], i[2] + 1,
             fnq[1] * (dx[1] * (1.f - xa[0]) * (xa[2]) - h));
    curr.add(1, i[0] + 1, i[1], i[2] + 1,
             fnq[1] * (dx[1] * (xa[0]) * (xa[2]) + h));

    curr.add(2, i[0], i[1], i[2],
             fnq[2] * (dx[2] * (1.f - xa[0]) * (1.f - xa[1]) + h));
    curr.add(2, i[0] + 1, i[1], i[2],
             fnq[2] * (dx[2] * (xa[0]) * (1.f - xa[1]) - h));
    curr.add(2, i[0], i[1] + 1, i[2],
             fnq[2] * (dx[2] * (1.f - xa[0]) * (xa[1]) - h));
    curr.add(2, i[0] + 1, i[1] + 1, i[2],
             fnq[2] * (dx[2] * (xa[0]) * (xa[1]) + h));
  };

private:
  real3_t fnqs_;
};

}; // namespace psc
