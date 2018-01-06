
#ifndef FIELDS3D_HXX
#define FIELDS3D_HXX

#include <type_traits>

template<bool AOS>
struct Layout
{
  using isAOS = std::integral_constant<bool, AOS>;
};

using LayoutAOS = Layout<true>;
using LayoutSOA = Layout<false>;

template<typename R, typename L=LayoutSOA>
struct fields3d {
  using real_t = R;
  using layout = L;

  real_t *data;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  int first_comp; // first component

  real_t  operator()(int m, int i, int j, int k) const { return data[index(m, i, j, k)];  }
  real_t& operator()(int m, int i, int j, int k)       { return data[index(m, i, j, k)];  }

  int index(int m, int i, int j, int k);

  int size()
  {
    return nr_comp * im[0] * im[1] * im[2];
  }

  int n_cells()
  {
    return im[0] * im[1] * im[2];
  }

  void zero(int m)
  {
    memset(&(*this)(m, ib[0], ib[1], ib[2]), 0, n_cells() * sizeof(real_t));
  }

  void zero(int mb, int me)
  {
    for (int m = mb; m < me; m++) {
      zero(m);
    }
  }

  void zero()
  {
    memset(data, 0, sizeof(real_t) * size());
  }
};

template<typename R, typename L>
int fields3d<R, L>::index(int m, int i, int j, int k)
{
#ifdef BOUNDS_CHECK
  assert(m >= first_comp_ && m < n_comp_);
  assert(i >= ib[0] && i < ib[0] + im[0]);
  assert(j >= ib[1] && j < ib[1] + im[1]);
  assert(k >= ib[2] && k < ib[2] + im[2]);
#endif

  if (L::isAOS::value) {
    return (((((k - ib[2])) * im[1] +
	      (j - ib[1])) * im[0] +
	     (i - ib[0])) * nr_comp + m);
  } else {
    return (((((m - first_comp) * im[2] +
	       (k - ib[2])) * im[1] +
	      (j - ib[1])) * im[0] +
	     (i - ib[0])));
  }
}

#endif

