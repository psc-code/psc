#pragma once

#include "Vec3.h"

/**
 * @brief Facilitates iteration over `N`-dimensional spaces in row-major order.
 *
 * For example,
 * ```c++
 * for (Int3 i3 : VecRange(mins, maxs)) {
 *     // ...
 * }
 * ```
 * is equivalent to
 * ```c++
 * for (int ix = mins[0]; ix < maxs[0]; ix++) {
 *     for (int iy = mins[1]; iy < maxs[1]; iy++) {
 *         for (int iz = mins[2]; iz < maxs[2]; iz++) {
 *             Int3 i3{ix, iy, iz};
 *             // ...
 *         }
 *     }
 * }
 * ```
 * @tparam T scalar type (e.g. `int`)
 * @tparam N number of dimensions
 */
template <typename T, std::size_t N>
class VecRange
{
public:
  using Vec = kg::Vec<T, N>;

  class Iterator
  {
  public:
    Iterator(const Vec& starts, const Vec& stops, bool end = false)
      : starts_(starts), stops_(stops), current_(starts), end_(end)
    {}

    Iterator& operator++()
    {
      for (int d = N - 1; d >= 0; d--) {
        ++current_[d];

        if (current_[d] < stops_[d]) {
          break;
        } else {
          current_[d] = starts_[d];

          if (d == 0) {
            end_ = true;
          }
        }
      }
      return *this;
    }

    Vec operator*() { return current_; }

    bool operator!=(const Iterator& other) const { return end_ != other.end_; }

  private:
    const Vec& starts_;
    const Vec& stops_;
    Vec current_;
    bool end_;
  };

  VecRange(const Vec& starts, const Vec& stops) : starts_(starts), stops_(stops)
  {}

  Iterator begin() const { return Iterator(starts_, stops_, false); }

  Iterator end() const { return Iterator(starts_, stops_, true); }

private:
  const Vec starts_;
  const Vec stops_;
};