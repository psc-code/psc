
#pragma once

namespace psc
{

template <typename T>
class atomic
{
  using self_type = atomic<T>;

public:
  atomic() = default;
  atomic(T value) : value_(value) {}

  operator T() const { return value_; }

  self_type& operator+=(T other)
  {
    value_ += other;
    return *this;
  }

private:
  T value_;
};

} // namespace psc
