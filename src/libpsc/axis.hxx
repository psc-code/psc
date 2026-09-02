#pragma once

namespace psc
{

struct Axis
{
  static Axis X;
  static Axis Y;
  static Axis Z;

  Axis(int axis) : axis{axis} {}

  operator int() const { return axis; }

  Axis& operator++()
  {
    this->axis += 1;
    return *this;
  }

  Axis operator++(int)
  {
    Axis temp = *this;
    ++(*this);
    return temp;
  }

  Axis next()
  {
    int next = (axis + 1) % 3;
    return next;
  }

  Axis prev()
  {
    int next = (axis + 2) % 3;
    return next;
  }

  int axis;
};

Axis Axis::X = Axis(0);
Axis Axis::Y = Axis(1);
Axis Axis::Z = Axis(2);

} // namespace psc