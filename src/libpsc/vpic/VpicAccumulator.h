
#ifndef VPIC_ACCUMULATOR_H
#define VPIC_ACCUMULATOR_H

template<class AccumulatorBase, class FA>
struct VpicAccumulator : AccumulatorBase
{
  typedef AccumulatorBase Base;
  typedef FA FieldArray;
  using typename Base::Grid;

  void clear()                      { ::clear_accumulator_array(this); }
  void reduce()                     { ::reduce_accumulator_array(this); }
  void unload(FieldArray& fa) const { ::unload_accumulator_array(&fa, this); }
};

#endif
