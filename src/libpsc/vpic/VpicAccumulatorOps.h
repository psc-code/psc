
#ifndef VPIC_ACCUMULATOR_OPS_H
#define VPIC_ACCUMULATOR_OPS_H

template<class A, class FA>
struct VpicAccumulatorOps {
  typedef A Accumulator;
  typedef FA FieldArray;

  void clear_accumulator_array(Accumulator *accumulator)
  {
    TIC ::clear_accumulator_array(accumulator); TOC(clear_accumulators, 1);
  }

  void reduce_accumulator_array(Accumulator *accumulator)
  {
    TIC ::reduce_accumulator_array(accumulator); TOC(reduce_accumulators, 1);
  }
  
  void unload_accumulator_array(FieldArray *fa, Accumulator *accumulator)
  {
    TIC ::unload_accumulator_array(fa, accumulator); TOC(unload_accumulator, 1);
  }
  
};


#endif
