
#ifndef VPIC_ACCUMULATOR_OPS_H
#define VPIC_ACCUMULATOR_OPS_H

template<class A>
struct VpicAccumulatorOps {
  typedef A Accumulator;

  void clear_accumulator_array(Accumulator *accumulator)
  {
    TIC ::clear_accumulator_array(accumulator); TOC(clear_accumulators, 1);
  }

  void reduce_accumulator_array(Accumulator *accumulator)
  {
    TIC ::reduce_accumulator_array(accumulator); TOC(reduce_accumulators, 1);
  }
  
};


#endif
