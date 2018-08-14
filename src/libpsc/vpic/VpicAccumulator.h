
#pragma once

template<typename MfieldsAccumulator, typename FieldArray>
struct VpicAccumulatorOps
{
  static void clear(MfieldsAccumulator& acc)                        { ::clear_accumulator_array(acc); }
  static void reduce(MfieldsAccumulator& acc)                       { ::reduce_accumulator_array(acc); }
  static void unload(const MfieldsAccumulator& acc, FieldArray& fa) { ::unload_accumulator_array(&fa, acc); }
};

