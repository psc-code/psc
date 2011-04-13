
#include "psc_cuda.h"

struct psc_ops psc_ops_cuda = {
  .name                   = "cuda",
  .push_part_yz_a         = cuda_push_part_yz_a,
  .push_part_yz_b         = cuda_push_part_yz_b2,
};
