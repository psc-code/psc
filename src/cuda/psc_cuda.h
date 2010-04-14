
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

#include "psc.h"

typedef struct {
  float x, y, z, w;
} float4;

struct psc_cuda {
  float4 *xi4;   // xi , yi , zi , qni_div_mni
  float4 *pxi4;  // pxi, pyi, pzi, qni_wni (if qni==0, then qni_wni = wni)
};

#endif
