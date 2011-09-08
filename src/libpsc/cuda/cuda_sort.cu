
#include <psc_cuda.h>
#include <thrust/sort.h>

EXTERN_C void
sort_pairs_host(unsigned int *keys, unsigned int *vals, int n)
{
  thrust::sort_by_key(keys, keys + n, vals);
}
