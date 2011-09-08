
#include <psc_cuda.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>

EXTERN_C void
sort_pairs_device(unsigned int *_d_keys, unsigned int *_d_vals, int n)
{
  thrust::device_ptr<unsigned int> d_keys(_d_keys);
  thrust::device_ptr<unsigned int> d_vals(_d_vals);
  thrust::sort_by_key(d_keys, d_keys + n, d_vals);
}

EXTERN_C void
sort_pairs_host(unsigned int *keys, unsigned int *vals, int n)
{
#if 0
  thrust::sort_by_key(keys, keys + n, vals);
#else
  thrust::device_vector<unsigned int> d_keys(keys, keys + n);
  thrust::device_vector<unsigned int> d_vals(vals, vals + n);
  sort_pairs_device(thrust::raw_pointer_cast(&d_keys[0]),
		    thrust::raw_pointer_cast(&d_vals[0]), n);
  //  thrust::copy(d_keys.begin(), d_keys.end(), keys);
  thrust::copy(d_vals.begin(), d_vals.end(), vals);
#endif
}
