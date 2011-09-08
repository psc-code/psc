
#include <psc_cuda.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

EXTERN_C void
sort_pairs_host(unsigned int *_d_keys, unsigned int *_d_vals, int n)
{
  thrust::device_ptr<unsigned int> d_keys(_d_keys);
  thrust::device_ptr<unsigned int> d_vals(_d_vals);

  thrust::host_vector<unsigned int> h_keys(d_keys, d_keys + n);
  thrust::host_vector<unsigned int> h_vals(d_vals, d_vals + n);

  thrust::sort_by_key(h_keys.begin(), h_keys.end(), h_vals.begin());

  thrust::copy(h_keys.begin(), h_keys.end(), d_keys);
  thrust::copy(h_vals.begin(), h_vals.end(), d_vals);
}

EXTERN_C void
sort_pairs_device(unsigned int *_d_keys, unsigned int *_d_vals, int n)
{
  thrust::device_ptr<unsigned int> d_keys(_d_keys);
  thrust::device_ptr<unsigned int> d_vals(_d_vals);
  thrust::sort_by_key(d_keys, d_keys + n, d_vals);

  thrust::host_vector<unsigned int> h1_vals(d_vals, d_vals + n);
}

