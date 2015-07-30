
#include "psc_cuda2.h"

#define check(a) do { int ierr = a; if (ierr != cudaSuccess) fprintf(stderr, "IERR = %d (%d)\n", ierr, cudaSuccess); assert(ierr == cudaSuccess); } while(0)

EXTERN_C void *
cuda_calloc(size_t nmemb, size_t size)
{
  // FIXME, doesn't zero the memory
  void *ptr;
  check(cudaMalloc(&ptr, nmemb * size));
  return ptr;
}

EXTERN_C void
cuda_free(void *ptr)
{
  check(cudaFree(ptr));
}

EXTERN_C void
cuda_memcpy_host_from_device(void *h_ptr, void *d_ptr, size_t n)
{
  check(cudaMemcpy(h_ptr, d_ptr, n, cudaMemcpyDeviceToHost));
}

EXTERN_C void
cuda_memcpy_device_from_host(void *d_ptr, void *h_ptr, size_t n)
{
  check(cudaMemcpy(d_ptr, h_ptr, n, cudaMemcpyHostToDevice));
}

