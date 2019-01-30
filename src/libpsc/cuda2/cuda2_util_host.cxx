
#include "psc_cuda2.h"

EXTERN_C void *
cuda_calloc(size_t nmemb, size_t size)
{
  return NULL;
}

EXTERN_C void
cuda_free(void *ptr)
{
}

EXTERN_C void
cuda_memcpy_host_from_device(void *h_ptr, void *d_ptr, size_t n)
{
  assert(0);
}

EXTERN_C void
cuda_memcpy_device_from_host(void *d_ptr, void *h_ptr, size_t n)
{
  assert(0);
}

