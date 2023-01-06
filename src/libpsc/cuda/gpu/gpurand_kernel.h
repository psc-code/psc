
#ifdef GTENSOR_DEVICE_HIP
#include <hiprand_kernel.h>
#else
#include <curand_kernel.h>
#define hiprand_init curand_init
#define hiprandState curandState
#define hiprand_uniform curand_uniform
#define hiprand_normal curand_normal
#define hiprand_normal2 curand_normal2
#endif
