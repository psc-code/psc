
#ifndef KG_MACROS_H
#define KG_MACROS_H

#ifdef __CUDACC__
#define KG_INLINE __device__ __host__ inline
#else
#define KG_INLINE inline
#endif

#endif
