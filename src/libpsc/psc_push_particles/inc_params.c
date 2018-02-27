
#pragma once

#if (DIM & DIM_X)
#define IF_DIM_X(s) s do{} while(0)
#define IF_NOT_DIM_X(s) do{} while(0)
#else
#define IF_DIM_X(s) do{} while(0)
#define IF_NOT_DIM_X(s) s do{} while(0)
#endif

#if (DIM & DIM_Y)
#define IF_DIM_Y(s) s do{} while(0)
#define IF_NOT_DIM_Y(s) do{} while(0)
#else
#define IF_DIM_Y(s) do{} while(0)
#define IF_NOT_DIM_Y(s) s do{} while(0)
#endif

#if (DIM & DIM_Z)
#define IF_DIM_Z(s) s do{} while(0)
#define IF_NOT_DIM_Z(s) do{} while(0)
#else
#define IF_DIM_Z(s) do{} while(0)
#define IF_NOT_DIM_Z(s) s do{} while(0)
#endif

#if CALC_J == CALC_J_1VB_SPLIT
using opt_calcj = opt_calcj_1vb_split;
#elif CALC_J == CALC_J_1VB_VAR1
using opt_calcj = opt_calcj_1vb_var1;
#elif CALC_J == CALC_J_1VB_2D
using opt_calcj = opt_calcj_1vb_2d;
#endif

#define CUDA_CONSTANT
#define CUDA_DEVICE
#define __forceinline__
#define atomicAdd(addr, val) \
  do { *(addr) += (val); } while (0)


