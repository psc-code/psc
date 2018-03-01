
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

#define CUDA_CONSTANT
#define CUDA_DEVICE
#define __forceinline__


