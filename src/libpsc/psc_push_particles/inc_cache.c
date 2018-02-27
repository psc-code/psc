
// ======================================================================
// field caching

#ifndef EM_CACHE
#define EM_CACHE EM_CACHE_NONE
#endif

// OPT, shouldn't we be able to do with less ghosts?
#if DIM == DIM_YZ
#define BLOCKBND_X 0
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#elif DIM == DIM_XYZ
#define BLOCKBND_X 2
#define BLOCKBND_Y 2
#define BLOCKBND_Z 2
#endif

#define BLOCKGSIZE_X (BLOCKSIZE_X + 2 * BLOCKBND_X)
#define BLOCKGSIZE_Y (BLOCKSIZE_Y + 2 * BLOCKBND_Y)
#define BLOCKGSIZE_Z (BLOCKSIZE_Z + 2 * BLOCKBND_Z)

