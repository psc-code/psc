
#pragma once

// ----------------------------------------------------------------------
// SCurr

// OPT: take i < cell_end condition out of load
// OPT: reduce two at a time
// OPT: try splitting current calc / measuring by itself

// OPT: don't need as many ghost points for current and EM fields (?)

#define BND_CURR_L (1)
#define BND_CURR_R (2)

#define NR_CBLOCKS 16
#define CBLOCK_ID (threadIdx.x & (NR_CBLOCKS - 1))
#define CBLOCK_SIZE_Y (BS_Y + BND_CURR_L + BND_CURR_R)
#define CBLOCK_SIZE_Z (BS_Z + BND_CURR_L + BND_CURR_R)
#define CBLOCK_SIZE (CBLOCK_SIZE_Y * CBLOCK_SIZE_Z * (NR_CBLOCKS))

#define CBLOCK_OFF(jy, jz, m, wid) ((((m) * CBLOCK_SIZE_Z + ((jz) + BND_CURR_L)) * CBLOCK_SIZE_Y + ((jy) + BND_CURR_L)) * (NR_CBLOCKS) + wid)

template<typename BS>
class SCurr
{
  static const int BS_X = BS::x::value, BS_Y = BS::y::value, BS_Z = BS::z::value;

public:
  static const int shared_size = 3 * CBLOCK_SIZE;

  float *scurr;
  DFields d_flds;

  __device__ SCurr(float *_scurr, DFields _d_flds) :
    scurr(_scurr), d_flds(_d_flds)
  {
    int i = threadIdx.x;
    while (i < shared_size) {
      scurr[i] = float(0.);
      i += THREADS_PER_BLOCK;
    }
  }

  __device__ void add_to_fld(const int *ci0)
  {
    __syncthreads();				\
    int i = threadIdx.x;
    int stride = (BS_Y + BND_CURR_L + BND_CURR_R) * (BS_Z + BND_CURR_L + BND_CURR_R);
    while (i < stride) {
      int rem = i;
      int jz = rem / (BS_Y + BND_CURR_L + BND_CURR_R);
      rem -= jz * (BS_Y + BND_CURR_L + BND_CURR_R);
      int jy = rem;
      jz -= BND_CURR_L;
      jy -= BND_CURR_L;
      for (int m = 0; m < 3; m++) {
	float val = float(0.);
	// FIXME, OPT
	for (int wid = 0; wid < NR_CBLOCKS; wid++) {
	  val += (*this)(wid, jy, jz, m);
	}
	d_flds(JXI+m, 0,jy+ci0[1],jz+ci0[2]) += val;
	i += THREADS_PER_BLOCK;
      }
    }
  }

  __device__ float operator()(int wid, int jy, int jz, int m) const
  {
    uint off = CBLOCK_OFF(jy, jz, m, wid);
    return scurr[off];
  }
  __device__ float& operator()(int wid, int jy, int jz, int m)
  {
    uint off = CBLOCK_OFF(jy, jz, m, wid);
    return scurr[off];
  }

  __device__ void add(int m, int jy, int jz, float val, const int *ci0)
  {
    float *addr = &(*this)(CBLOCK_ID, jy, jz, m);
    atomicAdd(addr, val);
  }
};

// ----------------------------------------------------------------------
// GCurr

template<typename BS>
class GCurr
{
public:
  static const int shared_size = 1;

  float *scurr;
  DFields d_flds;

  __device__ GCurr(float *_scurr, DFields _d_flds) :
    scurr(_scurr), d_flds(_d_flds)
  {
  }

  __device__ void add_to_fld(int *ci0)
  {
  }

  __device__ void add(int m, int jy, int jz, float val, int *ci0)
  {
    float *addr = &d_flds(JXI+m, 0,jy+ci0[1],jz+ci0[2]);
    atomicAdd(addr, val);
  }
};

struct BlockBase
{
  int bid;
  int p;
  int ci0[3];
};

struct Block : BlockBase
{
  template<typename BS>
  __device__
  int find_block_pos_patch(const DMparticlesCuda<BS>& dmprts, int *block_pos, int block_start)
  {
    p = dmprts.find_block_pos_patch(block_pos, ci0);
    return p;
  }
};

struct BlockQ : BlockBase
{
  template<typename BS>
  __device__
  int find_block_pos_patch(const DMparticlesCuda<BS>& dmprts, int *block_pos, int block_start)
  {
    p = dmprts.find_block_pos_patch_q(block_pos, ci0, block_start);
    return p;
  }
};

// ======================================================================
// CurrmemShared

struct CurrmemShared
{
  template<typename BS>
  using Curr = SCurr<BS>;

  using Block = BlockQ;
  
  static Range<int> block_starts() { return range(4);  }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx =  (cmprts.b_mx()[1] + 1) / 2;
    int gy = ((cmprts.b_mx()[2] + 1) / 2) * cmprts.n_patches;
    return dim3(gx, gy);
  }

  template<typename BS>
  __device__ static int find_bid(DMparticlesCuda<BS>& dmprts, int p, int *block_pos)
  {
    return dmprts.find_bid_q(p, block_pos);
  }
};

// ======================================================================
// CurrmemGlobal

struct CurrmemGlobal
{
  template<typename BS>
  using Curr = GCurr<BS>;

  using Block = Block;

  static Range<int> block_starts() { return range(1);  }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx = cmprts.b_mx()[1];
    int gy = cmprts.b_mx()[2] * cmprts.n_patches;
    return dim3(gx, gy);
  }

  template<typename BS>
  __device__ static int find_bid(DMparticlesCuda<BS>& dmprts, int p, int *block_pos)
  {
    return dmprts.find_bid();
  }
};

