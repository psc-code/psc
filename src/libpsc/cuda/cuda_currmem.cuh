
#pragma once

// #define DEBUG_CURRMEM

// ----------------------------------------------------------------------
// SCurr

// OPT: take i < cell_end condition out of load
// OPT: reduce two at a time
// OPT: try splitting current calc / measuring by itself

// OPT: don't need as many ghost points for current and EM fields (?)

template <typename BS, typename DIM>
class SCurr
{
  static const int BS_X = BS::x::value, BS_Y = BS::y::value,
                   BS_Z = BS::z::value;
  static const int N_GHOSTS_L = 1;
  static const int N_GHOSTS_R = 2;
  static const uint STRIDE_Y = BS_Y + N_GHOSTS_L + N_GHOSTS_R;
  static const uint STRIDE_Z = BS_Z + N_GHOSTS_L + N_GHOSTS_R;

  static const uint N_COPIES = 16;

public:
  static const int shared_size = 3 * STRIDE_Y * STRIDE_Z * N_COPIES;
  using real_t = float;

  float* scurr;

private:
  gt::gtensor_span_device<float, 4> gt_;
  Int3 ib_;

public:
  template <typename E>
  __device__ SCurr(float* _scurr, const E& gt, const Int3& ib)
    : scurr(_scurr), gt_(gt), ib_(ib)
  {
    int i = threadIdx.x;
    while (i < shared_size) {
      scurr[i] = float(0.);
      i += THREADS_PER_BLOCK;
    }
  }

  __device__ void add_to_fld(const int* ci0)
  {
    __syncthreads();
    int i = threadIdx.x;
    int stride = STRIDE_Y * STRIDE_Z;
    auto _d_flds = make_Fields3d<dim_xyz>(gt_, ib_);
    while (i < stride) {
      int rem = i;
      int jz = rem / STRIDE_Y;
      rem -= jz * STRIDE_Y;
      int jy = rem;
      jz -= N_GHOSTS_L;
      jy -= N_GHOSTS_L;
      for (int m = 0; m < 3; m++) {
        float val = float(0.);
        // FIXME, OPT
        for (int wid = 0; wid < N_COPIES; wid++) {
          val += (*this)(wid, jy, jz, m);
        }
        _d_flds(JXI + m, 0, jy + ci0[1], jz + ci0[2]) += val;
        i += THREADS_PER_BLOCK;
      }
    }
  }

  __device__ float operator()(int wid, int jy, int jz, int m) const
  {
    uint off = index(jy, jz, m, wid);
    return scurr[off];
  }

  __device__ float& operator()(int wid, int jy, int jz, int m)
  {
    uint off = index(jy, jz, m, wid);
    return scurr[off];
  }

  __device__ void add(int m, int jy, int jz, float val, const int* ci0)
  {
    uint wid = threadIdx.x & (N_COPIES - 1);
    float* addr = &(*this)(wid, jy, jz, m);
    atomicAdd(addr, val);
  }

private:
  __device__ uint index(int jy, int jz, uint m, uint wid) const
  {
#ifdef DEBUG_CURRMEM
    if (jy < -N_GHOSTS_L || jy >= BS_Y + N_GHOSTS_R || jz < -N_GHOSTS_L ||
        jz >= BS_Z + N_GHOSTS_R || m >= 3 || wid >= N_COPIES) {
      printf("CUDA_ERROR jyz %d:%d BS %d:%d m %d/3 wid %d/%d\n", jy, jz, BS_Y,
             BS_Z, m, wid, N_COPIES);
    }
#endif
    uint off =
      ((((m)*STRIDE_Z + ((jz) + N_GHOSTS_L)) * STRIDE_Y + ((jy) + N_GHOSTS_L)) *
         (N_COPIES) +
       wid);
#ifdef DEBUG_CURRMEM
    if (off >= shared_size) {
      printf("CUDA_ERROR off %d %d wid %d %d:%d m %d\n", off, shared_size, wid,
             jy, jz, m);
    }
#endif
    return off;
  }
};

// ----------------------------------------------------------------------
// GCurr

template <typename BS>
class GCurr
{
public:
  using real_t = float;
  static const int shared_size = 1;

  float* scurr;

private:
  gt::gtensor_span_device<float, 4> gt_;
  Int3 ib_;

public:
  template <typename E>
  __device__ GCurr(float* _scurr, const E& gt, const Int3& ib)
    : scurr(_scurr), gt_(gt), ib_(ib)
  {}

  __device__ void add_to_fld(int* ci0) {}

  __device__ void add(int m, int jy, int jz, float val, const int* ci0)
  {
    auto _d_flds = make_Fields3d<dim_xyz>(gt_, ib_);
    float* addr = &_d_flds(JXI + m, 0, jy + ci0[1], jz + ci0[2]);
    atomicAdd(addr, val);
  }

  __device__ void add(int m, int jx, int jy, int jz, float val)
  {
    auto _d_flds = make_Fields3d<dim_xyz>(gt_, ib_);
    float* addr = &_d_flds(JXI + m, jx, jy, jz);
    atomicAdd(addr, val);
  }
};

// ======================================================================
// CurrmemShared

struct CurrmemShared
{
  template <typename BS, typename DIM>
  using Curr = SCurr<BS, DIM>;

  template <typename BS, typename DIM>
  using Block = BlockQ<BS, DIM>;
};

// ======================================================================
// CurrmemGlobal

struct CurrmemGlobal
{
  template <typename BS, typename DIM>
  using Curr = GCurr<BS>;

  template <typename BS, typename DIM>
  using Block = BlockSimple<BS, DIM>;
};
