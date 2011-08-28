
// ======================================================================
// field caching

#if DIM == DIM_Z
__shared__ real z_fld_cache[6 * (BLOCKSIZE_Z + 4) * 1 * 1];
#elif DIM == DIM_YZ
__shared__ real yz_fld_cache[6 * (BLOCKSIZE_Z + 4) * (BLOCKSIZE_Y + 4) * 1];
#else
#error unknown DIM
#endif

#if DIM == DIM_Z

#define F3C(fldnr, jx,jy,jz)						\
  (*({									\
      int off = ((((fldnr-EX)						\
		   *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))			\
		  *1 + (0))						\
		 *1 + (0));						\
      &(z_fld_cache[off]);						\
    }))

#elif DIM == DIM_YZ

#if 1
#define F3C(fldnr, jx,jy,jz)						\
  (*({									\
      int off = ((((fldnr-EX)						\
		   *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))			\
		  *(BLOCKSIZE_Y + 4) + ((jy)-(-2)))			\
		 *1 + ((jx)));						\
      &(yz_fld_cache[off]);						\
    }))
#else
#define F3C(fldnr, jx,jy,jz)						\
  (*({									\
      int off = ((((fldnr-EX)						\
		   *(BLOCKSIZE_Z+4) + ((jz)-(-2)))			\
		  *(BLOCKSIZE_Y+4) + ((jy)-(-2)))			\
		 *1 + ((jx)));						\
      if (jx != 0) printf("!!! jx\n");					\
      if (jy < -2 || jy >= BLOCKSIZE_Y + 2) printf("!!! jy %d\n", jy);	\
      if (jz < -2 || jz >= BLOCKSIZE_Z + 2) printf("!!! jz %d\n", jz);	\
      &(yz_fld_cache[off]);						\
    }))
#endif

#endif

#if DIM == DIM_Z
__device__ static void
cache_fields(real *d_flds, const int l[3])
{
#ifdef __CUDACC__
  int ti = threadIdx.x;
  int n = (BLOCKSIZE_Z + 4);
  while (ti < n) {
    int tmp = ti;
    int jx = 0;
    int jy = 0;
    int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
    //    tmp /= BLOCKSIZE_Z + 4;
    //    int m = tmp + EX;
    //    printf("n %d ti %d m %d, jx %d,%d,%d\n", n, ti, m, jx, jy, jz);
    // currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3C(m, jx,jy,jz) = F3(m, jx+l[0],jy+l[1],jz+l[2]);
    }
    ti += blockDim.x;
  }
  __syncthreads();
#else
  int ti = threadIdx.x;
  int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
  if (ti == 0) {
    for (ti = 0; ti < n; ti++) {
      int tmp = ti;
      int jx = 0;
      int jy = 0;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      //    tmp /= BLOCKSIZE_Z + 4;
      //    int m = tmp + EX;
      //    printf("n %d ti %d m %d, jx %d,%d,%d\n", n, ti, m, jx, jy, jz);
      // currently it seems faster to do the loop rather than do m by threadidx
      for (int m = EX; m <= HZ; m++) {
	F3C(m, jx,jy,jz) = F3(m, jx+l[0],jy+l[1],jz+l[2]);
      }
    }
  }
#endif
}

#elif DIM == DIM_YZ

__device__ static void
cache_fields(real *d_flds, const int l[3])
{
#ifdef __CUDACC__
  int ti = threadIdx.x;
  int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
  while (ti < n) {
    int tmp = ti;
    int jx = tmp % BLOCKSIZE_X;
    tmp /= BLOCKSIZE_X;
    int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
    tmp /= BLOCKSIZE_Y + 4;
    int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
    //    tmp /= BLOCKSIZE_Z + 4;
    //    int m = tmp + EX;
    //    printf("n %d ti %d m %d, jx %d,%d,%d\n", n, ti, m, jx, jy, jz);
    // currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3C(m, jx,jy,jz) = F3_DEV(m, jx+l[0],jy+l[1],jz+l[2]);
    }
    ti += blockDim.x;
  }
  __syncthreads();
#else
  int ti = threadIdx.x;
  int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
  if (ti == 0) {
    for (ti = 0; ti < n; ti++) {
      int tmp = ti;
      int jx = tmp % BLOCKSIZE_X;
      tmp /= BLOCKSIZE_X;
      int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
      tmp /= BLOCKSIZE_Y + 4;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      //    tmp /= BLOCKSIZE_Z + 4;
      //    int m = tmp + EX;
      //    printf("n %d ti %d m %d, jx %d,%d,%d\n", n, ti, m, jx, jy, jz);
      // currently it seems faster to do the loop rather than do m by threadidx
      int m;
      for (m = EX; m <= HZ; m++) {
	F3C(m, jx,jy,jz) = F3(m, jx+l[0],jy+l[1],jz+l[2]);
      }
    }
  }
#endif
}

#endif
