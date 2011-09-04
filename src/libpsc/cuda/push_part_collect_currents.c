
// ----------------------------------------------------------------------
// collect_currents

__global__ static void
collect_currents(real *d_flds, real *d_scratch, int nr_blocks)
{
#if DIM == DIM_Z
  int jy = 0;
  int jz = threadIdx.x - SW;
#elif DIM == DIM_YZ
  int jy = threadIdx.x - SW;
  int jz = threadIdx.y - SW;
#endif

  for (int b = 0; b < nr_blocks; b++) {
    real *scratch = d_scratch + b * 3 * BLOCKSTRIDE;

    int ci[3];
    blockIdx_to_blockCrd(b, ci);
    ci[0] *= BLOCKSIZE_X;
    ci[1] *= BLOCKSIZE_Y;
    ci[2] *= BLOCKSIZE_Z;
    ci[0] += d_ilo[0];
    ci[1] += d_ilo[1];
    ci[2] += d_ilo[2];

    if (threadIdx.x == 0) {
      for (int m = 0; m < 3; m++) {
	for (jz = -SW; jz < BLOCKSIZE_Z + SW; jz++) {
	  for (jy = -SW; jy < BLOCKSIZE_Y + SW; jy++) {
	    F3_DEV(JXI+m, 0+ci[0],jy+ci[1],jz+ci[2]) += scratch(m,jy,jz);
	  }
	}
      }
    }
  }
}

