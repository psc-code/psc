
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include <cstdio>
#include <cassert>

// ----------------------------------------------------------------------
// cuda_mfields_create

struct cuda_mfields *
cuda_mfields_create()
{
  struct cuda_mfields *cmflds = 
    (struct cuda_mfields *) calloc(1, sizeof(*cmflds));

  return cmflds;
}

// ----------------------------------------------------------------------
// cuda_mfields_destroy

void
cuda_mfields_destroy(struct cuda_mfields *cmflds)
{
  free(cmflds);
}

// ----------------------------------------------------------------------
// cuda_mfields_ctor

void
cuda_mfields_ctor(struct cuda_mfields *cmflds, mrc_json_t json)
{
  cudaError_t ierr;

  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  
  cmflds->n_patches = mrc_json_get_object_entry_integer(json_info, "n_patches");
  cmflds->n_fields = mrc_json_get_object_entry_integer(json_info, "n_fields");
  mrc_json_get_object_entry_int3(json_info, "ib", cmflds->ib);
  mrc_json_get_object_entry_int3(json_info, "im", cmflds->im);
  mrc_json_get_object_entry_int3(json_info, "ldims", cmflds->ldims);
  double dx[3];
  mrc_json_get_object_entry_double3(json_info, "dx", dx);
  for (int d = 0; d < 3; d++) {
    cmflds->dx[d] = dx[d];
  }

  cmflds->n_cells_per_patch = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
  cmflds->n_cells = cmflds->n_patches * cmflds->n_cells_per_patch;

  ierr = cudaMalloc((void **) &cmflds->d_flds,
		    cmflds->n_fields * cmflds->n_cells * sizeof(*cmflds->d_flds)); cudaCheck(ierr);

  cmflds->d_flds_by_patch = new fields_cuda_real_t *[cmflds->n_patches];
  for (int p = 0; p < cmflds->n_patches; p++) {
    cmflds->d_flds_by_patch[p] = cmflds->d_flds + p * cmflds->n_fields * cmflds->n_cells_per_patch;
  }
}

// ----------------------------------------------------------------------
// cuda_mfields_dtor

void
cuda_mfields_dtor(struct cuda_mfields *cmflds)
{
  cudaError_t ierr;

  ierr = cudaFree(cmflds->d_flds); cudaCheck(ierr);
  
  delete[] cmflds->d_flds_by_patch;
}

// ----------------------------------------------------------------------
// cuda_mfields_to_json

mrc_json_t
cuda_mfields_to_json(struct cuda_mfields *cmflds)
{
  mrc_json_t json = mrc_json_object_new(9);
  mrc_json_object_push_integer(json, "n_patches", cmflds->n_patches);
  mrc_json_object_push_integer(json, "n_fields", cmflds->n_fields);
  mrc_json_object_push_integer(json, "n_cells_per_patch", cmflds->n_cells_per_patch);
  mrc_json_object_push_integer(json, "n_cells", cmflds->n_cells);

  mrc_json_object_push(json, "ib", mrc_json_integer_array_new(3, cmflds->ib));
  mrc_json_object_push(json, "im", mrc_json_integer_array_new(3, cmflds->im));
  mrc_json_object_push(json, "ldims", mrc_json_integer_array_new(3, cmflds->ldims));
  double dx[3] = { cmflds->dx[0], cmflds->dx[1], cmflds->dx[2] };
  mrc_json_object_push(json, "dx", mrc_json_double_array_new(3, dx));

  mrc_json_t json_flds = mrc_json_object_new(2);
  mrc_json_object_push(json, "flds", json_flds);
  mrc_json_object_push_boolean(json_flds, "__field5d__", true);
  mrc_json_t json_flds_patches = mrc_json_array_new(cmflds->n_patches);
  mrc_json_object_push(json_flds, "data", json_flds_patches);

  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  for (int p = 0; p < cmflds->n_patches; p++) {
    cuda_mfields_copy_from_device(cmflds, p, flds, 0, cmflds->n_fields);

    mrc_json_t json_flds_comps = mrc_json_array_new(cmflds->n_fields);
    mrc_json_array_push(json_flds_patches, json_flds_comps);
    for (int m = 0; m < cmflds->n_fields; m++) {
      mrc_json_t json_fld_z = mrc_json_array_new(cmflds->im[2]);
      mrc_json_array_push(json_flds_comps, json_fld_z);
      for (int k = cmflds->ib[2]; k < cmflds->ib[2] + cmflds->im[2]; k++) {
	mrc_json_t json_fld_y = mrc_json_array_new(cmflds->im[1]);
	mrc_json_array_push(json_fld_z, json_fld_y);
	for (int j = cmflds->ib[1]; j < cmflds->ib[1] + cmflds->im[1]; j++) {
	  mrc_json_t json_fld_x = mrc_json_array_new(cmflds->im[0]);
	  mrc_json_array_push(json_fld_y, json_fld_x);
	  for (int i = cmflds->ib[0]; i < cmflds->ib[0] + cmflds->im[0]; i++) {
	    mrc_json_array_push_double(json_fld_x, _F3_S(flds, m, i,j,k));
	  }
	}
      }
    }
  }
  fields_single_t_dtor(&flds);

  return json;
}

// ----------------------------------------------------------------------
// cuda_mfields_dump

void
cuda_mfields_dump(struct cuda_mfields *cmflds, const char *filename)
{
  mrc_json_t json = cuda_mfields_to_json(cmflds);

  const char *buf = mrc_json_to_string(json);
  if (filename) {
    FILE *file = fopen(filename, "w");
    assert(file);
    fwrite(buf, 1, strlen(buf), file);
    fclose(file);
  } else {
    printf("cuda_mfields (json):\n%s\n", buf);
  }
  free((void *) buf);

  // FIXME free json
}

// ----------------------------------------------------------------------
// cuda_mfields_get_host_fields

fields_single_t
cuda_mfields_get_host_fields(struct cuda_mfields *cmflds)
{
  return fields_single_t_ctor(cmflds->ib, cmflds->im, cmflds->n_fields);
}

// ----------------------------------------------------------------------
// cuda_mfields_copy_to_device

void
cuda_mfields_copy_to_device(struct cuda_mfields *cmflds, int p, fields_single_t h_flds, int mb, int me)
{
  cudaError_t ierr;
  
  if (mb == me) {
    return;
  }
  assert(mb < me);

  unsigned int size = cmflds->n_cells_per_patch;
  ierr = cudaMemcpy(cmflds->d_flds_by_patch[p] + mb * size,
		    h_flds.data + mb * size,
		    (me - mb) * size * sizeof(float),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mfields_copy_from_device

void
cuda_mfields_copy_from_device(struct cuda_mfields *cmflds, int p, fields_single_t h_flds, int mb, int me)
{
  cudaError_t ierr;

  if (mb == me) {
    return;
  }
  assert(mb < me);

  unsigned int size = cmflds->n_cells_per_patch;
  ierr = cudaMemcpy(h_flds.data + mb * size,
		    cmflds->d_flds_by_patch[p] + mb * size,
		    (me - mb) * size * sizeof(float),
		    cudaMemcpyDeviceToHost); cudaCheck(ierr);
}

#define BND (2) // FIXME
#define X3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *mz + ((jz)+2))							\
    *my + ((jy)+2))							\
   *1 + (0))

#define F3_DDEV(d_flds, fldnr,ix,jy,jz)		\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]

#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

// ----------------------------------------------------------------------
// cuda_mfields_axpy_comp_yz

__global__ static void
axpy_comp_yz(float *y_flds, int ym, float a, float *x_flds, int xm,
	     int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= my || iz >= mz) {
    return;
  }

  iy -= BND;
  iz -= BND;

  F3_DDEV(y_flds, ym, 0,iy,iz) += a * F3_DDEV(x_flds, xm, 0,iy,iz);
}

void
cuda_mfields_axpy_comp_yz(struct cuda_mfields *cmflds_y, int ym, float a,
			  struct cuda_mfields *cmflds_x, int xm)
{
  int my = cmflds_y->im[1];
  int mz = cmflds_y->im[2];

  dim3 dimGrid((cmflds_y->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (cmflds_y->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);

  for (int p = 0; p < cmflds_y->n_patches; p++) {
    axpy_comp_yz<<<dimGrid, dimBlock>>>(cmflds_y->d_flds_by_patch[p], ym, a,
					cmflds_x->d_flds_by_patch[p], xm, my, mz);
  }
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// cuda_mfields_zero_comp_yz

// FIXME, this should be more easily doable by just a cudaMemset()

__global__ static void
zero_comp_yz(float *x_flds, int xm, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= my || iz >= mz) {
    return;
  }

  iy -= BND;
  iz -= BND;

  F3_DDEV(x_flds, xm, 0,iy,iz) = 0.f;
}

void
cuda_mfields_zero_comp_yz(struct cuda_mfields *cmflds, int xm)
{
  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  dim3 dimGrid((cmflds->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (cmflds->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);

  // OPT, should be done in a single kernel
  for (int p = 0; p < cmflds->n_patches; p++) {
    zero_comp_yz<<<dimGrid, dimBlock>>>(cmflds->d_flds_by_patch[p], xm,
					my, mz);
  }
  cuda_sync_if_enabled();
}

