
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include <mrc_profile.h>

// FIXME, hardcoding is bad, needs to be consistent, etc...

#define BND (2)

#define THREADS_PER_BLOCK 256

// ======================================================================
// cuda_mfields_bnd

// ----------------------------------------------------------------------
// cuda_mfields_bnd_create

struct cuda_mfields_bnd *
cuda_mfields_bnd_create()
{
  return (struct cuda_mfields_bnd *) calloc(1, sizeof(struct cuda_mfields_bnd));
}

// ----------------------------------------------------------------------
// cuda_mfields_bnd_destroy

void
cuda_mfields_bnd_destroy(struct cuda_mfields_bnd *cbnd)
{
  free(cbnd);
}

// ----------------------------------------------------------------------
// cuda_mfields_bnd_setup_d_nei_patch

static void
cuda_mfields_bnd_setup_d_nei_patch(struct cuda_mfields_bnd *cbnd, int n_recv_entries,
				   struct cuda_mfields_bnd_entry *recv_entry)
{
  cbnd->h_nei_patch = new int[9 * cbnd->n_patches];

  for (int p = 0; p < cbnd->n_patches; p++) {
    for (int dir1 = 0; dir1 < 9; dir1++) {
      cbnd->h_nei_patch[p * 9 + dir1] = -1;
    }
  }
  for (int i = 0; i < n_recv_entries; i++) {
    cbnd->h_nei_patch[recv_entry[i].patch * 9 + recv_entry[i].dir1 / 3] = recv_entry[i].nei_patch;
  }

  cudaError_t ierr;

  ierr = cudaMalloc((void **) &cbnd->d_nei_patch,
		    9 * cbnd->n_patches * sizeof(*cbnd->d_nei_patch)); cudaCheck(ierr);
  ierr = cudaMemcpy(cbnd->d_nei_patch, cbnd->h_nei_patch, 
		    9 * cbnd->n_patches * sizeof(*cbnd->d_nei_patch),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mfields_bnd_ctor

void
cuda_mfields_bnd_ctor(struct cuda_mfields_bnd *cbnd, struct cuda_mfields_bnd_params *prm)
{
  cbnd->n_patches = prm->n_patches;
  for (int d = 0; d < 3; d++) {
    cbnd->im[d] = prm->im[d];
    cbnd->ib[d] = prm->ib[d];
  }
  
  uint buf_size = 0;
  for (int p = 0; p < cbnd->n_patches; p++) {
    int B = 2*BND;
    if (cbnd->im[0] == 1) {// + 2*BND) {
      buf_size = 2*B * (cbnd->im[1] + cbnd->im[2] - 2*B);
      assert(buf_size ==  cbnd->im[1] * cbnd->im[2] -
	     (cbnd->im[1] - 2*B) * (cbnd->im[2] - 2*B)); 
    } else {
      buf_size = cbnd->im[0] * cbnd->im[1] * cbnd->im[2] -
	(cbnd->im[0] - 2*B) * (cbnd->im[1] - 2*B) * (cbnd->im[2] - 2*B); 
    }
  }
  cbnd->d_buf.resize(MAX_BND_COMPONENTS * cbnd->n_patches * buf_size);
  cbnd->h_buf.resize(MAX_BND_COMPONENTS * cbnd->n_patches * buf_size);

  cbnd->bnd_by_patch = new cuda_mfields_bnd_patch[cbnd->n_patches];
  for (int p = 0; p < cbnd->n_patches; p++) {
    struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
    int sz = 1;
    for (int d = 0; d < 3; d++) {
      cf->im[d] = cbnd->im[d];
      cf->ib[d] = cbnd->ib[d];
      sz *= cf->im[d];
    }
    cf->arr = new float [MAX_BND_COMPONENTS * sz];
    cf->arr_off = cf->arr 
      - ((cf->ib[2] * cf->im[1] + cf->ib[1]) * cf->im[0] + cf->ib[0]);
  }

  cuda_mfields_bnd_setup_d_nei_patch(cbnd, prm->n_recv_entries, prm->recv_entry);
}

// ----------------------------------------------------------------------
// cuda_mfields_bnd_dtor

void
cuda_mfields_bnd_dtor(struct cuda_mfields_bnd *cbnd)
{
  cudaError_t ierr;

  ierr = cudaFree(cbnd->d_nei_patch); cudaCheck(ierr);
  cbnd->h_buf.clear(); // FIXME, if we just called an actual dtor...
  cbnd->d_buf.clear(); // FIXME, if we just called an actual dtor...

  for (int n = 0; n < MAX_BND_FIELDS; n++) {
    struct cuda_mfields_bnd_map *map = &cbnd->map[n];
    ierr = cudaFree(map->d_map_out); cudaCheck(ierr);
    ierr = cudaFree(map->d_map_in); cudaCheck(ierr);
    delete[] map->h_map_out;
    delete[] map->h_map_in;
  }

  for (int p = 0; p < cbnd->n_patches; p++) {
    struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
    delete[] cf->arr;
  }
  delete[] cbnd->bnd_by_patch;
}

// ----------------------------------------------------------------------
// cuda_mfields_bnd_get_patch

struct cuda_mfields_bnd_patch *
cuda_mfields_bnd_get_patch(struct cuda_mfields_bnd *cbnd, int p)
{
  return &cbnd->bnd_by_patch[p];
}

// ----------------------------------------------------------------------
// cuda_mfields_bnd_get_map

static void cuda_mfields_bnd_setup_map(struct cuda_mfields_bnd *cbnd, int n_fields,
				       struct cuda_mfields_bnd_map *map);

static struct cuda_mfields_bnd_map *
cuda_mfields_bnd_get_map(struct cuda_mfields_bnd *cbnd, int n_fields)
{
  assert(n_fields - 1 < MAX_BND_FIELDS);
  struct cuda_mfields_bnd_map *map = &cbnd->map[n_fields - 1];

  if (!map->h_map_out) {
    cuda_mfields_bnd_setup_map(cbnd, n_fields, map);
  }
  
  return map;
}

// ======================================================================

enum {
  PACK,
  UNPACK,
};

// ======================================================================
// fields_device_pack

// FIXME/OPT: can probably be accelerated by making component the fast index

template<int B, int WHAT, int NR_COMPONENTS>
__global__ static void
k_fields_device_pack_yz(float *d_buf, DMFields d_flds, int mb, int gmy, int gmz,
			int nr_patches)
{
  uint buf_size = 2*B * (gmy + gmz - 2*B);
  int jx = 0;
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int n_threads = NR_COMPONENTS * buf_size;
  int p = tid / n_threads;
  if (p >= nr_patches)
    return;

  int n = tid - p * n_threads;
  int m = n / buf_size; n -= m * buf_size;
  int jz, jy;
  if (n < B * gmy) {
    jz = n / gmy; n -= jz * gmy;
    jy = n;
  } else if (n < B * gmy + (gmz - 2*B) * 2*B) {
    n -= B * gmy;
    jz = n / (2*B); n -= jz * 2*B;
    if (n < B) {
      jy = n;
    } else {
      jy = n + gmy - 2*B;
    }
    jz += B;
  } else {
    n -= B * gmy + (gmz - 2*B) * 2*B;
    jz = n / gmy; n -= jz * gmy;
    jy = n;
    jz += gmz - B;
  }
  
  if (WHAT == PACK) {
    d_buf[tid] = d_flds[p](mb + m, jx,jy-BND,jz-BND);
  } else if (WHAT == UNPACK) {
    d_flds[p](mb + m, jx,jy-BND,jz-BND) = d_buf[tid]; 
  }
}

template<int B, int WHAT, int NR_COMPONENTS>
__global__ static void
k_fields_device_pack2_yz(float *d_buf, DMFields d_flds, int mb, int *d_nei_patch_by_dir1,
			 int gmy, int gmz, int nr_patches, int nr_fields)
{
  uint nr_ghosts = 2*B * (gmy + gmz - 2*B);
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int p = tid / nr_ghosts;
  if (p >= nr_patches)
    return;

  int n = tid - p * nr_ghosts;
  int jy, jz;
  int diry, dirz;
  if (n < 2*B * gmy) {
    jz = n / gmy;
    jy = n - jz * gmy;
    if (jy < B) {
      diry = -1;
    } else if (jy < gmy - B) {
      diry = 0;
    } else {
      diry = 1;
    }
    if (jz < B) {
      dirz = -1;
    } else {
      dirz = 1;
      jz += gmz - 2*B;
    }
  } else {
    n -= 2*B * gmy;
    jz = n / (2*B) + B;
    jy = n % (2*B);
    dirz = 0;
    if (jy < B) {
      diry = -1;
    } else {
      diry = 1;
      jy += gmy - 2*B;
    }
  }

  int s_p = d_nei_patch_by_dir1[p*9 + 3*dirz + diry + 4];
  // copy only ghost areas that interface with remote patches
  if (1||s_p < 0) {
    for (int m = 0; m < NR_COMPONENTS; m++) {
      if (WHAT == PACK) {
	d_buf[m * nr_ghosts * nr_patches + tid] = d_flds[p](m + mb, 0,jy,jz);
      } else if (WHAT == UNPACK) {
	d_flds[p](m + mb, 0,jy,jz) = d_buf[m * nr_ghosts * nr_patches + tid]; 
      }
    }
  }
}

template<int B, int NR_COMPONENTS>
__global__ static void
k_fill_ghosts_local_yz(float *d_flds, int *d_nei_patch_by_dir1,
		       int gmy, int gmz, int nr_fields, int nr_patches)
{
  int nr_ghosts = 2*B * (gmy + gmz - 2*B);
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int r_p = tid / nr_ghosts;
  if (r_p >= nr_patches)
    return;

  int n = tid - r_p * nr_ghosts;

  int jy, jz;
  int diry, dirz;
  if (n < 2*B * gmy) {
    jz = n / gmy;
    jy = n - jz * gmy;
    if (jy < B) {
      diry = -1;
    } else if (jy < gmy - B) {
      diry = 0;
    } else {
      diry = 1;
    }
    if (jz < B) {
      dirz = -1;
    } else {
      dirz = 1;
      jz += gmz - 2*B;
    }
  } else {
    n -= 2*B * gmy;
    jz = n / (2*B) + B;
    jy = n % (2*B);
    dirz = 0;
    if (jy < B) {
      diry = -1;
    } else {
      diry = 1;
      jy += gmy - 2*B;
    }
  }
  int s_p = d_nei_patch_by_dir1[r_p*9 + 3*dirz + diry + 4];
  if (s_p >= 0) {
    float *r_f = &d_flds[((r_p * nr_fields) * gmz) * gmy];
    float *s_f = &d_flds[((s_p * nr_fields)
			  * gmz - dirz * (gmz - 2*2)) 
			 * gmy - diry * (gmy - 2*2)];
    for (int m = 0; m < NR_COMPONENTS; m++) {
      int i = (m * gmz + jz) * gmy + jy;
      r_f[i] = s_f[i];
    }
  }
}

void
cuda_fill_ghosts_local_gold(float *h_flds, int *nei_patch_by_dir1, int mb, int me,
			    int *im, int nr_fields, int nr_patches, int nr_ghosts, int threadIdx)
{
  int iy, iz;
  int diry, dirz;

  int r_p = threadIdx / nr_ghosts;
  if (r_p >= nr_patches)
    return;
  int tid = threadIdx - nr_ghosts * r_p;

  if (tid < 4 * im[1]) {
    iy = tid % im[1];
    iz = tid / im[1];
    if (iy < 2) {
      diry = -1;
    } else if (iy < im[1] - 2) {
	diry = 0;
    } else {
      diry = 1;
    }
    if (iz < 2) {
      dirz = -1;
    } else {
      dirz = 1;
      iz += im[2] - 2*2;
    }
  } else {
    int tid2 = tid - 4 * im[1];
    iy = tid2 % 4;
    iz = tid2 / 4 + 2;
    dirz = 0;
    if (iy < 2) {
      diry = -1;
    } else {
      diry = 1;
      iy += im[1] - 2*2;
    }
  }
  int s_p = nei_patch_by_dir1[r_p*9 + 3*dirz + diry + 4];
  if (s_p >= 0) {
    float *r_f = &h_flds[((r_p * nr_fields + mb) * im[2]) * im[1]];
    float *s_f = &h_flds[((s_p * nr_fields + mb)
			  * im[2] - dirz * (im[2] - 2*2)) 
			 * im[1] - diry * (im[1] - 2*2)];
    for (int m = 0; m < me - mb; m++) {
      int i = (m * im[2] + iz) * im[1] + iy;
      r_f[i] = s_f[i];
    }
  }
}

template<int B, bool pack>
static void
fields_device_pack_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
		      int mb, int me)
{
  uint size = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
  int gmy = cmflds->im[1], gmz = cmflds->im[2];
  uint buf_size = 2*B * (gmy + gmz - 2*B);
  int n_threads = buf_size * (me - mb) * cmflds->n_patches;

  dim3 dimGrid((n_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_buf = cbnd->d_buf.data().get();
  if (me - mb == 3) {
    k_fields_device_pack_yz<B, pack, 3> <<<dimGrid, dimBlock>>>
      (d_buf, *cmflds, mb, gmy, gmz, cmflds->n_patches);
  } else if (me - mb == 2) {
    k_fields_device_pack_yz<B, pack, 2> <<<dimGrid, dimBlock>>>
      (d_buf, *cmflds, mb, gmy, gmz, cmflds->n_patches);
  } else if (me - mb == 1) {
    k_fields_device_pack_yz<B, pack, 1> <<<dimGrid, dimBlock>>>
      (d_buf, *cmflds, mb, gmy, gmz, cmflds->n_patches);
  } else {
    printf("mb %d me %d\n", mb, me);
    assert(0);
  }
  cuda_sync_if_enabled();
}

template<int B, bool pack>
static void
fields_device_pack2_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
		       int mb, int me)
{
  const int NR_COMPONENTS = 3;
  assert(me - mb == NR_COMPONENTS);

  int *im = cmflds->im;
  assert(im[0] == 1);
  int n_fields = cmflds->n_fields;
  int n_patches = cmflds->n_patches;
  int n_ghosts = 2*B * (im[1] + im[2] - 2*B);
  int n_threads = n_ghosts * n_patches;

  dim3 dimGrid((n_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  k_fields_device_pack2_yz<B, pack, NR_COMPONENTS> <<<dimGrid, dimBlock>>>
    (cbnd->d_buf.data().get(), *cmflds, mb, cbnd->d_nei_patch,
     im[1], im[2], n_patches, n_fields);
  cuda_sync_if_enabled();
}

void
cuda_mfields_bnd_fill_ghosts_local_xyz(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
				       int mb, int me)
{
  const int B = 2;
  int *im = cmflds->im;
  int n_fields = cmflds->n_fields;
  int n_patches = cmflds->n_patches;
  int n_ghosts = im[0] * im[1] * im[2] -
    (im[0] - 2*B) * (im[1] - 2*B) * (im[2] - 2*B);
  int n_threads = n_ghosts * n_patches;

#if 0
  dim3 dimGrid((n_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_flds = cmflds->data() + mb * im[1] * im[2];
  if (me - mb == 3) {
    k_fill_ghosts_local_yz<B, 3> <<<dimGrid, dimBlock>>>
      (d_flds, cbnd->d_nei_patch, im[1], im[2], n_fields, n_patches);
  } else if (me - mb == 1) {
    k_fill_ghosts_local_yz<B, 1> <<<dimGrid, dimBlock>>>
      (d_flds, cbnd->d_nei_patch, im[1], im[2], n_fields, n_patches);
  } else {
    assert(0);
  }
  cuda_sync_if_enabled();

#else

  thrust::device_ptr<float> d_flds(cmflds->data());
  thrust::host_vector<float> h_flds(d_flds, d_flds + n_patches * n_fields * im[0] * im[1] * im[2]);

  for (int tid = 0; tid < n_threads; tid++) {
    cuda_fill_ghosts_local_gold(&h_flds[0], cbnd->d_nei_patch, mb, me, im, n_fields, 
				n_patches, n_ghosts, tid);
  }
  
  thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
#endif
}

void
cuda_mfields_bnd_fill_ghosts_local(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
				   int mb, int me)
{
  const int B = 2;
  int *im = cmflds->im;
  if (im[0] > 1) {
    cuda_mfields_bnd_fill_ghosts_local_xyz(cbnd, cmflds, mb, me);
  }
  assert(im[0] == 1);
  int n_fields = cmflds->n_fields;
  int n_patches = cmflds->n_patches;
  int n_ghosts = 2*B * (im[1] + im[2] - 2*B);
  int n_threads = n_ghosts * n_patches;

#if 1
  dim3 dimGrid((n_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_flds = cmflds->data() + mb * im[1] * im[2];
  if (me - mb == 3) {
    k_fill_ghosts_local_yz<B, 3> <<<dimGrid, dimBlock>>>
      (d_flds, cbnd->d_nei_patch, im[1], im[2], n_fields, n_patches);
  } else if (me - mb == 1) {
    k_fill_ghosts_local_yz<B, 1> <<<dimGrid, dimBlock>>>
      (d_flds, cbnd->d_nei_patch, im[1], im[2], n_fields, n_patches);
  } else {
    assert(0);
  }
  cuda_sync_if_enabled();

#else

  thrust::device_ptr<float> d_flds(cmflds->data());
  thrust::host_vector<float> h_flds(d_flds, d_flds + n_patches * n_fields * im[1] * im[2]);

  for (int tid = 0; tid < n_threads; tid++) {
    cuda_fill_ghosts_local_gold(&h_flds[0], cbnd->d_nei_patch, mb, me, im, n_fields, 
				n_patches, n_ghosts, tid);
  }
  
  thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
#endif
}

#ifdef F3_CF_BOUNDS_CHECK
#define F3_CF_0_OFF(cf, fldnr, jx,jy,jz) ({				\
  assert(jx == 0); /* FIXME yz only! */				        \
  assert(jx >= 0 && jx < (cf)->im[0]);					\
  assert(jy >= 0 && jy < (cf)->im[1]);					\
  assert(jz >= 0 && jz < (cf)->im[2]);					\
  int __off = (((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx); \
  __off; })
#else
#define F3_CF_0_OFF(cf, fldnr, jx,jy,jz)				\
  ((((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx))
#endif

#define F3_CF_0(cf, fldnr, jx,jy,jz)					\
  ((cf)->arr[F3_CF_0_OFF(cf, fldnr, jx,jy,jz)])

// ======================================================================
// fields_host_pack

#define WHAT do {					\
    if (what == PACK) {					\
      h_buf[tid++] = F3_CF_0(cf, m, 0,jy,jz);		\
    } else if (what == UNPACK) {			\
      F3_CF_0(cf, m, 0,jy,jz) = h_buf[tid++];		\
    }							\
  } while(0)

template<int B, int what>
static void
fields_host_pack_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
		    int mb, int me)
{
  int gmy = cmflds->im[1], gmz = cmflds->im[2];
  uint buf_size = 2*B * (gmy + gmz - 2*B);

  for (int p = 0; p < cmflds->n_patches; p++) {
    struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
    float* h_buf = &cbnd->h_buf[p * buf_size * (me - mb)];
    
    int gmy = cf->im[1], gmz = cf->im[2];
    int tid = 0;
    for (int m = 0; m < me - mb; m++) {
      for (int jz = 0; jz < B; jz++) {
	for (int jy = 0; jy < gmy; jy++) {
	  WHAT;
	}
      }
      for (int jz = B; jz < gmz - B; jz++) {
	for (int jy = 0; jy < B; jy++) {
	  WHAT;
	}
	for (int jy = gmy - B; jy < gmy; jy++) {
	  WHAT;
	}
      }
      for (int jz = gmz - B; jz < gmz; jz++) {
	for (int jy = 0; jy < gmy; jy++) {
	  WHAT;
	}
      }
    }
  }
}

template<int B, int what>
static void
fields_host_pack2_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
		     int mb, int me)
{
  auto& h_buf = cbnd->h_buf;
  int tid = 0;
  for (int m = 0; m < me - mb; m++) {
    for (int p = 0; p < cmflds->n_patches; p++) {
      struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
      
      int gmy = cf->im[1], gmz = cf->im[2];
      for (int jz = 0; jz < B; jz++) {
	for (int jy = 0; jy < gmy; jy++) {
	  WHAT;
	}
      }
      for (int jz = gmz - B; jz < gmz; jz++) {
	for (int jy = 0; jy < gmy; jy++) {
	  WHAT;
	}
      }
      for (int jz = B; jz < gmz - B; jz++) {
	for (int jy = 0; jy < B; jy++) {
	  WHAT;
	}
	for (int jy = gmy - B; jy < gmy; jy++) {
	  WHAT;
	}
      }
    }
  }
}

#undef WHAT

template<int B, int what>
static void
fields_host_pack3_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
		     int mb, int me)
{
  int *im = cmflds->im;
  int n_fields = cmflds->n_fields;

  struct cuda_mfields_bnd_map *map = cuda_mfields_bnd_get_map(cbnd, n_fields);

  auto& h_buf = cbnd->h_buf;
  int nr_map;
  int *h_map;
  if (B == 2) {
    h_map = map->h_map_out;
    nr_map = map->nr_map_out;
  } else if (B == 4) {
    h_map = map->h_map_in;
    nr_map = map->nr_map_in;
  } else {
    assert(0);
  }
  for (int tid = 0; tid < nr_map; tid++) {
    int i = h_map[tid];
    int p = i / (n_fields * im[2] * im[1]);
    int off = i - p * (n_fields * im[2] * im[1]);
    struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
    for (int m = 0; m < me - mb; m++) {
      if (what == PACK) {
	h_buf[tid + m * nr_map] = cf->arr[off + m * im[2] * im[1]];
      } else if (what == UNPACK) {
	cf->arr[off + m * im[2] * im[1]] = h_buf[tid + m * nr_map];
      }
    }
  }
}

template<int B>
static void
__fields_cuda_from_device_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
			     int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_device_pack", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy", 1., 0, 0);
    pr3 = prof_register("field_host_unpack", 1., 0, 0);
  }
    
  int gmy = cmflds->im[1], gmz = cmflds->im[2];
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(cmflds->ib[1] == -BND);
  assert(cmflds->im[1] >= 2 * B);
  assert(cmflds->im[2] >= 2 * B);

  prof_start(pr1);
  fields_device_pack_yz<B, PACK>(cmflds, cbnd, mb, me);
  prof_stop(pr1);
  
  prof_start(pr2);
  cbnd->h_buf = cbnd->d_buf; // copy to host
  prof_stop(pr2);

  prof_start(pr3);
  fields_host_pack_yz<B, UNPACK>(cmflds, cbnd, mb, me);
  prof_stop(pr3);
}

template<int B>
static void
__fields_cuda_to_device_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
			   int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_host_pack", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy", 1., 0, 0);
    pr3 = prof_register("field_device_unpack", 1., 0, 0);
  }
  
  int gmy = cmflds->im[1], gmz = cmflds->im[2];
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(cmflds->ib[1] == -BND);
  assert(cmflds->im[1] >= 2 * B);
  assert(cmflds->im[2] >= 2 * B);

  prof_start(pr1);
  fields_host_pack_yz<B, PACK>(cmflds, cbnd, mb, me);
  prof_stop(pr1);

  prof_start(pr2);
  thrust::copy(cbnd->h_buf.begin(), cbnd->h_buf.end(),
	       cbnd->d_buf.begin());
  prof_stop(pr2);

  prof_start(pr3);
  fields_device_pack_yz<B, UNPACK>(cmflds, cbnd, mb, me);
  prof_stop(pr3);
}

template<int B, bool pack>
static void fields_device_pack3_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
				   int mb, int me);

template<int B>
static void
__fields_cuda_from_device3_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
			      int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_device_pack 3i", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy 3i", 1., 0, 0);
    pr3 = prof_register("field_host_unpack 3i", 1., 0, 0);
  }
    
  struct cuda_mfields_bnd_map *map = cuda_mfields_bnd_get_map(cbnd, cmflds->n_fields);

  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(cmflds->ib[1] == -BND);
  assert(cmflds->im[1] >= 2 * B);
  assert(cmflds->im[2] >= 2 * B);

  int nr_map;
  if (B == 4) {
    nr_map = map->nr_map_in;
  } else {
    assert(0);
  }

  prof_start(pr1);
  fields_device_pack3_yz<B, PACK>(cmflds, cbnd, mb, me);
  prof_stop(pr1);
  
  prof_start(pr2);
  assert(B == 4);
  thrust::copy(cbnd->d_buf.begin(), cbnd->d_buf.begin() + (me - mb) * nr_map, cbnd->h_buf.begin());
  prof_stop(pr2);

  prof_start(pr3);
  fields_host_pack3_yz<B, UNPACK>(cmflds, cbnd, mb, me);
  prof_stop(pr3);
}

template<int B>
static void
__fields_cuda_to_device3_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
			    int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_host_pack 3o", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy 3o", 1., 0, 0);
    pr3 = prof_register("field_device_unpack 3o", 1., 0, 0);
  }

  struct cuda_mfields_bnd_map *map = cuda_mfields_bnd_get_map(cbnd, cmflds->n_fields);

  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(cmflds->ib[1] == -BND);
  assert(cmflds->im[1] >= 2 * B);
  assert(cmflds->im[2] >= 2 * B);

  int nr_map;
  if (B == 2) {
    nr_map = map->nr_map_out;
  } else {
    assert(0);
  }

  prof_start(pr1);
  fields_host_pack3_yz<B, PACK>(cmflds, cbnd, mb, me);
  prof_stop(pr1);

  prof_start(pr2);
  thrust::copy(cbnd->h_buf.begin(), cbnd->h_buf.begin() + (me - mb) * nr_map,
	       cbnd->d_buf.begin());
  prof_stop(pr2);

  prof_start(pr3);
  fields_device_pack3_yz<B, UNPACK>(cmflds, cbnd, mb, me);
  prof_stop(pr3);
}

// ======================================================================

void
cuda_mfields_bnd_from_device_inside(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
				    int mb, int me)
{
  if (cmflds->im[0] == 2 * -cmflds->ib[0] + 1) {
    __fields_cuda_from_device_yz<2*BND>(cmflds, cbnd, mb, me);
  } else {
    for (int p = 0; p < cmflds->n_patches; p++) {
      struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
      uint size = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
      cudaError_t ierr;
      ierr = cudaMemcpy(cf->arr, (*cmflds)[p].data() + mb * size,
			(me - mb) * size * sizeof(float),
			cudaMemcpyDeviceToHost); cudaCheck(ierr);
    }
  }
}

void
cuda_mfields_bnd_from_device_inside_only(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
					 int mb, int me)
{
  if (cmflds->im[0] == 2 * -cmflds->ib[0] + 1) {
    __fields_cuda_from_device3_yz<2*BND>(cmflds, cbnd, mb, me);
  } else {
    for (int p = 0; p < cmflds->n_patches; p++) {
      struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
      uint size = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
      cudaError_t ierr;
      ierr = cudaMemcpy(cf->arr, (*cmflds)[p].data() + mb * size,
			(me - mb) * size * sizeof(float),
			cudaMemcpyDeviceToHost); cudaCheck(ierr);
    }
  }
}

void
cuda_mfields_bnd_to_device_outside(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
				   int mb, int me)
{
  if (cmflds->im[0] == 2 * -cmflds->ib[0] + 1) {
    __fields_cuda_to_device3_yz<BND>(cmflds, cbnd, mb, me);
  } else {
    for (int p = 0; p < cmflds->n_patches; p++) {
      struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
      uint size = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
      cudaError_t ierr;
      ierr = cudaMemcpy((*cmflds)[p].data() + mb * size,
			cf->arr,
			(me - mb) * size * sizeof(float),
			cudaMemcpyHostToDevice); cudaCheck(ierr);
    }
  }
}

void
cuda_mfields_bnd_to_device_inside(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
				int mb, int me)
{
  if (cmflds->im[0] == 2 * -cmflds->ib[0] + 1) {
    __fields_cuda_to_device_yz<2*BND>(cmflds, cbnd, mb, me);
  } else {
    for (int p = 0; p < cmflds->n_patches; p++) {
      struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
      uint size = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
      cudaError_t ierr;
      ierr = cudaMemcpy((*cmflds)[p].data() + mb * size,
			cf->arr,
			(me - mb) * size * sizeof(float),
			cudaMemcpyHostToDevice); cudaCheck(ierr);
    }
  }
}

template<int WHAT>
static void
g_fields_device_pack3_yz(int tid, float *d_buf, float *d_flds, int *d_map, int *d_nei_patch_by_dir1,
			 int gmy, int gmz, int nr_patches, int nr_components, int nr_map)
{
  if (tid >= nr_map)
    return;

  // copy only ghost areas that interface with remote patches
  int i = d_map[tid];
  for (int m = 0; m < nr_components; m++) {
    // FIXME, should use D_F3
    if (WHAT == PACK) {
      d_buf[m * nr_map + tid] = d_flds[i + m * gmz * gmy];
    } else if (WHAT == UNPACK) {
      d_flds[i + m * gmz * gmy] = d_buf[m * nr_map + tid]; 
    }
  }
}

template<int WHAT>
__global__ static void
k_fields_device_pack3_yz(float *d_buf, float *d_flds, int *d_map, int *d_nei_patch_by_dir1,
			 int gmy, int gmz, int nr_patches, int nr_components, int nr_map)
{
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  if (tid >= nr_map)
    return;

  // copy only ghost areas that interface with remote patches
  int i = d_map[tid];
  for (int m = 0; m < nr_components; m++) {
    // FIXME, should use D_F3
    if (WHAT == PACK) {
      d_buf[m * nr_map + tid] = d_flds[i + m * gmz * gmy];
    } else if (WHAT == UNPACK) {
      d_flds[i + m * gmz * gmy] = d_buf[m * nr_map + tid]; 
    }
  }
}

#undef check

#undef _GLIBCXX_USE_INT128

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

static void
fields_create_map_out_yz(struct cuda_mfields_bnd *cbnd, int n_fields,
			 int B, int *p_nr_map, int **p_h_map)
{
  bool remote_only = true;
  int *im = cbnd->im;
  int n_patches = cbnd->n_patches;

  int nr_map = 0;
  for (int p = 0; p < n_patches; p++) {
    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	if (dir[1] == 0 && dir[2] == 0) {
	  continue;
	}
	int s_p = cbnd->h_nei_patch[p*9 + 3*dir[2] + dir[1] + 4];
	if (!remote_only || s_p < 0) {
	  nr_map += ((dir[1] == 0 ? im[1] - 2*B : B) *
		     (dir[2] == 0 ? im[2] - 2*B : B));
	}
      }
    }
  }
  *p_nr_map = nr_map;
  *p_h_map = new int[nr_map];

  int tid = 0;
  for (int p = 0; p < n_patches; p++) {
    for (int jjz = 0; jjz < 2*B; jjz++) {
      int jz = jjz;
      int dirz;
      if (jz < B) {
	dirz = -1;
      } else {
	dirz = 1;
	jz += im[2] - 2*B;
      }
      for (int jy = 0; jy < im[1]; jy++) {
	int diry;
	if (jy < B) {
	  diry = -1;
	} else if (jy < im[1] - B) {
	  diry = 0;
	} else {
	  diry = 1;
	}
	int s_p = cbnd->h_nei_patch[p*9 + 3*dirz + diry + 4];
	if (!remote_only || s_p < 0) {
	  (*p_h_map)[tid++] = ((p * n_fields + 0) * im[2] + jz) * im[1] + jy;
	}
      }
    }
    for (int jz = B; jz < im[2] - B; jz++) {
      int dirz = 0;
      for (int jjy = 0; jjy < 2*B; jjy++) {
	int jy = jjy;
	int diry;
	if (jy < B) {
	  diry = -1;
	} else {
	  diry = 1;
	  jy += im[1] - 2*B;
	}
	int s_p = cbnd->h_nei_patch[p*9 + 3*dirz + diry + 4];
	if (!remote_only || s_p < 0) {
	  (*p_h_map)[tid++] = ((p * n_fields + 0) * im[2] + jz) * im[1] + jy;
	}
      }
    }
  }
  assert(tid == nr_map);
}

static void
fields_create_map_in_yz(struct cuda_mfields_bnd *cbnd, int n_fields,
			int B, int *p_nr_map, int **p_h_map)
{
  bool remote_only = true;
  int *im = cbnd->im;
  int ldims[3] = { 1, im[1] - 2*2, im[2] - 2*2 };
  int n_patches = cbnd->n_patches;
  bool *has_nei = new bool[9 * n_patches];
  
  // FIXME, we don't need the ghosts here...

  for (int p = 0; p < n_patches; p++) {
    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	int dir1 = p*9 + 3*dir[2] + dir[1] + 4;
	has_nei[dir1] = false;
	int p_mm = cbnd->h_nei_patch[p*9 + -1 + 3*-1 + 4];
	int p_0m = cbnd->h_nei_patch[p*9 +  0 + 3*-1 + 4];
	int p_pm = cbnd->h_nei_patch[p*9 + +1 + 3*-1 + 4];
	int p_m0 = cbnd->h_nei_patch[p*9 + -1 + 3* 0 + 4];
	int p_p0 = cbnd->h_nei_patch[p*9 + +1 + 3* 0 + 4];
	int p_mp = cbnd->h_nei_patch[p*9 + +1 + 3*+1 + 4];
	int p_0p = cbnd->h_nei_patch[p*9 +  0 + 3*+1 + 4];
	int p_pp = cbnd->h_nei_patch[p*9 + -1 + 3*+1 + 4];
	if (dir[1] == 0 && dir[2] == 0) {
	} else if (dir[1] == -1 && dir[2] == -1) {
	  if (p_mm < 0 || p_m0 < 0 || p_0m < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] ==  0 && dir[2] == -1) {
	  if (p_0m < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == +1 && dir[2] == -1) {
	  if (p_pm < 0 || p_0m < 0 || p_p0 < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == -1 && dir[2] ==  0) {
	  if (p_m0 < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == +1 && dir[2] ==  0) {
	  if (p_p0 < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == -1 && dir[2] == +1) {
	  if (p_mp < 0 || p_m0 < 0 || p_0p < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] ==  0 && dir[2] == +1) {
	  if (p_0p < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == +1 && dir[2] == +1) {
	  if (p_pp < 0 || p_0p < 0 || p_p0 < 0) {
	    has_nei[dir1] = true;
	  }
	}
      }
    }
  }

  int nr_map = 0;
  for (int p = 0; p < n_patches; p++) {
    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	if (dir[1] == 0 && dir[2] == 0) {
	  continue;
	}
	int dir1 = p*9 + 3*dir[2] + dir[1] + 4;
	if (!remote_only || has_nei[dir1]) {
	  nr_map += ((dir[1] == 0 ? ldims[1] - 2*B : B) *
		     (dir[2] == 0 ? ldims[2] - 2*B : B));
	}
      }
    }
  }
  *p_nr_map = nr_map;
  *p_h_map = new int[nr_map];

  int tid = 0;
  for (int p = 0; p < n_patches; p++) {
    for (int jjz = 0; jjz < 2*B; jjz++) {
      int jz = jjz;
      int dirz;
      if (jz < B) {
	dirz = -1;
      } else {
	dirz = 1;
	jz += ldims[2] - 2*B;
      }
      for (int jy = 0; jy < ldims[1]; jy++) {
	int diry;
	if (jy < B) {
	  diry = -1;
	} else if (jy < ldims[1] - B) {
	  diry = 0;
	} else {
	  diry = 1;
	}
	int dir1 = p*9 + 3*dirz + diry + 4;
	if (!remote_only || has_nei[dir1]) {
	  (*p_h_map)[tid++] = ((p * n_fields + 0) * im[2] + (jz + B)) * im[1] + (jy + B);
	}
      }
    }
    for (int jz = B; jz < ldims[2] - B; jz++) {
      int dirz = 0;
      for (int jjy = 0; jjy < 2*B; jjy++) {
	int jy = jjy;
	int diry;
	if (jy < B) {
	  diry = -1;
	} else {
	  diry = 1;
	  jy += ldims[1] - 2*B;
	}
	int dir1 = p*9 + 3*dirz + diry + 4;
	if (!remote_only || has_nei[dir1]) {
	  (*p_h_map)[tid++] = ((p * n_fields + 0) * im[2] + (jz + B)) * im[1] + (jy + B);
	}
      }
    }
  }
  assert(tid == nr_map);
  delete[] has_nei;
}

// ----------------------------------------------------------------------
// cuda_mfields_bnd_setup_map

void
cuda_mfields_bnd_setup_map(struct cuda_mfields_bnd *cbnd, int n_fields,
			   struct cuda_mfields_bnd_map *map)
{
  cudaError_t ierr;

  fields_create_map_out_yz(cbnd, n_fields, 2, &map->nr_map_out, &map->h_map_out);
  //printf("map_out %d\n", map->nr_map_out);
  
  ierr = cudaMalloc((void **) &map->d_map_out,
		    map->nr_map_out * sizeof(*map->d_map_out)); cudaCheck(ierr);
  ierr = cudaMemcpy(map->d_map_out, map->h_map_out, 
		    map->nr_map_out * sizeof(*map->d_map_out),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
  
  fields_create_map_in_yz(cbnd, n_fields, 2, &map->nr_map_in, &map->h_map_in);
  //printf("map_in %d\n", map->nr_map_in);
  
  ierr = cudaMalloc((void **) &map->d_map_in,
		    map->nr_map_in * sizeof(*map->d_map_in)); cudaCheck(ierr);
  ierr = cudaMemcpy(map->d_map_in, map->h_map_in, 
		    map->nr_map_in * sizeof(*map->d_map_in),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
}

template<int B, bool pack>
static void
fields_device_pack3_yz(struct cuda_mfields *cmflds, struct cuda_mfields_bnd *cbnd,
		       int mb, int me)
{
  struct cuda_mfields_bnd_map *map = cuda_mfields_bnd_get_map(cbnd, cmflds->n_fields);

  int *im = cmflds->im;
  assert(im[0] == 1);
  int n_patches = cmflds->n_patches;
  int n_map;
  int *d_map;
  if (B == 2) {
    n_map = map->nr_map_out;
    d_map = map->d_map_out;
  } else if (B == 4) {
    n_map = map->nr_map_in;
    d_map = map->d_map_in;
  } else {
    assert(0);
  }

  if (n_map == 0) {
    return;
  }
  
#if 1
  dim3 dimGrid((n_map + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_flds = cmflds->data() + mb * im[1] * im[2];
  if (me - mb == 3) {
    k_fields_device_pack3_yz<pack> <<<dimGrid, dimBlock>>>
      (cbnd->d_buf.data().get(), d_flds, d_map, cbnd->d_nei_patch,
       im[1], im[2], n_patches, 3, n_map);
  } else if (me - mb == 1) {
    k_fields_device_pack3_yz<pack> <<<dimGrid, dimBlock>>>
      (cbnd->d_buf.data().get(), d_flds, d_map, cbnd->d_nei_patch,
       im[1], im[2], n_patches, 1, n_map);
  } else {
    assert(0);
  }
  cuda_sync_if_enabled();
#else
  thrust::device_ptr<float> d_bnd_buf(cmflds->d_bnd_buf);
  thrust::device_ptr<float> d_flds(cmflds->d_flds);
  thrust::host_vector<float> h_bnd_buf(d_bnd_buf, d_bnd_buf + nr_patches * nr_ghosts * NR_COMPONENTS);
  thrust::host_vector<float> h_flds(d_flds, d_flds + nr_patches * nr_fields * im[1] * im[2]);

  for (int tid = 0; tid < nr_map; tid++) {
    g_fields_device_pack3_yz<pack>
      (tid, &h_bnd_buf[0], &h_flds[mb * im[1] * im[2]], &h_map[0], &h_nei_patch_by_dir1[0],
       im[1], im[2], nr_patches, NR_COMPONENTS, nr_map);
  }

  thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
#endif
}

