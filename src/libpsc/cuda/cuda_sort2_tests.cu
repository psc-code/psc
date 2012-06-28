
template <int RADIX_BITS>
struct my_sort {
  thrust::device_vector<K> &d_keys;
  thrust::device_vector<int> &d_offsets;
  
  int nr_blocks;
  int spine_elements;
  thrust::device_vector<K> d_alt_keys;
  thrust::device_vector<V> d_values;
  thrust::device_vector<V> d_alt_values;
  thrust::device_vector<int> d_spine;

  bool debug;

  static int get_spine_elements(int nr_blocks) {
    int spine_cycle_elements = B40C_RADIXSORT_SPINE_CYCLE_ELEMENTS;
    int spine_cycles = (nr_blocks * (1 << RADIX_BITS) + spine_cycle_elements - 1) /
      spine_cycle_elements;
    int spine_elements = spine_cycles * B40C_RADIXSORT_SPINE_CYCLE_ELEMENTS;
    return spine_elements;
  }
  
  my_sort(thrust::device_vector<K> &d_keys_, thrust::device_vector<int> &d_offsets_) :
    d_keys(d_keys_),
    d_offsets(d_offsets_),
    nr_blocks(d_offsets.size() - 1),
    spine_elements(get_spine_elements(nr_blocks)),
    d_alt_keys(d_keys.size()),
    d_values(d_keys.size()),
    d_alt_values(d_keys.size()),
    d_spine(spine_elements)
  {
    thrust::sequence(d_values.begin(), d_values.end());
  }

  void sort() {
    const int threads = B40C_RADIXSORT_THREADS;
    const int RADIX_DIGITS = 1 << RADIX_BITS;

    RakingReduction2<K, V, 0, RADIX_BITS, 0, NopFunctor<K> > <<<nr_blocks, threads>>>
      (thrust::raw_pointer_cast(&d_spine[0]),
       thrust::raw_pointer_cast(&d_keys[0]),
       thrust::raw_pointer_cast(&d_offsets[0]));
    
    if (debug) {
      thrust::host_vector<int> h_spine(&d_spine[0], &d_spine[nr_blocks * RADIX_DIGITS]);
      print("h_spine", h_spine);
    }
    
    SrtsScanSpine<void><<<nr_blocks, B40C_RADIXSORT_SPINE_THREADS>>>
      (thrust::raw_pointer_cast(&d_spine[0]),
       thrust::raw_pointer_cast(&d_spine[0]),
       spine_elements);
    
    if (debug) {
      thrust::host_vector<int> h_spine(&d_spine[0], &d_spine[nr_blocks * RADIX_DIGITS]);
      print("h_spine scan", h_spine);
    }
    
    ScanScatterDigits2<K, V, 0, RADIX_BITS, 0, NopFunctor<K>, NopFunctor<K> > <<<nr_blocks, threads>>>
      (thrust::raw_pointer_cast(&d_spine[0]),
       thrust::raw_pointer_cast(&d_keys[0]),
       thrust::raw_pointer_cast(&d_alt_keys[0]),
       thrust::raw_pointer_cast(&d_values[0]),
       thrust::raw_pointer_cast(&d_alt_values[0]),
       thrust::raw_pointer_cast(&d_offsets[0]));

    if (debug) {
      thrust::host_vector<K> h_keys = d_keys;
      print("h_keys  ", h_keys);
      thrust::host_vector<K> h_values = d_values;
      print("h_values", h_values);
      thrust::host_vector<K> h_alt_keys = d_alt_keys;
      print("h_alt_keys  ", h_alt_keys);
      thrust::host_vector<K> h_alt_values = d_alt_values;
      print("h_alt_values", h_alt_values); 
      printf("\n");
    }
  }
};
  
// blockIdx to rel offset
struct PreShiftFunctor {
  __device__ __host__ __forceinline__ void operator()(unsigned int &converted_key) {
    converted_key -= blockIdx.x - 1;
  }
  __device__ __host__ __forceinline__ static bool MustApply(){ return true;}
};

// rel offset to blockIdx
struct PostShiftFunctor {
  __device__ __host__ __forceinline__ void operator()(unsigned int &converted_key) {
    converted_key += blockIdx.x - 1;
  }
  __device__ __host__ __forceinline__ static bool MustApply(){ return true;}
};

template <int RADIX_BITS>
struct my_sort_b {
  thrust::device_vector<K> &d_keys;
  thrust::device_vector<int> &d_offsets;
  
  int nr_blocks;
  int spine_elements;
  thrust::device_vector<K> d_alt_keys;
  thrust::device_vector<V> d_values;
  thrust::device_vector<V> d_alt_values;
  thrust::device_vector<int> d_spine;

  bool debug;

  static int get_spine_elements(int nr_blocks) {
    int spine_cycle_elements = B40C_RADIXSORT_SPINE_CYCLE_ELEMENTS;
    int spine_cycles = (nr_blocks * (1 << RADIX_BITS) + spine_cycle_elements - 1) /
      spine_cycle_elements;
    int spine_elements = spine_cycles * B40C_RADIXSORT_SPINE_CYCLE_ELEMENTS;
    return spine_elements;
  }
  
  my_sort_b(thrust::device_vector<K> &d_keys_, thrust::device_vector<int> &d_offsets_) :
    d_keys(d_keys_),
    d_offsets(d_offsets_),
    nr_blocks(d_offsets.size() - 1),
    spine_elements(get_spine_elements(nr_blocks)),
    d_alt_keys(d_keys.size()),
    d_values(d_keys.size()),
    d_alt_values(d_keys.size()),
    d_spine(spine_elements),
    debug(false)
  {
    thrust::sequence(d_values.begin(), d_values.end());
  }

  // 1d +-1 only
  void top_scan_1d()
  {
    const int RADIX_DIGITS = 1 << RADIX_BITS;
    thrust::host_vector<int> h_spine = d_spine;

    int sum = 0;
    for (int b = 0; b < nr_blocks + RADIX_DIGITS; b++) {
      for (int row = RADIX_DIGITS - 1; row >= 0; row--) {
	int col = b - row;
	if (col >= 0 && col < nr_blocks) {
	  int val = h_spine[row * nr_blocks + col];
	  h_spine[row * nr_blocks + col] = sum;
	  sum += val;
	}
      }
    }
    d_spine = h_spine;
  }

  void top_scan()
  {
    const int RADIX_DIGITS = 1 << RADIX_BITS;
    thrust::host_vector<int> h_spine = d_spine;

    thrust::host_vector<int> h_matrix(nr_blocks * nr_blocks);

    for (int b = 0; b < nr_blocks; b++) {
      for (int d = 0; d < RADIX_DIGITS; d++) {
	int val = h_spine[d * nr_blocks + b];
	if (!val)
	  continue;

	int b2 = dir1_to_block_idx(b, d);
	h_matrix[b2 * nr_blocks + b] = val;
      }
    }

    if (debug) {
      print_matrix(h_matrix, nr_blocks);
    }

    int sum = 0;
    for (int b2 = 0; b2 < nr_blocks; b2++) {
      for (int b = 0; b < nr_blocks; b++) {
	int val = h_matrix[b2 * nr_blocks + b];
	h_matrix[b2 * nr_blocks + b] = sum;
	sum += val;
      }
    }

    if (debug) {
      print_matrix(h_matrix, nr_blocks);
    }

    // can skip the ones where h_spine == 0 (?)
    for (int b = 0; b < nr_blocks; b++) {
      for (int b2 = 0; b2 < nr_blocks; b2++) {
	int d = block_idx_to_dir1_oob(b, b2);
	if (d >= 0) {
	  int val = h_matrix[b2 * nr_blocks + b];
	  h_spine[d * nr_blocks + b] = val;
	}
      }
    }
    d_spine = h_spine;
  }

  void sort()
  {
    const int threads = B40C_RADIXSORT_THREADS;
    const int RADIX_DIGITS = 1 << RADIX_BITS;

    RakingReduction2<K, V, 0, RADIX_BITS, 0, PreShiftFunctor2> <<<nr_blocks, threads>>>
      (thrust::raw_pointer_cast(&d_spine[0]),
       thrust::raw_pointer_cast(&d_keys[0]),
       thrust::raw_pointer_cast(&d_offsets[0]));
    
    if (debug) {
      thrust::host_vector<int> h_spine(&d_spine[0], &d_spine[nr_blocks * RADIX_DIGITS]);
      print("h_spine", h_spine);
    }

    top_scan();

    if (debug) {
      thrust::host_vector<int> h_spine(&d_spine[0], &d_spine[nr_blocks * RADIX_DIGITS]);
      print("h_spine scan", h_spine);
    }

    ScanScatterDigits2<K, V, 0, RADIX_BITS, 0, PreShiftFunctor2, PostShiftFunctor2> <<<nr_blocks, threads>>>
      (thrust::raw_pointer_cast(&d_spine[0]),
       thrust::raw_pointer_cast(&d_keys[0]),
       thrust::raw_pointer_cast(&d_alt_keys[0]),
       thrust::raw_pointer_cast(&d_values[0]),
       thrust::raw_pointer_cast(&d_alt_values[0]),
       thrust::raw_pointer_cast(&d_offsets[0]));

    if (debug) {
      thrust::host_vector<K> h_keys = d_keys;
      print("h_keys", h_keys);
      thrust::host_vector<K> h_alt_keys = d_alt_keys;
      print("h_alt_keys", h_alt_keys);
      thrust::host_vector<K> h_alt_values = d_alt_values;
      print("h_alt_values", h_alt_values);
    }
  }
};
  
static void
my_sort_3()
{
  static const int RADIX_BITS = 4;

  K keys[] = { 0, 1, 0, 0,  0, 1, 1, 0,
	       1, 0, 1, 1,  0, 1, 1, 1 };
  int N = 16;
  int nr_blocks = 2;
  int offsets[] = { 0, 8, 16 };

  thrust::device_vector<K> d_keys(keys, keys + N);
  thrust::device_vector<int> d_offsets(offsets, offsets + nr_blocks + 1);

  struct my_sort<RADIX_BITS> ms(d_keys, d_offsets);
  ms.sort();
}

static void
my_sort_4()
{
  static const int RADIX_BITS = 4;

  const int N = 16;
  const int nr_blocks = 4;
  K keys[N] = { 0, 0, 0, 1,  1, 1, 0, 2,
		2, 2, 1, 3,  3, 3, 2, 3 };
  int offsets[nr_blocks + 1] = { 0, 4, 8, 12, 16 };

  thrust::device_vector<K> d_keys(keys, keys + N);
  thrust::device_vector<int> d_offsets(offsets, offsets + nr_blocks + 1);

  struct my_sort<RADIX_BITS> ms(d_keys, d_offsets);
  ms.sort();
}

static void
my_sort_5()
{
  static const int RADIX_BITS = 4;

  const int N = 16;
  const int nr_blocks = 4;
  K keys[N] = { 0, 0, 0, 1,  1, 1, 0, 2,
		2, 2, 1, 3,  3, 3, 2, 3 };
  int offsets[nr_blocks + 1] = { 0, 4, 8, 12, 16 };

  thrust::device_vector<K> d_keys(keys, keys + N);
  thrust::device_vector<int> d_offsets(offsets, offsets + nr_blocks + 1);

  struct my_sort_b<RADIX_BITS> ms(d_keys, d_offsets);
  ms.sort();
}

static void
my_sort_6()
{
  static const int RADIX_BITS = 4;

  const int nr_blocks = 4;
  const int ppb = 4;
  const int N = ppb * nr_blocks;

  thrust::device_vector<K> d_keys(N);
  thrust::device_vector<int> d_offsets(nr_blocks + 1);
  for (int b = 0; b <= nr_blocks; b++) {
    d_offsets[b] = ppb * b;
  }
  for (int b = 0; b < nr_blocks; b++) {
    for (int n = d_offsets[b]; n < d_offsets[b+1]; n++) {
      if (n - d_offsets[b] == 2 && b > 0) {
	d_keys[n] = b - 1;
      } else if (n - d_offsets[b] == 3 && b < nr_blocks - 1) {
	d_keys[n] = b + 1;
      } else {
	d_keys[n] = b;
      }
    }
  }

  struct my_sort_b<RADIX_BITS> ms(d_keys, d_offsets);
  ms.sort();
}

static void
my_sort_7()
{
  const int b_mx[3] = { NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z };
  const int nr_blocks = NBLOCKS_Y * NBLOCKS_Z;
  const int ppb = 512;

  thrust::host_vector<int> h_offsets(nr_blocks + 1);
  for (int b = 0; b <= nr_blocks; b++) {
    h_offsets[b] = ppb * b + b;
  }
  int N = h_offsets[nr_blocks];
  thrust::host_vector<K> h_keys(N);
  for (int b = 0; b < nr_blocks; b++) {
    for (int n = h_offsets[b]; n < h_offsets[b+1]; n++) {
      int block_pos[3];
      block_idx_to_block_pos(b, block_pos);
#if 0
      if (n - h_offsets[b] == 2) {
	block_pos[1]--;
      }
      if (n - h_offsets[b] == 3) {
	block_pos[1]++;
      }
      if (n - h_offsets[b] == 6) {
	block_pos[2]--;
      }
      if (n - h_offsets[b] == 7) {
	block_pos[2]++;
      }
#endif
      block_pos[1] &= NBLOCKS_Y - 1;
      block_pos[2] &= NBLOCKS_Z - 1;
      h_keys[n] = block_pos_to_block_idx(block_pos);
    }
  }

  thrust::device_vector<K> d_keys = h_keys;
  thrust::device_vector<int> d_offsets = h_offsets;

#if 1
  thrust::device_vector<V> d_values(N);
  thrust::device_vector<K> d_alt_keys(N);
  thrust::device_vector<V> d_alt_values(N);
  thrust::sequence(d_values.begin(), d_values.end());

  sort_pairs_device_2(raw_pointer_cast(&d_keys[0]), raw_pointer_cast(&d_values[0]),
		      raw_pointer_cast(&d_alt_keys[0]), raw_pointer_cast(&d_alt_values[0]),
		      N,
		      raw_pointer_cast(&d_offsets[0]), b_mx);

  thrust::host_vector<K> h_keys_res = d_keys;
  thrust::host_vector<V> h_values = d_values;
  int last = h_keys_res[0];
  for (int i = 1; i < h_keys_res.size(); i++) {
    if (h_keys_res[i] < last) {
      printf("order! i %d last %d this %d\n", i, last, h_keys_res[i]);
    }
    last = h_keys_res[i];
    assert(h_keys_res[i] == h_keys[h_values[i]]);
  }
#else
  static const int RADIX_BITS = 4;

  struct my_sort_b<RADIX_BITS> ms(d_keys, d_offsets);
  ms.sort();

  thrust::host_vector<unsigned int> h_alt_keys = ms.d_alt_keys;
  thrust::host_vector<unsigned int> h_alt_values = ms.d_alt_values;
  int last = h_alt_keys[0];
  for (int n = 1; n < h_alt_keys.size(); n++) {
    if (h_alt_keys[n] < last) {
      printf("order! n %d last %d this %d\n", n, last, h_alt_keys[n]);
    }
    last = h_alt_keys[n];
    assert(h_alt_keys[n] == h_keys[h_alt_values[n]]);
  }
#endif
}

#define BLOCKSIZE_X 1
#define NO_CHECKERBOARD

struct sort_info {
  int bdims[3]; // # blocks per direction
};

static inline void
block_idx_to_block_pos(int block_idx, int block_pos[3], struct sort_info *si)
{
#if BLOCKSIZE_X == 1
  block_pos[2] = block_idx / si->bdims[1];
  block_pos[1] = block_idx - si->bdims[1] * block_pos[2];
#else
#error TBD
  //  cell_map_1to3(map, block_idx, block_pos);
#endif
}

static inline int
block_pos_to_block_idx(int block_pos[3], struct sort_info *si)
{
#if BLOCKSIZE_X == 1
  return block_pos[2] * si->bdims[1] + block_pos[1];
#else
#error TBD
  //  return cell_map_3to1(map, block_pos);
#endif
}

static inline int
get_dir1(int block_idx, int block_pos_0[3], struct sort_info *si)
{
  int dir[3];
  block_idx_to_block_pos(block_idx, dir, si);
  for (int d = 0; d < 3; d++) {
    dir[d] -= block_pos_0[d];
    if (dir[d] == 1 - si->bdims[d]) 
      dir[d] = 1;
    if (dir[d] == -1 + si->bdims[d]) 
	dir[d] = -1;
  }
  int dir1 = (1 - dir[2]) * 3 + (1 - dir[1]);
  assert(dir1 >= 0 && dir1 < 9);
  return dir1;
}

// ----------------------------------------------------------------------
// my_sort

void
my_sort_1(int n, int *keys, int *ids, int nr_blocks)
{
  int *cnts = (int *) calloc(nr_blocks, sizeof(*cnts));
  // count
  for (int i = 0; i < n; i++) {
    cnts[keys[i]]++;
  }
  // pfx sum
  int sum = 0;
  for (int b = 0; b < nr_blocks; b++) {
    int cnt = cnts[b];
    cnts[b] = sum;
    sum += cnt;
  }
  // scatter
  for (int i = 0; i < n; i++) {
    int key = keys[i];
    ids[cnts[key]++] = i;
  }
  free(cnts);
}

static void
my_sort_2_block_count(int start, int end, int *keys, int *cnts2, int b,
		      struct sort_info *si)
{
  int b_block_pos[3];
  block_idx_to_block_pos(b, b_block_pos, si);

  for (int i = start; i < end; i++) {
    int dir1 = get_dir1(keys[i], b_block_pos, si);
    cnts2[dir1]++;
  }
}

static void
my_sort_2_block_scatter(int start, int end, int *keys, int *ids, int *cnts2,
			int b, struct sort_info *si)
{
  int b_block_pos[3];
  block_idx_to_block_pos(b, b_block_pos, si);

  for (int i = start; i < end; i++) {
    int dir1 = get_dir1(keys[i], b_block_pos, si);
    ids[cnts2[dir1]++] = i;
  }
}

void
my_sort_2(int n, int *keys, int *ids, int nr_blocks, int *offsets, struct sort_info *si)
{
  int *b2_cnts = (int *) calloc(9 * nr_blocks, sizeof(*b2_cnts));
  // count per block
  for (int b = 0; b < nr_blocks; b++) {
    my_sort_2_block_count(offsets[b], offsets[b+1], keys, b2_cnts + 9*b, b,
			  si);
  }

  // all block prefix sum
  // the different order in dir makes the sort not stable -- this could be fixed
  // if the block indices weren't doing the bit games
  int sum = 0;
  for (int b = 0; b < nr_blocks; b++) {
    int b_block_pos[3];
    block_idx_to_block_pos(b, b_block_pos, si);
    int dir[3] = {};
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	int dir1 = (1 + dir[2]) * 3 + (1 + dir[1]);
	int c_block_pos[3];
	for (int d = 0; d < 3; d++) {
	  c_block_pos[d] = b_block_pos[d] + dir[d];
	  if (c_block_pos[d] < 0)
	    c_block_pos[d] +=  si->bdims[d];
	  if (c_block_pos[d] >= si->bdims[d])
	    c_block_pos[d] -= si->bdims[d];
	}
	int c = block_pos_to_block_idx(c_block_pos, si);
	int cnt = b2_cnts[9*c + dir1];
	b2_cnts[9*c + dir1] = sum;
	sum += cnt;
      }
    }
  }

  // scatter
  for (int b = 0; b < nr_blocks; b++) {
    my_sort_2_block_scatter(offsets[b], offsets[b+1], keys, ids, b2_cnts + 9*b,
			    b, si);
  }
}

void
my_sort_check(particles_cuda_t *pp, int *keys, int *ids_ref)
{
  if (1) {
    my_sort_7();
    exit(0);
  }

#if 0
  int *offsets = (int *) malloc((pp->nr_blocks + 1) * sizeof(*offsets));
  cuda_copy_offsets_from_dev(pp, offsets);

  int *ids = (int *) malloc(pp->n_part * sizeof(*ids));
  //  my_sort_1(pp->n_part, keys, ids, pp->nr_blocks);
  struct sort_info si = {
    .bdims = { pp->b_mx[0], pp->b_mx[1], pp->b_mx[2] },
  };
  my_sort_2(pp->n_part, keys, ids, pp->nr_blocks, offsets, &si);

  free(offsets);

  // check result
  int last = -1;
  for (int i = 0; i < pp->n_part; i++) {
#if 0 // only if stable
    if (ids[i] != ids_ref[i]) {
      printf("%d: %d ref %d\n", i, ids[i], ids_ref[i]);
    }
#endif
    if (keys[ids[i]] < last) {
      printf("i %d last %d now %d\n", i, keys[ids[i]], last);
    }
    last = keys[ids[i]];
  }
  free(ids);
#endif
}

void
print(const char *title, const thrust::host_vector<unsigned int> &h_vec)
{
  printf("%s:", title);
  for (int i = 0; i < h_vec.size(); i++) {
    printf(" %2u", h_vec[i]);
  }
  printf("\n");
}

void
print_matrix(const thrust::host_vector<int> &h_matrix, int nr_blocks)
{
  printf("\n");
  for (int iy = 0; iy < nr_blocks; iy++) {
    for (int ix = 0; ix < nr_blocks; ix++) {
      int val = h_matrix[iy * nr_blocks + ix];
#if 1
      if (val) {
	printf("%2d ", val);
      } else {
	printf(" . ");
      }
#else
      if (val) {
	printf("*", val);
      } else {
	printf(".");
      }
#endif
    }
    printf("\n");
  }
  printf("\n");
}

#if 0
EXTERN_C void
sort_pairs_device_2(unsigned int *d_keys, unsigned int *d_values, int n,
		      int *d_offsets, const int b_mx[3])
{
  const int threads = B40C_RADIXSORT_THREADS;
  const int RADIX_BITS = 4;
  const int RADIX_DIGITS = 1 << RADIX_BITS;
  const bool debug = false;
  const bool check = false;

  static int pr_A, pr_B, pr_C;
  if (!pr_A) {
    pr_A = prof_register("sort_bottom_sum", 1., 0, 0);
    pr_B = prof_register("sort_top_scan", 1., 0, 0);
    pr_C = prof_register("sort_bottom_scan", 1., 0, 0);
  }

  assert(b_mx[0] == NBLOCKS_X);
  assert(b_mx[1] == NBLOCKS_Y);
  assert(b_mx[2] == NBLOCKS_Z);
  int nr_blocks = b_mx[0] * b_mx[1] * b_mx[2];

  thrust::device_vector<K> d_alt_keys(n, -1);
  thrust::device_vector<V> d_alt_values(n);

  int spine_cycle_elements = B40C_RADIXSORT_SPINE_CYCLE_ELEMENTS;
  int spine_cycles = (nr_blocks * (1 << RADIX_BITS) + spine_cycle_elements - 1) /
    spine_cycle_elements;
  int spine_elements = spine_cycles * B40C_RADIXSORT_SPINE_CYCLE_ELEMENTS;

  thrust::device_vector<int> d_spine(spine_elements);

#if 1
  prof_start(pr_A);
  RakingReduction2<K, V, 0, RADIX_BITS, 0, PreShiftFunctor2> <<<nr_blocks, threads>>>
    (thrust::raw_pointer_cast(&d_spine[0]), d_keys, d_offsets);
  prof_stop(pr_A);
#else
  thrust::host_vector<K> h_keys_tr(n);
  {
     thrust::host_vector<K> h_keys(n);
     thrust::device_ptr<K> d_keys_ptr(d_keys);
     copy(d_keys_ptr, d_keys_ptr + n, h_keys.begin());

     thrust::host_vector<int> h_offsets(nr_blocks + 1);
     thrust::device_ptr<int> d_offsets_ptr(d_offsets);
     copy(d_offsets_ptr, d_offsets_ptr + nr_blocks + 1, h_offsets.begin());

     for (int b = 0; b < nr_blocks; b++) {
       for (int i = h_offsets[b]; i < h_offsets[b+1]; i++) {
	 h_keys_tr[i] = block_idx_to_dir1(b, h_keys[i]);
       }
     }
  }
  thrust::device_vector<K> d_keys_tr = h_keys_tr;
  
  RakingReduction2<K, V, 0, RADIX_BITS, 0, NopFunctor<K> > <<<nr_blocks, threads>>>
    (thrust::raw_pointer_cast(&d_spine[0]), raw_pointer_cast(&d_keys_tr[0]), d_offsets);
#endif
  
  if (0&&debug) {
    thrust::host_vector<int> h_spine(&d_spine[0], &d_spine[nr_blocks * RADIX_DIGITS]);
    print("h_spine", h_spine);
  }

  thrust::host_vector<int> h_matrix(nr_blocks * nr_blocks);
  prof_start(pr_B);
  {
    thrust::host_vector<int> h_spine = d_spine;

    for (int b = 0; b < nr_blocks; b++) {
      for (int d = 0; d < RADIX_DIGITS; d++) {
	int val = h_spine[d * nr_blocks + b];
	if (!val)
	  continue;

	int b2 = dir1_to_block_idx(b, d);
	h_matrix[b2 * nr_blocks + b] = val;
      }
    }

    if (debug) {
      print_matrix(h_matrix, nr_blocks);
    }

    if (check) {
      thrust::host_vector<K> h_keys(n);
      thrust::device_ptr<K> d_keys_ptr(d_keys);
      copy(d_keys_ptr, d_keys_ptr + n, h_keys.begin());

      thrust::host_vector<int> h_offsets(nr_blocks + 1);
      thrust::device_ptr<int> d_offsets_ptr(d_offsets);
      copy(d_offsets_ptr, d_offsets_ptr + nr_blocks + 1, h_offsets.begin());

      for (int b = 0; b < nr_blocks; b++) {
	int block_offset = h_offsets[b];
	int block_elements = h_offsets[b+1] - h_offsets[b];
	int cycles = block_elements / 512;
	printf("b %d cyc %d left %d\n", b, cycles, block_elements - cycles * 512);
      }

      thrust::host_vector<int> h_matrix2(nr_blocks * nr_blocks);
      for (int b = 0; b < nr_blocks; b++) {
	for (int i = h_offsets[b]; i < h_offsets[b+1]; i++) {
	  int b2 = h_keys[i];
	  h_matrix2[b2 * nr_blocks + b]++;
	}
      }
      for (int b = 0; b < nr_blocks; b++) {
	for (int b2 = 0; b2 < nr_blocks; b2++) {
	  assert(h_matrix[b2 * nr_blocks + b] == h_matrix2[b2 * nr_blocks + b]);
	}
      }
    }

    int sum = 0;
    for (int b2 = 0; b2 < nr_blocks; b2++) {
      for (int b = 0; b < nr_blocks; b++) {
	int val = h_matrix[b2 * nr_blocks + b];
	h_matrix[b2 * nr_blocks + b] = sum;
	sum += val;
      }
    }

    if (0&&debug) {
      print_matrix(h_matrix, nr_blocks);
    }

    // OPT can skip the ones where h_spine == 0 (?)
    for (int b = 0; b < nr_blocks; b++) {
      for (int b2 = 0; b2 < nr_blocks; b2++) {
	int d = block_idx_to_dir1_oob(b, b2);
	if (d >= 0) {
	  int val = h_matrix[b2 * nr_blocks + b];
	  h_spine[d * nr_blocks + b] = val;
	}
      }
    }
    d_spine = h_spine;
  }
  prof_stop(pr_B);

  if (0&&debug) {
    thrust::host_vector<int> h_spine(&d_spine[0], &d_spine[nr_blocks * RADIX_DIGITS]);
    print("h_spine scan", h_spine);
  }

  prof_start(pr_C);
#if 1
  ScanScatterDigits2<K, V, 0, RADIX_BITS, 0, PreShiftFunctor2, PostShiftFunctor2> <<<nr_blocks, threads>>>
    (thrust::raw_pointer_cast(&d_spine[0]),
     d_keys,
     thrust::raw_pointer_cast(&d_alt_keys[0]),
     d_values,
     thrust::raw_pointer_cast(&d_alt_values[0]),
     d_offsets);
#else
  ScanScatterDigits2<K, V, 0, RADIX_BITS, 0, NopFunctor<K>, PostShiftFunctor2> <<<nr_blocks, threads>>>
    (thrust::raw_pointer_cast(&d_spine[0]),
     thrust::raw_pointer_cast(&d_keys_tr[0]),
     thrust::raw_pointer_cast(&d_alt_keys[0]),
     d_values,
     thrust::raw_pointer_cast(&d_alt_values[0]),
     d_offsets);
#endif
  prof_stop(pr_C);

  if (check) {
    thrust::host_vector<K> h_keys(n);
    thrust::device_ptr<K> d_keys_ptr(d_keys);
    copy(d_keys_ptr, d_keys_ptr + n, h_keys.begin());
    
    thrust::host_vector<int> h_offsets(nr_blocks + 1);
    thrust::device_ptr<int> d_offsets_ptr(d_offsets);
    copy(d_offsets_ptr, d_offsets_ptr + nr_blocks + 1, h_offsets.begin());
    
    thrust::host_vector<K> h_alt_keys2(n, -1);
    thrust::host_vector<K> h_alt_values2(n, -1);
    for (int b = 0; b < nr_blocks; b++) {
      printf("h_matrix[b,b] = %d\n", h_matrix[b * nr_blocks + b]);
      for (int i = h_offsets[b]; i < h_offsets[b+1]; i++) {
	int b2 = h_keys[i];
	int dir1 = block_idx_to_dir1(b, b2);
	int b2_ = dir1_to_block_idx(b, dir1);
	assert(b2 == b2_);

	int pos = h_matrix[b2 * nr_blocks + b]++;
	if (0&&pos == 6399) {
	  printf("6399: %d -> %d b %d b2 %d\n", h_alt_keys2[pos], b2, b, b2);
	}
	if (0&&pos == 9599) {
	  printf("9599: %d -> %d b %d b2 %d\n", h_alt_keys2[pos], b2, b, b2);
	}
	h_alt_keys2[pos] = b2;
	h_alt_values2[pos] = i;
      }
    }

    thrust::host_vector<K> h_alt_keys = d_alt_keys;
    thrust::host_vector<K> h_alt_values = d_alt_values;
    bool bug = false;
    for (int i = 0; i < n; i++) {
      if (h_alt_keys[i] != h_alt_keys2[i]) {
	printf("i %d: res %d ref %d\n", i, h_alt_keys[i], h_alt_keys2[i]);
	bug = true;
      }
      if (h_alt_values[i] != h_alt_values2[i]) {
	printf("i %d: val %d ref %d\n", i, h_alt_values[i], h_alt_values2[i]);
	bug = true;
      }
    }
    assert(!bug);
  }

  if (0&&debug) {
    thrust::host_vector<K> h_alt_keys = d_alt_keys;
    print("h_alt_keys", h_alt_keys);
    thrust::host_vector<K> h_alt_values = d_alt_values;
    print("h_alt_values", h_alt_values);
  }

  if (check) {
    thrust::host_vector<K> h_alt_keys = d_alt_keys;
    int last = h_alt_keys[0];
    for (int i = 0; i < n; i++) {
      int key = h_alt_keys[i];
      assert(key >= 0 && key < nr_blocks);
      if (key < last) {
	printf("i %d, last %d now %d\n", i, last, key);
      }
      assert(key >= last);
      last = key;
    }
  }

  thrust::device_ptr<K> d_keys_ptr(d_keys);
  thrust::copy(d_alt_keys.begin(), d_alt_keys.end(), d_keys_ptr);
  thrust::device_ptr<K> d_values_ptr(d_values);
  thrust::copy(d_alt_values.begin(), d_alt_values.end(), d_values_ptr);
}
#endif
