#include "psc.h"

// ======================================================================

int
cell_map_init(struct cell_map *map, const int dims[3], const int blocksize[3])
{
  map->b_bits_max = 0;
  for (int d = 0; d < 3; d++) {
    assert(dims[d] % blocksize[d] == 0);
    map->dims[d] = dims[d] / blocksize[d];
    int bits = 0;
    while (blocksize[d] > (1 << bits))
      bits++;
    assert(blocksize[d] == (1 << bits));
    map->b_bits[d] = bits;
    if (bits > map->b_bits_max) {
      map->b_bits_max = bits;
    }
  }
  map->N = dims[0] * dims[1] * dims[2];
  return map->N;
}

void
cell_map_free(struct cell_map *map)
{
}

int
cell_map_3to1(struct cell_map *map, int i[3])
{
  // FIXME, don't change i[]
  for (int d = 0; d < 3; d++) {
    assert(i[d] < (map->dims[d] << map->b_bits[d]));
  }
  int cnt = 0;
  unsigned int idx = 0;
  for (int b = 0; b < map->b_bits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b < map->b_bits[d]) {
	if (i[d] & 1) {
	  idx |= (1 << cnt);
	}
	i[d] >>= 1;
	cnt++;
      }
    }
  }
  idx |= (((i[2]) * (map->dims[1]) + i[1]) * map->dims[0] + i[0]) << cnt;
#if 1
  if (idx >= map->N) {
    printf("idx %d N %d %d %d %d dims %d %d %d\n", idx, map->N, i[0], i[1], i[2],
	   map->dims[0], map->dims[1], map->dims[2]);
  }
#endif
  assert(idx < map->N);
  return idx;
}

void
cell_map_1to3(struct cell_map *map, int idx, int i[3])
{
  for (int d = 0; d < 3; d++) {
    i[d] = 0;
  }
  for (int b = 0; b < map->b_bits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b < map->b_bits[d]) {
	if (idx & 1) {
	  i[d] |= (1 << b);
	}
	idx >>= 1;
      }
    }
  }

  i[0] |= (idx % map->dims[0]) << map->b_bits[0];
  idx /= map->dims[0];
  i[1] |= (idx % map->dims[1]) << map->b_bits[1];
  idx /= map->dims[1];
  i[2] |= idx << map->b_bits[2];
  for (int d = 0; d < 3; d++) {
    assert(i[d] < (map->dims[d] << map->b_bits[d]));
  }
}


