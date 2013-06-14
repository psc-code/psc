#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <bitfield3d.h>

void bitfield3d_create(struct bitfield3d* field, unsigned int size[3])
{
  unsigned int length = ((size[0]*size[1]*size[2])+(sizeof(bitfield3d_type) * 8 - 1))/ (sizeof(bitfield3d_type) * 8);
  field->ystride = size[0];
  field->zstride = size[0] * size[1];

  field->values = malloc(length*sizeof(bitfield3d_type));
  memset(field->values, 0, length*sizeof(bitfield3d_type));
  field->length = length;
}

int bitfield3d_isset(const struct bitfield3d* field, int idx[3])
{
  unsigned int i = idx[0] + idx[1] * field->ystride + idx[2] * field->zstride;
  unsigned int pi = i % (sizeof(unsigned int) * 8);
  i /= (sizeof(unsigned int) * 8);
  bitfield3d_type v = field->values[i];
  bitfield3d_type mask = 1 << pi;
  
  return ((v & mask) != 0);
}

void bitfield3d_set(struct bitfield3d* field, int idx[3], int value)
{
  unsigned int i = idx[0] + idx[1] * field->ystride + idx[2] * field->zstride;
  unsigned int pi = i % (sizeof(unsigned int) * 8);
  i /= (sizeof(unsigned int) * 8);
  assert( i <= field->length && i >= 0);  //You have tried to set a value out of bounds
  
  bitfield3d_type* v = &field->values[i];  
  bitfield3d_type mask = 1 << pi;
  if(value)
  {
    *v |= mask;
  }
  else
  {
    *v &= ~mask;
  }
}

void bitfield3d_fill(struct bitfield3d* field, int value)
{
  bitfield3d_type v = 0;
  if(value != 0)
  {
    v = ~0;
  }
  
  for(int i=0; i < field->length; ++i) field->values[i] = v;
}

int bitfield3d_compare(const struct bitfield3d* field1, const struct bitfield3d* field2)
{
  for(int i=0; i < field1->length; ++i)
  {
    if(field1->values[i] != field2->values[i]) return 1;
  }
  return 0;
}

void bitfield3d_merge(struct bitfield3d* out, struct bitfield3d* field2)
{
  //(Incomplete) test if the dimensions match
  assert(out->length == field2->length);
  
  for(int i=0; i<out->length; ++i)
  {
    out->values[i] |= field2->values[i];
  }
}

void bitfield3d_destroy(struct bitfield3d* field)
{
  free(field->values);
}

int bitfield3d_count_bits_set(const struct bitfield3d* field)
{
// FIXME! This is a hack, builtin_popcount doesn't necessarily
// have to do with SSE2, but rather should be tested for by configure,
// and an alternate implementation should be offered when needed.
#ifdef USE_SSE2
  int bits = 0;
  for(int i=0; i<field->length; ++i)
  {
#ifdef BITFIELD64
    bits += __builtin_popcountll (field->values[i]);	//Will compile to popcnt if available
#else
    bits += __builtin_popcount(field->values[i]);
#endif
  }
  return bits;
#else
  assert(0);
#endif
}

void bitfield3d_copy(struct bitfield3d* out, struct bitfield3d* src)
{
  out->ystride = src->ystride;
  out->zstride = src->zstride;
  out->length = src->length;
  out->values = malloc(src->length * sizeof(bitfield3d_type));
  memcpy(out->values, src->values, sizeof(bitfield3d_type) * src->length);
}
