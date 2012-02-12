///@file bitfield3d.h \brief Implements a 3D bitfield (collection of booleans) at 32-bit/64-bit granularity
#ifndef _BITFIELD3D
#define _BITFIELD3D

#if __x86_64__ || __SIZEOF__SIZE_T__ == 8
#define BITFIELD64
typedef unsigned long long bitfield3d_type;
#else
typedef unsigned long bitfield3d_type;		//Must be unsigned
#endif

struct bitfield3d
{
  bitfield3d_type *values;
  unsigned int ystride, zstride;
  unsigned int length;	//Number of (unsigned int)s allocated
};

///\brief Creates a bitfield
///@param size The dimensions of the field
///\note Please destroy all created bitfields manually via bitfield3d_destroy
void bitfield3d_create(struct bitfield3d* field, unsigned int size[3]);

///\brief Creates a deep copy of the bitfield \a src
///Use bitfield3d_destroy(out) to cleanup the newly created bitfield
///\note Does NOT cleanup out. You'll have to do this yourself before calling bitfield3d_copy
void bitfield3d_copy(struct bitfield3d* out, struct bitfield3d* src);

///\brief Checks if the specified bit is set
///@return the value of the bit at \a idx
int bitfield3d_isset(const struct bitfield3d* field, int idx[3]);

///\brief Sets the bit at \a idx to \a value
///@param idx the index to set at
///@param value the value to set
void bitfield3d_set(struct bitfield3d* field, int idx[3], int value);

///Cleans up all associated memory of the bitfield
void bitfield3d_destroy(struct bitfield3d* field);

///Merges the two bitfields in-place into \a out
///Acts like a binary OR
///\note Both fields must have the same dimensions
void bitfield3d_merge(struct bitfield3d* out, struct bitfield3d* field2);

///Fills the field with value
///\param value if this is set to 0, the field will be filled with zeroes, otherwise it will be filled with 1s
void bitfield3d_fill(struct bitfield3d* field, int value);

///Compares 2 bitfields binarily
///\return nonzero if the fields are equal, 0 otherwise
int bitfield3d_compare(const struct bitfield3d* field1, const struct bitfield3d* field2);

///Returns the number of bits set in the whole field
///\note While this does not need to iterate all bits, it needs to iterate the whole array of 32/64bit integers, so use cautiously nontheless
int bitfield3d_count_bits_set(const struct bitfield3d* field);

#endif

