
#ifndef CUDA_IFACE_BND_H
#define CUDA_IFACE_BND_H

// FIXME, better call it cuda_mfields_real_t
typedef float fields_cuda_real_t;

// ----------------------------------------------------------------------
// cuda_mfields_bnd_patch

struct cuda_mfields_bnd_patch {
  fields_cuda_real_t *arr_off;
  int im[3];
  int ib[3];
  fields_cuda_real_t *arr;
};

// ----------------------------------------------------------------------
// cuda_mfields_bnd

struct cuda_mfields_bnd;

// FIXME This is the alternative to doing it with json, I guess...
// Should settle on one or the other eventually
struct cuda_mfields_bnd_entry {
  int patch;
  int nei_patch;
  int dir1;
};

struct cuda_mfields_bnd_params {
  int n_patches;
  int im[3];
  int ib[3];

  int n_recv_entries; // how many entries we have in the following neighbor list
  struct cuda_mfields_bnd_entry *recv_entry;
};


struct cuda_mfields_bnd *cuda_mfields_bnd_create(void);
void cuda_mfields_bnd_destroy(struct cuda_mfields_bnd *cbnd);
void cuda_mfields_bnd_ctor(struct cuda_mfields_bnd *cbnd, struct cuda_mfields_bnd_params *prm);
void cuda_mfields_bnd_dtor(struct cuda_mfields_bnd *cbnd);
struct cuda_mfields_bnd_patch *cuda_mfields_bnd_get_patch(struct cuda_mfields_bnd *cbnd, int p);


void cuda_mfields_bnd_from_device_inside(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
					 int mb, int me);
void cuda_mfields_bnd_from_device_inside_only(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
					      int mb, int me);
void cuda_mfields_bnd_to_device_outside(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
					int mb, int me);
void cuda_mfields_bnd_to_device_inside(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
					int mb, int me);
void cuda_mfields_bnd_fill_ghosts_local(struct cuda_mfields_bnd *cbnd, struct cuda_mfields *cmflds,
					int mb, int me);


// special cases for single GPU, single patch, which can be handled
// entirely on the GPU

void cuda_fill_ghosts_periodic_yz(struct cuda_mfields *cmflds, int p, int mb, int me);
void cuda_fill_ghosts_periodic_z(struct cuda_mfields *cmflds, int p, int mb, int me);
void cuda_add_ghosts_periodic_yz(struct cuda_mfields *cmflds, int p, int mb, int me);
void cuda_add_ghosts_periodic_z(struct cuda_mfields *cmflds, int p, int mb, int me);

// These fields are for field boundary exchange, so they could, eventually, miss the
// interior part

#undef F3_CF_BOUNDS_CHECK

#ifdef F3_CF_BOUNDS_CHECK
#define F3_CF_OFF(cf, fldnr, jx,jy,jz) ({				\
  assert(jx == 0); /* FIXME yz only! */				        \
  assert(jx >= (cf)->ib[0] && jx < (cf)->ib[0] + (cf)->im[0]);		\
  assert(jy >= (cf)->ib[1] && jy < (cf)->ib[1] + (cf)->im[1]);		\
  assert(jz >= (cf)->ib[2] && jz < (cf)->ib[2] + (cf)->im[2]);		\
  int __off = (((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx); \
  __off; })

#else
#define F3_CF_OFF(cf, fldnr, jx,jy,jz)					\
  ((((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx))
#endif

#define F3_CF(cf, fldnr, jx,jy,jz)					\
  ((cf)->arr_off[F3_CF_OFF(cf, fldnr, jx,jy,jz)])


// ----------------------------------------------------------------------
// routines for actual domain boundaries

void cuda_conducting_wall_H_lo_hi_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_E_lo_hi_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_J_lo_hi_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_H_lo_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_H_hi_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_E_lo_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_E_hi_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_J_lo_y(struct cuda_mfields *cmflds, int p);
void cuda_conducting_wall_J_hi_y(struct cuda_mfields *cmflds, int p);


#endif
