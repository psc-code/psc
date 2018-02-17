
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include "mrc_json.h"
#include "psc_fields_single.h"

#include <particles.hxx>
#include <grid.hxx>

struct BS
{
  using x = std::integral_constant<unsigned int, 1>;
  using y = std::integral_constant<unsigned int, 4>;
  using z = std::integral_constant<unsigned int, 4>;
};

// ----------------------------------------------------------------------
// float_3 etc

typedef float float_3[3];
typedef double double_3[3];
#ifndef __CUDACC__
typedef struct { float x; float y; float z; float w; } float4;
#endif

// ----------------------------------------------------------------------
// cuda_base

void cuda_base_init(void);

// ----------------------------------------------------------------------
// psc_mparticles_cuda

using particle_cuda_real_t = float;

struct particle_cuda_t : psc_particle<particle_cuda_real_t> {};

using psc_particle_cuda_buf_t = std::vector<particle_cuda_t>;

// ----------------------------------------------------------------------
// cuda_mparticles_prt

struct cuda_mparticles_prt {
  float xi[3];
  float pxi[3];
  int kind;
  float qni_wni;
};

struct cuda_mparticles;

struct psc_mparticles_cuda : psc_mparticles_base
{
  using particle_t = particle_cuda_t;
  
  psc_mparticles_cuda(Grid_t& grid);
  psc_mparticles_cuda(const psc_mparticles_cuda&) = delete;
  ~psc_mparticles_cuda();

  int get_n_prts() const override;
  void get_size_all(uint *n_prts_by_patch) const override;
  void reserve_all(const uint *n_prts_by_patch) override;
  void resize_all(const uint *n_prts_by_patch) override;
  void to_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  void from_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  void setup_internals();
  void inject_buf(cuda_mparticles_prt *buf, uint *buf_n_by_patch);

  static void copy_from_single(struct psc_mparticles *mprts_cuda,
			       struct psc_mparticles *mprts, uint flags);
  static void copy_to_single(struct psc_mparticles *mprts_cuda,
			     struct psc_mparticles *mprts, uint flags);
  static void copy_from_double(struct psc_mparticles *mprts_cuda,
			       struct psc_mparticles *mprts, uint flags);
  static void copy_to_double(struct psc_mparticles *mprts_cuda,
			     struct psc_mparticles *mprts, uint flags);
  
  const int *patch_get_b_mx(int p);
  
  cuda_mparticles* cmprts() { return cmprts_; }

private:
  cuda_mparticles* cmprts_;

  template<typename MP>
  friend struct bnd_particles_policy_cuda;
};

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields;

void cuda_mfields_calc_dive_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf, int p);

void cuda_push_fields_E_yz(struct cuda_mfields *cmflds, float dt);
void cuda_push_fields_H_yz(struct cuda_mfields *cmflds, float dt);
void cuda_marder_correct_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			    int p, float fac[3],
			    int ly[3], int ry[3], int lz[3], int rz[3]);

// ----------------------------------------------------------------------
// cuda_push_mprts

void cuda_push_mprts_yz(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds,
			const int bs[3], bool ip_ec, bool deposit_vb_3d, bool currmem_global);

void cuda_push_mprts_xyz(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds);

// ----------------------------------------------------------------------
// cuda_moments

void cuda_moments_yz_rho_1st_nc(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres);
void cuda_moments_yz_n_1st(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres);

// ----------------------------------------------------------------------
// cuda_heating_run_foil

struct cuda_heating_foil {
  // params
  float zl;
  float zh;
  float xc;
  float yc;
  float rH;
  float T;
  float Mi;
  int kind;

  // state (FIXME, shouldn't be part of the interface)
  float fac;
  float heating_dt;
};

void cuda_heating_setup_foil(struct cuda_heating_foil *foil);
void cuda_heating_run_foil(struct cuda_mparticles *cmprts);

// FIXME, mv elsewhere
#define HERE printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
