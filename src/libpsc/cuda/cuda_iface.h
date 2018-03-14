
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include "mrc_json.h"
#include "psc_fields_single.h"

#include <particles_simple.hxx> // FIXME?
#include <grid.hxx>

#if 1
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

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

struct particle_cuda_t : psc_particle<float> {};

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

struct MparticlesCuda : MparticlesBase
{
  using Self = MparticlesCuda;
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;
  using particle_buf_t = psc_particle_cuda_buf_t;
  
  MparticlesCuda(const Grid_t& grid);

  MparticlesCuda(const MparticlesCuda&) = delete;
  ~MparticlesCuda();

  int get_n_prts() const override;
  void get_size_all(uint *n_prts_by_patch) const override;
  void reserve_all(const uint *n_prts_by_patch) override;
  void resize_all(const uint *n_prts_by_patch) override;
  void to_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  void from_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  void setup_internals();
  void inject_buf(cuda_mparticles_prt *buf, uint *buf_n_by_patch);

  template<typename MP>
  static void copy_to(MparticlesCuda& mprts_cuda, MP& mprts);

  template<typename MP>
  static void copy_from(MparticlesCuda& mprts_cuda, MP& mprts);

  template<typename MP>
  static void copy_to(MparticlesBase& mprts_cuda, MparticlesBase& mprts);

  template<typename MP>
  static void copy_from(MparticlesBase& mprts_cuda, MparticlesBase& mprts);

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  const int *patch_get_b_mx(int p);
  
  cuda_mparticles* cmprts() { return cmprts_; }

  struct patch_t
  {
    using buf_t = particle_buf_t;
    
    patch_t(MparticlesCuda& mp, int p)
      : mp_(mp), p_(p), pi_(mp.grid())
    {}

    particle_buf_t& get_buf() { assert(0); static particle_buf_t fake{}; return fake; } // FIXME

    int blockPosition(real_t xi, int d) const { return pi_.blockPosition(xi, d); }
    Int3 blockPosition(const Real3& xi) const { return pi_.blockPosition(xi); }
    int validCellIndex(const particle_t& prt) const { return pi_.validCellIndex(&prt.xi); }
  
    const int* get_b_mx() const;

  private:
    MparticlesCuda& mp_;
    int p_;
    ParticleIndexer<real_t> pi_;
  };

  const patch_t& operator[](int p) const { return patches_[p]; }
  patch_t&       operator[](int p)       { return patches_[p]; }

private:
  cuda_mparticles* cmprts_;
  std::vector<patch_t> patches_;

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
