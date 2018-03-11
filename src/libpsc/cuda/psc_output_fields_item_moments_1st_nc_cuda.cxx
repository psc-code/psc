
#include "psc_output_fields_item_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "fields_item.hxx"

// ======================================================================
// Moment_rho_1st_nc_cuda

struct Moment_rho_1st_nc_cuda
{
  using mfields_t = PscMfieldsCuda;
  using mparticles_t = PscMparticlesCuda;
  
  constexpr static const char* name = "rho_1st_nc";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "rho_nc_cuda" }; } // FIXME
  constexpr static int flags = 0;

  static void run(mfields_t mres, mparticles_t mprts)
  {
    cuda_mparticles *cmprts = mprts->cmprts();
    cuda_mfields *cmres = mres->cmflds;
    
    mres->zero();
    cuda_moments_yz_rho_1st_nc(cmprts, cmres);
  }
};

// ======================================================================
// n_1st_cuda

struct Moment_n_1st_cuda
{
  using mfields_t = PscMfieldsCuda;
  using mparticles_t = PscMparticlesCuda;

  constexpr static const char* name = "n_1st";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "n_1st_cuda" }; } // FIXME
  constexpr static int flags = 0;

  static void run(mfields_t mres, mparticles_t mprts)
  {
    cuda_mparticles *cmprts = mprts->cmprts();
    cuda_mfields *cmres = mres->cmflds;
    
    mres->zero();
    cuda_moments_yz_n_1st(cmprts, cmres);
  }
};

template<typename Moment_t>
struct MomentCuda : ItemMomentCRTP<MomentCuda<Moment_t>, PscMfieldsCuda>
{
  using Base = ItemMomentCRTP<MomentCuda<Moment_t>, PscMfieldsCuda>;
  using mfields_t = PscMfieldsCuda;
  using mparticles_t = PscMparticlesCuda;
  
  constexpr static const char* name() { return Moment_t::name; }
  constexpr static int n_comps = Moment_t::n_comps;
  constexpr static fld_names_t fld_names() { return Moment_t::fld_names(); }
  constexpr static int flags = Moment_t::flags;

  MomentCuda(MPI_Comm comm, PscBndBase bnd)
    : Base(comm),
      bnd_(bnd)
  {}

  void run(mparticles_t mprts)
  {
    mfields_t mres{this->mres_};

    Moment_t::run(mres, mprts);
    bnd_.add_ghosts(mres.mflds(), 0, mres->n_comps());
  }

private:
  PscBndBase bnd_;
};

FieldsItemOps<ItemMoment2<MomentCuda<Moment_rho_1st_nc_cuda>>> psc_output_fields_item_rho_1st_nc_cuda_ops;
FieldsItemOps<ItemMoment2<MomentCuda<Moment_n_1st_cuda>>> psc_output_fields_item_n_1st_cuda_ops;

