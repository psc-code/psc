
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
struct MomentCuda
{
  using mfields_t = PscMfieldsCuda;
  using mparticles_t = PscMparticlesCuda;
  
  constexpr static const char* name() { return Moment_t::name; }

  MomentCuda(MPI_Comm comm, PscBndBase bnd)
    : bnd_(bnd)
  {
    auto n_comps = Moment_t::n_comps;
    auto fld_names = Moment_t::fld_names();
    assert(n_comps <= POFI_MAX_COMPS);

    if (!Moment_t::flags & POFI_BY_KIND) {
      mres_ = mfields_t::create(comm, ppsc->grid(), n_comps).mflds();
      for (int m = 0; m < n_comps; m++) {
	psc_mfields_set_comp_name(mres_, m, fld_names[m]);
      }
    } else {
      mres_ = mfields_t::create(comm, ppsc->grid(), n_comps * ppsc->nr_kinds).mflds();
      for (int k = 0; k < ppsc->nr_kinds; k++) {
	for (int m = 0; m < n_comps; m++) {
	  auto s = std::string(fld_names[m]) + "_" + ppsc->kinds[k].name;
	  psc_mfields_set_comp_name(mres_, k * n_comps + m, s.c_str());
	}
      }
    }
    psc_mfields_list_add(&psc_mfields_base_list, &mres_);
  }

  ~MomentCuda()
  {
    psc_mfields_list_del(&psc_mfields_base_list, &mres_);
    psc_mfields_destroy(mres_);
  }

  void run(mparticles_t mprts)
  {
    mfields_t mres{mres_};

    Moment_t::run(mres, mprts);
    bnd_.add_ghosts(mres.mflds(), 0, mres->n_comps());
  }

  psc_mfields*& mres() { return mres_; }

private:
  psc_mfields* mres_;
  PscBndBase bnd_;
};

FieldsItemOps<ItemMoment2<MomentCuda<Moment_rho_1st_nc_cuda>>> psc_output_fields_item_rho_1st_nc_cuda_ops;
FieldsItemOps<ItemMoment2<MomentCuda<Moment_n_1st_cuda>>> psc_output_fields_item_n_1st_cuda_ops;

