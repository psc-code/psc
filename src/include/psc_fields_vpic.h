
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "fields3d.hxx"
#include "fields_traits.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

struct fields_vpic_t : fields3d<float, LayoutAOS>
{
  using Base = fields3d<float, LayoutAOS>;

  using Base::Base;
};

struct MfieldsVpic : MfieldsBase
{
  using real_t = float;
  using fields_t = fields_vpic_t;
  using Base = MfieldsBase;

  using Base::Base;

  MfieldsVpic(const Grid_t& grid, int n_fields, Int3 ibn);
  ~MfieldsVpic();

  inline fields_vpic_t operator[](int p)
  {
    // FIXME hacky...
    if (n_comps() == VPIC_MFIELDS_N_COMP) {
      int ib[3], im[3];
      float* data = vmflds_fields_->getData(ib, im);
      return fields_vpic_t(ib, im, VPIC_MFIELDS_N_COMP, data);
    } else if (n_comps() == VPIC_HYDRO_N_COMP) {
      int ib[3], im[3];
      float* data = vmflds_hydro->getData(ib, im);
      return fields_vpic_t(ib, im, VPIC_HYDRO_N_COMP, data);
    } else {
      assert(0);
    }
  }

  FieldArray& vmflds() { return *vmflds_fields_; }
  Simulation* sim() { return sim_; }

  void zero_comp(int m) override { assert(0); }
  void set_comp(int m, double val) override { assert(0); }
  void scale_comp(int m, double val) override { assert(0); }
  void axpy_comp(int m_y, double alpha, MfieldsBase& x_base, int m_x) override { assert(0); }
  double max_comp(int m) override { assert(0); }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

 private:
  Simulation* sim_;
  FieldArray* vmflds_fields_;
 public:
  HydroArray *vmflds_hydro;
};

template<>
struct Mfields_traits<MfieldsVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

struct MfieldsStateVpic
{
  using fields_t = MfieldsVpic::fields_t;
  using real_t = MfieldsVpic::real_t;
  
  MfieldsStateVpic(const Grid_t& grid, int n_fields, Int3 ibn)
    : mflds_{grid, n_fields, ibn}
  {}

  int n_patches() const { return mflds_.n_patches(); }
  const Grid_t& grid() { return mflds_.grid(); }
  fields_vpic_t operator[](int p) { return mflds_[p]; }
  Int3 ibn() const { return mflds_.ibn(); }
  FieldArray& vmflds() { return mflds_.vmflds(); }
  Simulation* sim() { return mflds_.sim(); }
  MfieldsBase& base() { return mflds_; }

  void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names)
  {
    int n_comps = comp_names.size(); // FIXME, hack around the material ids, which aren't even float...
    
    mrc_fld* fld = mrc_domain_m3_create(ppsc->mrc_domain_);
    mrc_fld_set_name(fld, name.c_str());
    mrc_fld_set_param_int(fld, "nr_ghosts", 0);
    mrc_fld_set_param_int(fld, "nr_comps", n_comps);
    mrc_fld_setup(fld);
    for (int m = 0; m < n_comps; m++) {
      mrc_fld_set_comp_name(fld, m, comp_names[m].c_str());
    }

    for (int p = 0; p < n_patches(); p++) {
      mrc_fld_patch *m3p = mrc_fld_patch_get(fld, p);
      mrc_fld_foreach(fld, i,j,k, 0,0) {
	for (int m = 0; m < n_comps; m++) {
	  MRC_M3(m3p ,m , i,j,k) = (*this)[p](m, i,j,k);
	}
      } mrc_fld_foreach_end;
      mrc_fld_patch_put(fld);
    }
  
    mrc_fld_write(fld, io);
    mrc_fld_destroy(fld);
  }

  
private:
  MfieldsVpic mflds_;
};

template<>
struct Mfields_traits<MfieldsStateVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
