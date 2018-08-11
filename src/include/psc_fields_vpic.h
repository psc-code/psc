
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "fields3d.hxx"
#include "fields_traits.hxx"

#include <psc_method.h>

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

struct fields_vpic_t : fields3d<float, LayoutAOS>
{
  using Base = fields3d<float, LayoutAOS>;

  using Base::Base;
};

struct MfieldsHydroVpic : MfieldsBase
{
  using real_t = float;
  using fields_t = fields_vpic_t;

  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic(const Grid_t& grid, int n_fields, Int3 ibn)
    : MfieldsBase(grid, n_fields, ibn)
  {
    assert(grid.n_patches() == 1);
    assert((ibn == Int3{ 1, 1, 1 }));
    assert(n_fields == N_COMP);

    Simulation* sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    vmflds_hydro = sim->hydro_array_;
  }
  
  fields_vpic_t operator[](int p)
  {
    int ib[3], im[3];
    float* data = vmflds_hydro->getData(ib, im);
    return fields_vpic_t(ib, im, N_COMP, data);
  }

  void zero_comp(int m) override { assert(0); }
  void set_comp(int m, double val) override { assert(0); }
  void scale_comp(int m, double val) override { assert(0); }
  void axpy_comp(int m_y, double alpha, MfieldsBase& x_base, int m_x) override { assert(0); }
  void copy_comp(int mto, MfieldsBase& from, int mfrom) override { assert(0); }
  double max_comp(int m) override { assert(0); }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

 public:
  HydroArray *vmflds_hydro;
};

template<>
struct Mfields_traits<MfieldsHydroVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

struct MfieldsStateVpic : MfieldsStateBase
{
  using real_t = float;
  using fields_t = fields_vpic_t;
  
  enum {
    EX = 0,
    EY = 1,
    EZ = 2,
    DIV_E_ERR = 3,
    BX = 4,
    BY = 5,
    BZ = 6,
    DIV_B_ERR = 7,
    TCAX = 8,
    TCAY = 9,
    TCAZ = 10,
    RHOB = 11,
    JFX = 12,
    JFY = 13,
    JFZ = 14,
    RHOF = 15,
    N_COMP = 20,
  };

  MfieldsStateVpic(const Grid_t& grid)
    : MfieldsStateBase{grid, N_COMP, {1, 1, 1}}
  {
    assert(grid.n_patches() == 1);

    Simulation* sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
    sim->field_array_ = FieldArray::create(sim->grid_, sim->material_list_, 0.);
    vmflds_fields_ = sim->field_array_;
  }

  fields_vpic_t operator[](int p)
  {
    assert(p == 0);
    int ib[3], im[3];
    float* data = vmflds().getData(ib, im);
    return fields_vpic_t(ib, im, N_COMP, data);
  }

  FieldArray& vmflds() { return *vmflds_fields_; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

private:
  FieldArray* vmflds_fields_;
};

template<>
struct Mfields_traits<MfieldsStateVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
