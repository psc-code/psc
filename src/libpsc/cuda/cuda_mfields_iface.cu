
#include "cuda_iface.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include "psc_fields_cuda.h"

#undef dprintf
#if 0
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

MfieldsCuda::MfieldsCuda(const Grid_t& grid, int n_fields, Int3 ibn)
  : MfieldsBase(grid, n_fields, ibn)
{
  dprintf("CMFLDS: ctor\n");
  cmflds_ = new cuda_mfields(grid, n_fields, ibn);
}

MfieldsCuda::~MfieldsCuda()
{
  dprintf("CMFLDS: dtor\n");
  delete cmflds_;
}

void MfieldsCuda::reset(const Grid_t& new_grid)
{
  dprintf("CMFLDS: reset\n");
  MfieldsBase::reset(new_grid);
  Int3 ibn = -cmflds()->ib;
  int n_fields = cmflds()->n_fields;
  delete cmflds_;
  cmflds_ = new cuda_mfields(new_grid, n_fields, ibn);
}

fields_single_t MfieldsCuda::get_host_fields()
{
  dprintf("CMFLDS: get_host_fields\n");
  return cmflds()->get_host_fields();
}

void MfieldsCuda::copy_to_device(int p, fields_single_t h_flds, int mb, int me)
{
  dprintf("CMFLDS: copy_to_device\n");
  cmflds()->copy_to_device(p, h_flds, mb, me);
}

void MfieldsCuda::copy_from_device(int p, fields_single_t h_flds, int mb, int me)
{
  dprintf("CMFLDS: copy_from_device\n");
  cmflds()->copy_from_device(p, h_flds, mb, me);
}

void MfieldsCuda::axpy_comp_yz(int ym, float a, MfieldsCuda& mflds_x, int xm)
{
  dprintf("CMFLDS: axpy_comp_yz\n");
  cmflds()->axpy_comp_yz(ym, a, mflds_x.cmflds(), xm);
}

void MfieldsCuda::zero_comp(int m)
{
  dprintf("CMFLDS: zero_comp\n");
  assert(!grid_->isInvar(1));
  assert(!grid_->isInvar(2));
  if (grid_->isInvar(0)) {
    cmflds()->zero_comp(m, dim_yz{});
  } else {
    cmflds()->zero_comp(m, dim_xyz{});
  }
}

void MfieldsCuda::zero()
{
  dprintf("CMFLDS: zero\n");
  for (int m = 0; m < cmflds()->n_fields; m++) {
    zero_comp(m);
  }
}

int MfieldsCuda::index(int m, int i, int j, int k, int p) const
{
  return cmflds_->index(m, i, j, k, p);
}

void MfieldsCuda::write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names)
{
  auto& mflds = get_as<MfieldsSingle>(0, n_comps());
  mflds.write_as_mrc_fld(io, name, comp_names);
  put_as(mflds, 0, 0);
}

MfieldsCuda::Patch::Patch(MfieldsCuda& mflds, int p)
  : mflds_(mflds), p_(p)
{}

MfieldsCuda::Accessor MfieldsCuda::Patch::operator()(int m, int i, int j, int k)
{
  return { mflds_, mflds_.index(m, i,j,k, p_) };
}

MfieldsCuda::Accessor::Accessor(MfieldsCuda& mflds, int idx)
  : mflds_(mflds), idx_(idx)
{}

MfieldsCuda::Accessor::operator real_t() const
{
  return mflds_.cmflds_->get_value(idx_);
}

MfieldsCuda::real_t MfieldsCuda::Accessor::operator=(real_t val)
{
  mflds_.cmflds_->set_value(idx_, val);
  return val;
}

MfieldsCuda::real_t MfieldsCuda::Accessor::operator+=(real_t val)
{
  val += mflds_.cmflds_->get_value(idx_);
  mflds_.cmflds_->set_value(idx_, val);
  return val;
}  