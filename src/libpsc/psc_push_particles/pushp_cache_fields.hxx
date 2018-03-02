
#pragma once

template<typename fields_t, typename dim>
struct CacheFields;

template<typename fields_t>
struct CacheFields<fields_t, dim_yz>
{
  static fields_t from_em(fields_t flds)
  {
    fields_t fld = fields_t(flds.ib(), flds.im(), 9);
    Fields3d<fields_t> F(flds), F_CACHE(fld);
    // FIXME, can do -1 .. 2? NO!, except maybe for 1st order
    // Has to be at least -2 .. +3 because of staggering
    // FIXME, get rid of caching since it's no different from the actual
    // fields...
    for (int iz = fld.ib_[2]; iz < fld.ib_[2] + fld.im_[2]; iz++) {
      for (int iy = fld.ib_[1]; iy < fld.ib_[1] + fld.im_[1]; iy++) {
	F_CACHE(EX, 0,iy,iz) = F(EX, 0,iy,iz);
	F_CACHE(EY, 0,iy,iz) = F(EY, 0,iy,iz);
	F_CACHE(EZ, 0,iy,iz) = F(EZ, 0,iy,iz);
	F_CACHE(HX, 0,iy,iz) = F(HX, 0,iy,iz);
	F_CACHE(HY, 0,iy,iz) = F(HY, 0,iy,iz);
	F_CACHE(HZ, 0,iy,iz) = F(HZ, 0,iy,iz);
      }
    }
    return fld;
  }

  static void to_j(fields_t fld, fields_t flds)
  {
    Fields3d<fields_t> F(flds), F_CACHE(fld);
    for (int iz = fld.ib_[2]; iz < fld.ib_[2] + fld.im_[2]; iz++) {
      for (int iy = fld.ib_[1]; iy < fld.ib_[1] + fld.im_[1]; iy++) {
	F(JXI, 0,iy,iz) += F_CACHE(JXI, 0,iy,iz);
	F(JYI, 0,iy,iz) += F_CACHE(JYI, 0,iy,iz);
	F(JZI, 0,iy,iz) += F_CACHE(JZI, 0,iy,iz);
      }
    }
    fld.dtor();
  }
};

template<typename fields_t, typename dim>
struct CacheFieldsNone
{
  static fields_t from_em(fields_t flds)
  {
    return flds;
  }

  static void to_j(fields_t flds_cache, fields_t flds)
  {}
};

