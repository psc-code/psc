
#include "sort_vpic.hxx"

#include "psc_sort_private.h"

#include "vpic_iface.h"
#include "sort.hxx"

// ----------------------------------------------------------------------
// psc_sort: subclass "vpic"

struct psc_sort_ops_vpic : psc_sort_ops {
  using PscSort = PscSortWrapper<SortConvert<SortVpic>>;
  psc_sort_ops_vpic() {
    name                  = "vpic";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
  }
} psc_sort_vpic_ops;


