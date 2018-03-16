
#include "psc_sort_impl.hxx"
#include "sort.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"

// ======================================================================
// PscSortNone
//
// we could wrap SortNone as usual, but there's really no point to
// get_as() / put_as() ever

class PscSortNone : public SortBase
{
public:
  void run(PscMparticlesBase mprts_base) override {}
};

// ======================================================================
// psc_sort: subclass "countsort_single"

struct psc_sort_ops_countsort_single : psc_sort_ops {
  using PscSort = PscSortWrapper<SortConvert<SortCountsort<MparticlesSingle>>>;
  psc_sort_ops_countsort_single() {
    name                  = "countsort_single";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
  }
} psc_sort_countsort_single_ops;

// ======================================================================
// psc_sort: subclass "countsort2_single"

struct psc_sort_ops_countsort2_single : psc_sort_ops {
  using PscSort = PscSortWrapper<SortConvert<SortCountsort2<MparticlesSingle>>>;
  psc_sort_ops_countsort2_single() {
    name                  = "countsort2_single";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
  }
} psc_sort_countsort2_single_ops;

// ======================================================================
// psc_sort: subclass "countsort_double"

struct psc_sort_ops_countsort_double : psc_sort_ops {
  using PscSort = PscSortWrapper<SortConvert<SortCountsort<MparticlesDouble>>>;
  psc_sort_ops_countsort_double() {
    name                  = "countsort_double";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
  }
} psc_sort_countsort_double_ops;

// ======================================================================
// psc_sort: subclass "countsort2_double"

struct psc_sort_ops_countsort2_double : psc_sort_ops {
  using PscSort = PscSortWrapper<SortConvert<SortCountsort2<MparticlesDouble>>>;
  psc_sort_ops_countsort2_double() {
    name                  = "countsort2_double";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
  }
} psc_sort_countsort2_double_ops;

// ======================================================================
// psc_sort: subclass "none"

struct psc_sort_ops_none : psc_sort_ops {
  using PscSort = PscSortWrapper<PscSortNone>;
  psc_sort_ops_none() {
    name                  = "none";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
  }
} psc_sort_none_ops;



