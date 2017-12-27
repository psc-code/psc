
#ifndef TEST_INTERPOLATOR_BASE_H
#define TEST_INTERPOLATOR_BASE_H

#include "test_GridBase.h"
#include "test_InterpolatorBase.h"

template<typename InterpolatorBase>
InterpolatorBase* test_InterpolatorBase_create(typename InterpolatorBase::Grid *g)
{
  return InterpolatorBase::create(g);
}

template<typename InterpolatorBase>
void test_InterpolatorBase_methods(InterpolatorBase* interpolator)
{
  typename InterpolatorBase::Grid* g = interpolator->grid();
  (void) g;
}

template<typename InterpolatorBase>
void test_InterpolatorBase_destroy(InterpolatorBase* fa)
{
  // FIXME
}

template<typename InterpolatorBase>
void test_InterpolatorBase()
{
  auto *g = test_GridBase_create<typename InterpolatorBase::Grid>();
  
  InterpolatorBase *interpolator = test_InterpolatorBase_create<InterpolatorBase>(g);
  test_InterpolatorBase_methods(interpolator);
  test_InterpolatorBase_destroy(interpolator);

  test_GridBase_destroy(g);
}

#endif
