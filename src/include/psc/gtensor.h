
#pragma once

#ifdef PSC_HAVE_RMM
#include <rmm/mr/device/thrust_allocator_adaptor.hpp>

#define GTENSOR_DEFAULT_DEVICE_ALLOCATOR(T) rmm::mr::thrust_allocator<T>
#endif

#include <gtensor/gtensor.h>
