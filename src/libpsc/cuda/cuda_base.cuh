
#pragma once

// ----------------------------------------------------------------------
// cuda_base

void cuda_base_init(void);

void* myCudaMalloc(size_t len);
void myCudaFree(void *ptr);
