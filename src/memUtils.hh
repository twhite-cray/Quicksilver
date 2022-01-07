/// \file
/// Wrappers for memory allocation.

#ifndef MEMUTILS_HH
#define MEMUTILS_HH

#include <cstdlib>

static void* qsCalloc(size_t num, size_t iSize)
{
   return std::calloc(num, iSize);
}

static void qsFree(void* ptr)
{
   std::free(ptr);
}
#endif
