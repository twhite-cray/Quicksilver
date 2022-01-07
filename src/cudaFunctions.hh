#ifndef CUDAFUNCTIONS_HH
#define CUDAFUNCTIONS_HH

#include "cudaUtils.hh"
#include "DeclareMacro.hh"

#if defined (HAVE_HIP)
void warmup_kernel();
#endif

#endif
