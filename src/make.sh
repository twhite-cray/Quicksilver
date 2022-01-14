#!/bin/bash
module load craype-accel-amd-gfx90a
module load rocm
module list
set -x
export CXX='hipcc'
export CXXFLAGS='-g -ggdb -O3 -Wall --offload-arch=gfx90a -Werror -ferror-limit=1 -fgpu-rdc -munsafe-fp-atomics'
export CPPFLAGS="-DHAVE_MPI -DHAVE_ASYNC_MPI -I${MPICH_DIR}/include"
export CPPFLAGS="-DHAVE_CUDA -Icuda2hip ${CPPFLAGS}"
export LDFLAGS="-L${MPICH_DIR}/lib -lmpi ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}"
make clean
make -j
