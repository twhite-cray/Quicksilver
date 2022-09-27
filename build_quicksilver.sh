#!/bin/bash
module load rocm/5.1.0
module list
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
set -x
export CXX='hipcc'
export CXXFLAGS='-g -O3 -Wall --offload-arch=gfx90a -Werror -ferror-limit=1 -munsafe-fp-atomics -ffast-math'
export CPPFLAGS="-I${MPICH_DIR}/include"
export LDFLAGS="-L${MPICH_DIR}/lib -lmpi ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}"
cd src
make clean
make -j
