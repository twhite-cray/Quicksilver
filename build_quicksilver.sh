#!/bin/bash
module load rocm/5.3.0
module list
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
set -x
export QUIP="${PWD}/quip"
export CXX='hipcc'
export CXXFLAGS='-g -O3 --offload-arch=gfx90a:xnack+ -fgpu-rdc -std=c++11 -ferror-limit=1 -munsafe-fp-atomics -ffast-math'
export CPPFLAGS="-DHAVE_CUDA -DHAVE_MPI -I${MPICH_DIR}/include -I${QUIP}/cuda2hip"
export LDFLAGS="-L${MPICH_DIR}/lib -lmpi ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}"
cd src
make clean
make -j
