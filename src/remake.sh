#!/bin/bash
module load rocm
module list
set -x
export CXX='hipcc'
export CXXFLAGS='-g -O3 -Wall --offload-arch=gfx90a -Werror -ferror-limit=1'
export CPPFLAGS="-I${MPICH_DIR}/include"
export LDFLAGS="-L${MPICH_DIR}/lib -lmpi ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}"
make
