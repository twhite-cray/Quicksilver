#!/bin/bash
module load rocm
module list
set -x
I=${1}
J=${2}
K=${3}
P=40
E=8
X=$(( I * E ))
Y=$(( J * E ))
Z=$(( K * E ))
TASKS=$(( I * J * K ))
NODES=$(( ( TASKS + 7 ) / 8 ))
ELEMENTS=$(( X * Y * Z ))
PARTICLES=$(( ELEMENTS * P ))
export MPICH_GPU_SUPPORT_ENABLED=1
ulimit -c unlimited
rm -f core
srun -l --unbuffered -t 5:00 --exclusive -N ${NODES} -n ${TASKS} -c 8 ../src/qs -i Coral2_P1.inp -X ${X} -Y ${Y} -Z ${Z} -x ${X} -y ${Y} -z ${Z} -I ${I} -J ${J} -K ${K} --nParticles ${PARTICLES} 
