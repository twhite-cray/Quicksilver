#!/bin/bash
#SBATCH --exclusive
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH -S 8
#SBATCH --time=15:00
#SBATCH --job-name=quicksilver
#SBATCH -o output/job-%j.out

module load rocm/5.1.0

set -x

NODES=${SLURM_NNODES}
DIVS=( ${NODES} )
DMIN=$(( NODES * NODES ))
XYZ=( ${NODES} 1 1 )
SQRT=$( echo "sqrt(${NODES})+1" | bc )
for I in $(seq ${SQRT} -1 1)
do
  if [ "$((NODES%I))" -eq 0 ]
  then
    DIVS+=( ${I} )
  fi
done
for I in "${DIVS[@]}"
do
  for J in "${DIVS[@]}"
  do
    if [ "$(( NODES % (I * J) ))" -eq 0 ]
    then
      K=$(( NODES / (I * J) ))
      D=$(( I * I + J * J + K * K ))
      if [ ${D} -lt ${DMIN} ]
      then
        DMIN=${D}
        XYZ=( ${I} ${J} ${K} )
      fi
    fi
  done
done
echo ${XYZ[@]}
IFS=$'\n'
SORTED=( $( sort -nr <<<"${XYZ[*]}" ) )
unset IFS
px=$(( SORTED[0] * 2 ))
py=$(( SORTED[1] * 2 ))
pz=$(( SORTED[2] * 2 ))
export MPICH_VERSION_DISPLAY=1
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=NUMA

ulimit -c 0
EXE=../../src/qs

mkdir output/job-${SLURM_NNODES}n-${SLURM_JOBID}
cd output/job-${SLURM_NNODES}n-${SLURM_JOBID}

I=${px}
J=${py}
K=${pz}
P=16000
E=10
X=$(( I * E ))
Y=$(( J * E ))
Z=$(( K * E ))
TASKS=$(( 8 * NODES ))
CORES=${SLURM_CPUS_PER_TASK}
ELEMENTS=$(( X * Y * Z ))
PARTICLES=$(( ELEMENTS * P ))
srun -u -N ${NODES} -n ${TASKS} -c ${CORES} --gpus-per-node=8 --gpu-bind=closest bash -c "export HIP_VISIBLE_DEVICES=\${ROCR_VISIBLE_DEVICES} ; unset ROCR_VISIBLE_DEVICES ; ${EXE} -i ../Coral2_P1.inp -X ${X} -Y ${Y} -Z ${Z} -x ${X} -y ${Y} -z ${Z} -I ${I} -J ${J} -K ${K} --nParticles ${PARTICLES}" | tee quicksilver.out 2>&1

module list |& tee launch.modules.log
env |& tee launch.env.log
ulimit -a |& tee launch.ulimits.log
ldd "${EXE}" |& tee launch.ldd.log
