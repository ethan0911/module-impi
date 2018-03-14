#!/bin/bash

#
# For Stampede2 (requires Intel MPI)
#
#
# 1st argument: <# of nodes>
# 2nd argument: <# of threads/node>
# 3rd argument: <ospray build directory>
#

NN=$1 # number of nodes
NT=$2 # number of threads per node

OSPRAY_DIR=$3
SCRIPT_DIR=$4

JOBID=$(date +%s)
echo "job id = ${JOBID}"

cat > submit-n${NN}.sh <<EOF
#!/bin/bash
#SBATCH -n ${NN}
#SBATCH -N ${NN}
#SBATCH --time=4:00:00 # walltime, abbreviated by -t

module load gcc/7.1.0
module load impi
module load hdf5

export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=${NT}
export KMP_AFFINITY=verbose,granularity=core,compact,1,0
export OSPRAY_THREADS=${NT}

date

echo "----------------------------------------------"
echo " Method octant"
echo "----------------------------------------------"
echo 
IMPI_AMR_METHOD=octant \
OSPRAY_AMR_METHOD=current \
mpirun -np $NN \
bash $SCRIPT_DIR $OSPRAY_DIR \
-o octant.${JOBID} \
--osp:mpi
echo 
echo "----------------------------------------------"

echo "----------------------------------------------"
echo " Method current"
echo "----------------------------------------------"
echo 
IMPI_AMR_METHOD=current \
OSPRAY_AMR_METHOD=current \
mpirun -np $NN \
bash $SCRIPT_DIR $OSPRAY_DIR \
-o current.${JOBID} \
--osp:mpi
echo 
echo "----------------------------------------------"

echo "----------------------------------------------"
echo " Method octant pt"
echo "----------------------------------------------"
echo 
IMPI_AMR_METHOD=octant \
OSPRAY_AMR_METHOD=current \
mpirun -np $NN \
bash $SCRIPT_DIR $OSPRAY_DIR \
-o octant.pt.${JOBID} -renderer pt \
--osp:mpi
echo 
echo "----------------------------------------------"

echo "----------------------------------------------"
echo " Method current pt"
echo "----------------------------------------------"
echo 
IMPI_AMR_METHOD=current \
OSPRAY_AMR_METHOD=current \
mpirun -np $NN \
bash $SCRIPT_DIR $OSPRAY_DIR \
-o current.pt.${JOBID} -renderer pt \
--osp:mpi
echo 
echo "----------------------------------------------"

EOF

sbatch -p normal -A OSPRay submit-n${NN}.sh
