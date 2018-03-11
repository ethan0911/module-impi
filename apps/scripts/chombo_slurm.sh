#!/bin/bash

#
# For SCI kepler (requires Intel MPI)
#

# 1st argument: <# of nodes>
# 2nd argument: <# of threads/node>
# 3rd argument: <ospray build directory>

NN=$1
NT=$2
OSPRAY_DIR=$3
JOBID=$(date +%s)
cat > submit-n${NN}.sh <<EOF
#!/bin/bash
#SBATCH -n ${NN}
#SBATCH -N ${NN}
#SBATCH --time=4:00:00 # walltime, abbreviated by -t

export OSPRAY_THREADS=${NT}

date

echo "----------------------------------------------"
echo " Method octant"
echo "----------------------------------------------"
echo 
IMPI_AMR_METHOD=octant \
mpirun -np $NN ${OSPRAY_DIR}/ospImplicitIsoSurfaceBench \
/home/sci/feng/Desktop/ws/data/chombo/chombo_amr.osp --osp:mpi \
-iso 0.7 \
-vp 24.684797 17.313093 -10.046009 \
-vu 0.000000 1.000000 0.000000 \
-vi 8.000000 7.999991 7.999997 \
-frames 50 200 \
-o octant.${JOBID}
echo 
echo "----------------------------------------------"

echo "----------------------------------------------"
echo " Method nearest"
echo "----------------------------------------------"
echo 
IMPI_AMR_METHOD=nearest \
mpirun -np $NN ${OSPRAY_DIR}/ospImplicitIsoSurfaceBench \
/home/sci/feng/Desktop/ws/data/chombo/chombo_amr.osp --osp:mpi \
-iso 0.7 \
-vp 24.684797 17.313093 -10.046009 \
-vu 0.000000 1.000000 0.000000 \
-vi 8.000000 7.999991 7.999997 \
-frames 50 200 \
-o nearest.${JOBID}
echo 
echo "----------------------------------------------"

echo "----------------------------------------------"
echo " Method current"
echo "----------------------------------------------"
echo 
IMPI_AMR_METHOD=current \
mpirun -np $NN ${OSPRAY_DIR}/ospImplicitIsoSurfaceBench \
/home/sci/feng/Desktop/ws/data/chombo/chombo_amr.osp --osp:mpi \
-iso 0.7 \
-vp 24.684797 17.313093 -10.046009 \
-vu 0.000000 1.000000 0.000000 \
-vi 8.000000 7.999991 7.999997 \
-frames 50 200 \
-o current.${JOBID}
echo 
echo "----------------------------------------------"

EOF

sbatch submit-n${NN}.sh
