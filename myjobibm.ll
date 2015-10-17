#!/bin/bash
##===============================================================================
## Specifies the name of the shell to use for the job 
# @ job_name = DuffingIBM
# @ job_type = parallel
# @ class = verylong
# @ environment = COPY_ALL; MEMORY_AFFINITY=MCM; MP_SYNC_QP=YES; \
#                MP_RFIFO_SIZE=16777216; MP_SHM_ATTACH_THRESH=500000; \
#                MP_EUIDEVELOP=min; MP_USE_BULK_XFER=yes; \
#                MP_RDMA_MTU=4K; MP_BULK_MIN_MSG_SIZE=64k; MP_RC_MAX_QP=8192; \
#                PSALLOC=early; NODISCLAIM=true
# @ node = 1
# @ tasks_per_node = 1
# @ node_usage = not_shared
# @ output = $(jobid).out
# @ error = $(jobid).err
# @ wall_clock_limit = 24:00:00
# @ queue

cd $SCRATCH/UQTK/examples_cpp/Duffing/

# # next variable is for OpenMP
export OMP_NUM_THREADS=128
# # next variable is for ccsm_launch
export THRDS_PER_TASK=128

mkdir 30_IBM/
./Duffing.x -n 30
mv *.dat 30_IBM
mkdir 50_IBM/
./Duffing.x -n 50 -G 1
mv *.dat 50_IBM
mkdir 70_IBM/
./Duffing.x -n 70 -G 1 -P 2
mv *.dat 70_IBM

