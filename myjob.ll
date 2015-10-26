#!/bin/bash
##===================================
## P7 Load Leveler Submission Script
##===================================
##
## Don't change these parameters unless you really know what you are doing
##
#@ environment = MP_INFOLEVEL=0; MP_USE_BULK_XFER=yes; MP_BULK_MIN_MSG_SIZE=64K; \
##                MP_EAGER_LIMIT=64K; MP_DEBUG_ENABLE_AFFINITY=no
##
##===================================
## Avoid core dumps
# @ core_limit   = 0
##===================================
## Job specific
##===================================
#
# @ job_name = Duffing
# @ job_type = parallel
# @ class = verylong
# @ output = $(jobid).out
# @ error = $(jobid).err
# @ wall_clock_limit = 24:00:00
# @ node = 2
# @ tasks_per_node = 128
# @ queue
#
#===================================

module load gcc
rm *.dat
./Duffing.x -n 10 -e 0.1
mv *.dat ./10_gnu/
./Duffing.x -n 20 -e 0.1
mv *.dat ./20_gnu/
./Duffing.x -n 30 -A 2 -e 0.1
mv *.dat ./30_gnu/
./Duffing.x -n 40 -G 1 -A 2 -e 0.1
mv *.dat ./40_gnu/
./Duffing.x -n 50 -G 1 -A 2 -e 0.1
mv *.dat ./50_gnu/
./Duffing.x -n 60 -G 1 -A 2 -e 0.1
mv *.dat ./60_gnu/
./Duffing.x -n 70 -G 1 -A 2 -e 0.1
mv *.dat ./70_gnu/

