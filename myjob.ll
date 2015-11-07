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
./Duffing.x -s 0
mv *.dat ./test_Oct31/n10/s00/
./Duffing.x -s 0.1
mv *.dat ./test_Oct31/n10/s01/
./Duffing.x -s 0.2
mv *.dat ./test_Oct31/n10/s02/
./Duffing.x -s 0.3
mv *.dat ./test_Oct31/n10/s03/
./Duffing.x -s 0.4
mv *.dat ./test_Oct31/n10/s04/
./Duffing.x -s 0.5
mv *.dat ./test_Oct31/n10/s05/
./Duffing.x -s 0.6
mv *.dat ./test_Oct31/n10/s06/
./Duffing.x -s 0.7
mv *.dat ./test_Oct31/n10/s07/
./Duffing.x -s 0.8
mv *.dat ./test_Oct31/n10/s08/
./Duffing.x -s 0.9
mv *.dat ./test_Oct31/n10/s09/
./Duffing.x -s 1.0
mv *.dat ./test_Oct31/n10/s10/
