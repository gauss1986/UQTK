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
# @ node = 1
# @ tasks_per_node = 128
# @ queue
#
#===================================

module load gcc
./nDuffing.x
