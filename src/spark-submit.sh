#!/bin/bash
# Snippet for submitting jobs on the cluster; edit as necessary
# First argument is script to run
# Second argument is memory usage in GiB (optional)
# afterwards options to pass to spark-submit
MemUsageG=${2:-1}
let DrvOverhead=$MemUsageG*1000+300
let GdalCache=$MemUsageG*1000
# +1000
#  --num-executors 20 \
SPARK_HOME=/usr/hdp/current/spark2-client/ SPARK_MAJOR_VERSION=2 SPARKR_BACKEND_CONNECTION_TIMEOUT=1209600 \
  spark-submit --master yarn --executor-memory 600m --driver-memory 600m --executor-cores=4 \
  --conf spark.memory.fraction=0.01 --conf spark.driver.memoryOverhead=${DrvOverhead}m \
  --conf spark.executor.memoryOverhead=${DrvOverhead}m --deploy-mode cluster \
  --conf spark.shuffle.service.enabled=true --conf spark.dynamicAllocation.enabled=true \
  --conf spark.yarn.maxAppAttempts=2 \
  --conf spark.executorEnv.SPARKR_BACKEND_CONNECTION_TIMEOUT=1209600 \
  --conf spark.yarn.appMasterEnv.SPARKR_BACKEND_CONNECTION_TIMEOUT=1209600 \
  --conf spark.executorEnv.CPL_DEBUG=ON --conf spark.executorEnv.GDAL_CACHEMAX=${GdalCache} \
  --conf spark.executorEnv.GDAL_READDIR_LIMIT_ON_OPEN=2000 $1 ${@:3}
