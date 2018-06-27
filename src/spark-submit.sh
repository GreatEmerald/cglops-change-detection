#!/bin/sh
# Snippet for submitting jobs on the cluster; edit as necessary
SPARK_HOME=/usr/hdp/current/spark2-client/ SPARK_MAJOR_VERSION=2 SPARKR_BACKEND_CONNECTION_TIMEOUT=1209600 spark-submit --master yarn --executor-memory 600m --driver-memory 600m --executor-cores=1 --conf spark.memory.fraction=0.01 --conf spark.yarn.driver.memoryOverhead=2300m --conf spark.yarn.executor.memoryOverhead=400m --deploy-mode cluster --conf spark.shuffle.service.enabled=true --conf spark.dynamicAllocation.enabled=true --conf spark.executorEnv.CPL_DEBUG=ON --conf spark.executorEnv.GDAL_CACHEMAX=400 $@
