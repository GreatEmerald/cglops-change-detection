#!/bin/sh
# Snippet for submitting jobs on the cluster; edit as necessary
SPARK_HOME=/usr/hdp/current/spark2-client/ SPARK_MAJOR_VERSION=2 SPARKR_BACKEND_CONNECTION_TIMEOUT=1209600 spark-submit --master yarn --executor-memory 1920m --driver-memory 1920m --executor-cores=4 --conf spark.memory.fraction=0.01 --deploy-mode cluster --conf spark.shuffle.service.enabled=true --conf spark.dynamicAllocation.enabled=true $@
