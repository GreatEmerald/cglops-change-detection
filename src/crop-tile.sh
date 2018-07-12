#!/bin/sh
SPARK_HOME=/usr/hdp/current/spark2-client/ SPARK_MAJOR_VERSION=2 spark-submit --master=local[4] detect-breaks.R --crop-only $@
