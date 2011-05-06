#!/bin/bash
if [ $# != 8 ]; then
    echo "usage:"
    echo -n "$0 <basedir> <instances> <hosts> <tasks-list> <slack-list>"
    echo " <memvar-list> <cpumedian-list> <cpuvar-list>"
    exit 1
fi

BASEDIR=$1
INSTANCES=$2
HOSTS=$3
TASKS_LIST=$4
SLACK_LIST=$5
MEMVAR_LIST=$6
CPUMEDIAN_LIST=$7
CPUVAR_LIST=$8

if [ -e "$BASEDIR" ]; then
    echo "$BASEDIR already exists!"
    exit 1
fi

mkdir $BASEDIR
cd $BASEDIR

SEED=0

for TASKS in $TASKS_LIST; do
    for SLACK in $SLACK_LIST; do
        for MEMVAR in $MEMVAR_LIST; do
            for CPUMEDIAN in $CPUMEDIAN_LIST; do
                for CPUVAR in $CPUVAR_LIST; do
                    for (( INSTANCE=1; $INSTANCE <= $INSTANCES; INSTANCE++ )); do
                        PREFIX="problem-${HOSTS}-${TASKS}-${SLACK}-${MEMVAR}-${CPUMEDIAN}-${CPUVAR}-${INSTANCE}-${SEED}"
                        ../genparallelprob $SEED $HOSTS $TASKS $SLACK $MEMVAR $CPUMEDIAN $CPUVAR > ${PREFIX}.prob.txt
                        SEED="$(($SEED + 1))"
                    done
                done
            done
        done
    done
done
