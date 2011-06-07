#!/bin/bash
if [ $# != 13 ]; then
    echo "usage:"
    echo -n "$0 <basedir> <num instances> <num rigid> <num fluid> <num servers>"
    echo " <num service list> <rigid res. sigma list> <fluid res. sigma list> <rigid slack list> <fluid mean list> <need cv list> <prob sla list> <slavalue list>"
    exit 1
fi

BASEDIR=$1
NUMINSTANCES=$2
NUMRIGID=$3
NUMFLUID=$4
NUMSERVERS=$5
NUMSERVICES_LIST=$6
RIGIDRESSIGMA_LIST=$7
FLUIDRESSIGMA_LIST=$8
RIGIDSLACK_LIST=$9
FLUIDMEAN_LIST=${10}
NEEDCV_LIST=${11}
PROBSLA_LIST=${12}
SLAVALUE_LIST=${13}

mkdir -p $BASEDIR
cd $BASEDIR

SEED=0

for NUMSERVICES in $NUMSERVICES_LIST; do
    for RIGIDRESSIGMA in $RIGIDRESSIGMA_LIST; do
        for FLUIDRESSIGMA in $FLUIDRESSIGMA_LIST; do
            for RIGIDSLACK in $RIGIDSLACK_LIST; do
                for FLUIDMEAN in $FLUIDMEAN_LIST; do
                    for NEEDCV in $NEEDCV_LIST; do
                        for PROBSLA in $PROBSLA_LIST; do
                            for SLAVALUE in $SLAVALUE_LIST; do
                                for (( INSTANCE=1; $INSTANCE <= $NUMINSTANCES; INSTANCE++ )); do
                                    PREFIX="problem-${NUMRIGID}-${RIGIDRESSIGMA}-${NUMFLUID}-${FLUIDRESSIGMA}-${NUMSERVERS}-${NUMSERVICES}-${RIGIDSLACK}-${FLUIDMEAN}-${NEEDCV}-${PROBSLA}-${SLAVALUE}-${SEED}-${INSTANCE}"
                                    ../../scripts/genprob.py $SEED $NUMRIGID $NUMFLUID $NUMSERVERS $NUMSERVICES $RIGIDRESSIGMA $FLUIDRESSIGMA $RIGIDSLACK $FLUIDMEAN $NEEDCV $PROBSLA $SLAVALUE > ${PREFIX}.prob.txt
                                    SEED="$(($SEED + 1))"
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done
