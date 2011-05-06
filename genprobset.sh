#!/bin/bash
if [ $# != 11 ]; then
    echo "usage:"
    echo -n "$0 <basedir> <num instances> <num rigid> <num fluid> <num servers>"
    echo " <num service list> <rigid slack list> <fluid mean list> <need cv list> <prob sla list> <slavalue list>"
    exit 1
fi

BASEDIR=$1
NUMINSTANCES=$2
NUMRIGID=$3
NUMFLUID=$4
NUMSERVERS=$5
NUMSERVICES_LIST=$6
RIGIDSLACK_LIST=$7
FLUIDMEAN_LIST=$8
NEEDCV_LIST=${9}
PROBSLA_LIST=${10}
SLAVALUE_LIST=${11}

cd $BASEDIR

SEED=0


for NUMSERVICES in $NUMSERVICES_LIST; do
    for RIGIDSLACK in $RIGIDSLACK_LIST; do
        for FLUIDMEAN in $FLUIDMEAN_LIST; do
            for NEEDCV in $NEEDCV_LIST; do
                for PROBSLA in $PROBSLA_LIST; do
                    for SLAVALUE in $SLAVALUE_LIST; do
                        for (( INSTANCE=1; $INSTANCE <= $NUMINSTANCES; INSTANCE++ )); do
                          PREFIX="problem-${NUMRIGID}-${NUMFLUID}-${NUMSERVERS}-${NUMSERVICES}-${RIGIDSLACK}-${FLUIDMEAN}-${NEEDCV}-${PROBSLA}-${SLAVALUE}-${SEED}-${INSTANCE}"
			../genprob.py $SEED $NUMRIGID $NUMFLUID $NUMSERVERS $NUMSERVICES $RIGIDSLACK $FLUIDMEAN $NEEDCV $PROBSLA $SLAVALUE > ${PREFIX}.prob.txt
                        SEED="$(($SEED + 1))"
                        done
                    done
                done
            done
        done
    done
done
