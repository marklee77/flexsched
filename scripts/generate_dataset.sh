#!/bin/bash

NUMSERVERS_LIST="64"
NUMSERVICES_LIST="100"
SERVERCPUS_LIST="4"
SERVERMEAN_LIST="0.5"
SERVERCOV_LIST="0.000 0.025 0.050 0.075 0.100 0.125 0.150 0.175 0.200 0.225"
SERVERCOV_LIST="${SERVERCOV_LIST} 0.250 0.275 0.300 0.325 0.350 0.375 0.400"
SERVERCOV_LIST="${SERVERCOV_LIST} 0.425 0.450 0.475 0.500 0.525 0.550 0.575"
SERVERCOV_LIST="${SERVERCOV_LIST} 0.600 0.625 0.650 0.675 0.700 0.725 0.750"
SERVERCOV_LIST="${SERVERCOV_LIST} 0.775 0.800 0.825 0.850 0.875 0.900 0.925"
SERVERCOV_LIST="${SERVERCOV_LIST} 0.950 0.975 1.000"
SLACK_LIST="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"

if [ $# != 1 ]; then
  echo "usage: $0 <basedir name>"
  exit 1
fi

BASEDIR=$1
mkdir -p ${BASEDIR}

SCRIPTDIR=$(dirname $0)

# generate base files
for SEED in {101..200}; do
    for NUMSERVERS in ${NUMSERVERS_LIST}; do
        for NUMSERVICES in ${NUMSERVICES_LIST}; do
            for SERVERCPUS in ${SERVERCPUS_LIST}; do
                for SERVERMEAN in ${SERVERMEAN_LIST}; do
                    for SERVERCOV in ${SERVERCOV_LIST}; do
                        PREFIX="${BASEDIR}/prob-${NUMSERVERS}-${NUMSERVICES}-${SEED}-${SERVERCOV}"
                        BASEFILE="${PREFIX}.base"
                        ${SCRIPTDIR}/genprob-from-specs-synth-hw.py ${SEED} ${NUMSERVERS} ${NUMSERVICES} ${SERVERCPUS} ${SERVERMEAN} ${SERVERCOV} ${SCRIPTDIR}/../notes/task_spec_dist.txt > ${BASEFILE}
                        for SLACK in ${SLACK_LIST}; do
                            PROBFILE="${PREFIX}-${SLACK}.prob"
                            ${SCRIPTDIR}/rescale_problem.py ${BASEFILE} ${PROBFILE} "1.0 ${SLACK}" "1.0 0.0"
                        done
                    done
                done
            done
        done
    done
done

#echo "Number of base files (each file = one problem instance): $NUMINSTANCES">> $BASEDIR/README
#echo "">> $BASEDIR/README
#echo "File Names:">> $BASEDIR/README
#echo "Each file is named problem-a-b-c-d-e-f-g-h-i-j-k-l.txt with:">> $BASEDIR/README
#echo "a = number of rigid resources ($NUMRIGID)">> $BASEDIR/README
#echo "b = number of fluid resources ($NUMFLUID)">> $BASEDIR/README
#echo "c = number of servers ($NUMSERVERS)">> $BASEDIR/README
#echo "d = number of services ($NUMSERVICES)">> $BASEDIR/README
#echo "e = rigid slack ($RIGIDSLACK)">> $BASEDIR/README
#echo "f = rigid c.v. ($RIGIDCV)">> $BASEDIR/README
#echo "g = mean fluid need ($FLUIDMEAN)">> $BASEDIR/README
#echo "h = fluid c.v. ($FLUIDCV)">> $BASEDIR/README
#echo "i = sla prob  ($SLAPROB)">> $BASEDIR/README
#echo "j = sla value  ($SLAVALUE)">> $BASEDIR/README
#echo "k = rng seed used">> $BASEDIR/README
#echo "l = instance number">> $BASEDIR/README
#echo "">> $BASEDIR/README
#echo "">> $BASEDIR/README
#echo "">> $BASEDIR/README
#echo "File Content:">> $BASEDIR/README
#echo "<number of rigid needs>">> $BASEDIR/README
#echo "<number of fluid needs>">> $BASEDIR/README
#echo "<number of servers>">> $BASEDIR/README
#echo "<number of services>">> $BASEDIR/README
#echo "<sla> <rigid 1 of svc 1> <rigid 2 of svc 1> ... <fluid 1 of svc 1> <fluid 2 of svc 1>...">> $BASEDIR/README
#echo "<sla> <rigid 1 of svc 2> <rigid 2 of svc 2> ... <fluid 1 of svc 2> <fluid 2 of svc 2>...">> $BASEDIR/README
#echo "....">> $BASEDIR/README
