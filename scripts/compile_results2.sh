#/bin/bash
INDIR=$1

INDIR=$1

if [ ! -d "${INDIR}" ]; then
    INDIR=.
fi

SCRIPTDIR=$(dirname $0)

HEADERFLAG=0

for RESULT_FILE in ${INDIR}/prob-*.result; do

    if [ "${HEADERFLAG}" -eq 0 ]; then
        HEADERFLAG=1
        echo -n -e "key\thosts\tjobs\tindex\tslack\tcov\trcov"
        perl -F"\|" -ane 'print "$F[0]\n";' < ${RESULT_FILE} | sort | perl -ne 'chomp; print "\t$_";'
        echo
    fi

    RESULT_PREFIX=$(basename ${RESULT_FILE} .result)
    IFS='-' read JUNK NUMHOSTS NUMJOBS INDEX COV SLACK <<< "${RESULT_PREFIX}"
    KEY="${SLACK:2}${COV:0:1}${COV:2}${INDEX}"
    REALCOV=$(${SCRIPTDIR}/calc_cpu_mem_covs.py ${INDIR}/${RESULT_PREFIX}.prob)
    echo -n -e "${KEY}\t${NUMHOSTS}\t${NUMJOBS}\t${INDEX}\t${SLACK}\t${COV}\t${REALCOV}"
    perl -F"\|" -ane 'print "$F[0]\t$F[1]\n";' < ${RESULT_FILE} | sort -k1 | perl -ane 'print "\t$F[1]";'
    echo

done
