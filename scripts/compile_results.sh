#!/bin/bash

INDIR=$1
OUTDIR=$2

if [ ! -d "${INDIR}" ]; then
    INDIR=.
fi

if [ ! -d "${OUTDIR}" ]; then
    OUTDIR=${INDIR}
fi

SCRIPTDIR=$(dirname $0)

for RESULT_FILE in ${INDIR}/prob-*.result; do
    RESULT_PREFIX=$(basename ${RESULT_FILE} .result)
    IFS='-' read JUNK NUMHOSTS NUMJOBS INDEX COV SLACK <<< "${RESULT_PREFIX}"
    KEY="${COV:0:1}${COV:2}${INDEX}"
    REALCOV=$(${SCRIPTDIR}/calc_cpu_mem_covs.py ${INDIR}/${RESULT_PREFIX}.prob)
    RCOV_FILE=${OUTDIR}/realcovs-${NUMHOSTS}-${NUMJOBS}-${SLACK}.txt
    echo -e "${KEY}\t${COV}\t${INDEX}\t${REALCOV}" >> ${RCOV_FILE}
    while IFS="|" read ALGO MINYIELD JUNK; do
        MINYIELD_FILE=${OUTDIR}/minyields-${NUMHOSTS}-${NUMJOBS}-${ALGO}-${SLACK}.txt
        echo -e "${KEY}\t${MINYIELD}" >> ${MINYIELD_FILE}
    done < ${RESULT_FILE}
done
for RCOV_FILE in ${OUTDIR}/realcovs-*.txt; do
    RCOV_PREFIX=$(basename ${RCOV_FILE} .txt)
    RCOV_SORT_FILE=${OUTDIR}/${RCOV_PREFIX}-sorted.txt
    IFS='-' read JUNK NUMHOSTS NUMJOBS SLACK <<< "${RCOV_PREFIX}"
    echo -e "key\tcov\tindex\trcov" > ${RCOV_SORT_FILE}
    sort ${RCOV_FILE} >> ${RCOV_SORT_FILE}
    rm ${RCOV_FILE}
    ITER=0
    FIELDS="1.1 1.2 1.3"
    JOIN_FILE=${RCOV_SORT_FILE}
    for MINYIELD_FILE in $(ls ${OUTDIR}/minyields-${NUMHOSTS}-${NUMJOBS}-*-${SLACK}.txt | sort); do
        MINYIELD_PREFIX=$(basename ${MINYIELD_FILE} .txt)
        MINYIELD_SORT_FILE=${OUTDIR}/${MINYIELD_PREFIX}-sorted.txt
        IFS='-' read JUNK NUMHOSTS NUMJOBS ALGO JUNK2 <<< "${MINYIELD_PREFIX}"
        echo -e "key\t${ALGO}" > ${MINYIELD_SORT_FILE}
        sort ${MINYIELD_FILE} >> ${MINYIELD_SORT_FILE}
        FIELDS="${FIELDS} 1.$((${ITER}+4))"
        OUT_FILE="${OUTDIR}/minyield-vs-cov-${NUMHOSTS}-${NUMJOBS}-${SLACK}-join-${ITER}.txt"
        join -o "${FIELDS} 2.2" ${JOIN_FILE} ${MINYIELD_SORT_FILE} > ${OUT_FILE}
        rm ${MINYIELD_FILE} ${MINYIELD_SORT_FILE} ${JOIN_FILE}
        ITER=$((${ITER}+1))
        JOIN_FILE=${OUT_FILE}
    done        
    mv ${JOIN_FILE} minyield-vs-cov-${NUMHOSTS}-${NUMJOBS}-${SLACK}.txt
done
