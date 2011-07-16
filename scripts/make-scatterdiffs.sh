#!/bin/bash

DIR=$1
NUMHOSTS=$2
NUMJOBS=$3
ALGO_LIST="METAHVP RRNZ METAGREEDY METAGREEDYLIGHT METAVP VP_PP_DSUM_W1 HVP_PP_DSUM_ASUM_W1"

ALGOS=(${ALGO_LIST})
NUM_ALGOS=${#ALGOS[*]}
SCRIPTDIR=$(dirname $0)

for FILE in minyield-vs-cov-*-0.?.txt; do
    PREFIX=$(basename ${FILE} .txt)
    IFS='-' read JUNK1 JUNK2 JUNK3 NUMHOSTS NUMJOBS SLACK <<< "${PREFIX}"
    read -a HEADERS <<< "$(head -1 ${FILE})"
    ALGO_COLS=()
    for ALGO in "${ALGOS[@]}"; do
        IDX=1
        for HEADER in "${HEADERS[@]}"; do
            if [ "${ALGO}" = "${HEADER}" ]; then
                ALGO_COLS=("${ALGO_COLS[@]}" "${IDX}")
                break
            fi
            IDX=$((${IDX}+1))
        done
    done
    BASE_COL=${ALGO_COLS[0]}
    LAST_COL=${ALGO_COLS[$((${#ALGO_COLS[*]}-1))]}
    unset ALGO_COLS[0]
    grep -v "\-1\.000" ${FILE} > ${PREFIX}-nofails.txt
    perl -pe 's/-1\.000/0.000/g;' < ${FILE} > ${PREFIX}-zerofails.txt
    (
        echo "set xlabel \"coefficient of variation\""
        echo "set ylabel \"minimum yield difference\""
        echo "set key bottom outside horizontal width +1" 
        echo "set term png size 1024,768"
        if [ $(wc -l < ${PREFIX}-nofails.txt) -gt 1 ]; then
            echo "set title \"Minyield Difference from ${ALGOS[0]} vs. Actual Coefficient of Variation for ${NUMHOSTS} hosts, ${NUMJOBS} jobs, slack = ${SLACK}, no failures\""
            echo "set output \"minyield-vs-rcov-hosts-${NUMHOSTS}-jobs-${NUMJOBS}-slack-${SLACK}-nofails.png\""
            echo -n "plot "
            for COL in "${ALGO_COLS[@]}"; do
                echo -n "\"${PREFIX}-nofails.txt\" using 4:(column(${BASE_COL})-column(${COL})) with points title columnheader(${COL}),"
                echo -n "\"${PREFIX}-nofails.txt\" using (floor(20.0*column(4)+0.5)/20.0):(column(${BASE_COL})-column(${COL})) with linespoints smooth unique title columnheader(${COL})"
                if [ ${COL} -ne ${LAST_COL} ]; then
                    echo -n ","
                fi
                echo "\\"
            done
            echo
        fi
        if [ $(wc -l < ${PREFIX}-zerofails.txt) -gt 1 ]; then
            echo "set title \"Minyield Difference from ${ALGOS[0]} vs. Actual Coefficient of Variation for ${NUMHOSTS} hosts, ${NUMJOBS} jobs, slack = ${SLACK}, failure = 0.000\""
            echo "set output \"minyield-vs-rcov-hosts-${NUMHOSTS}-jobs-${NUMJOBS}-slack-${SLACK}-zerofails.png\""
            echo -n "plot "
            for COL in "${ALGO_COLS[@]}"; do
                echo -n "\"${PREFIX}-zerofails.txt\" using 4:(column(${BASE_COL})-column(${COL})) with points title columnheader(${COL}),"
                echo -n "\"${PREFIX}-zerofails.txt\" using (floor(20.0*column(4)+0.5)/20.0):(column(${BASE_COL})-column(${COL})) with linespoints smooth unique title columnheader(${COL})"
                if [ ${COL} -ne ${LAST_COL} ]; then
                    echo -n ","
                fi
                echo "\\"
            done
            echo
        fi
    ) | gnuplot
done
