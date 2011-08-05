#!/bin/bash

DIR=$1
ALGO_LIST="METAHVP RRNZ METAGREEDY METAGREEDYLIGHT METAVP VP_PP_DSUM_W1 HVP_PP_DSUM_ASUM_W1"
ALGO_LIST="METAHVP METAGREEDY METAGREEDYLIGHT METAVP VP_PP_DSUM_W1 HVP_PP_DSUM_ASUM_W1"

ALGOS=(${ALGO_LIST})
MAX_ALGO_IDX=$((${#ALGOS[*]}-1))
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
    grep -v "\-1\.000" ${FILE} > ${PREFIX}-nofails.txt
    perl -pe 's/-1\.000/0.000/g;' < ${FILE} > ${PREFIX}-zerofails.txt
    (
        echo "set xlabel \"classes\""
        echo "set ylabel \"minimum yield difference\""
        echo "set key bottom outside horizontal width +1" 
        echo "set term png size 1024,768"
        if [ $(wc -l < ${PREFIX}-nofails.txt) -gt 1 ]; then
            echo -n "set title \"Minyield Difference from ${ALGOS[0]} vs. "
            echo -n "Machine Classes for ${NUMHOSTS} hosts, "
            echo "${NUMJOBS} jobs, slack = ${SLACK}, no failures\""
            echo -n "set output "
            echo "\"minyield-vs-classes-hosts-${NUMHOSTS}-jobs-${NUMJOBS}-slack-${SLACK}-nofails.png\""
            echo -n "plot "
            for IDX in $(eval "echo {1..${MAX_ALGO_IDX}}"); do
                echo -n "\"${PREFIX}-nofails.txt\" "
                echo -n "using 5:(column(${ALGO_COLS[0]})-column(${ALGO_COLS[${IDX}]})) "
                echo -n "with points title \"${ALGOS[${IDX}]}\","
                echo -n "\"${PREFIX}-nofails.txt\" using "
                echo -n "5:(column(${ALGO_COLS[0]})-column(${ALGO_COLS[${IDX}]})) "
                echo -n "with linespoints smooth unique title \"${ALGOS[${IDX}]}\""
                if [ ${IDX} -lt ${MAX_ALGO_IDX} ]; then
                    echo ", \\"
                fi
            done
            echo
        fi
        if [ $(wc -l < ${PREFIX}-zerofails.txt) -gt 1 ]; then
            echo -n "set title \"Minyield Difference from ${ALGOS[0]} vs. "
            echo -n "Machine Classes for ${NUMHOSTS} hosts, "
            echo "${NUMJOBS} jobs, slack = ${SLACK}, failure = 0.000\""
            echo -n "set output "
            echo "\"minyield-vs-classes-hosts-${NUMHOSTS}-jobs-${NUMJOBS}-slack-${SLACK}-zerofails.png\""
            echo -n "plot "
            for IDX in $(eval "echo {1..${MAX_ALGO_IDX}}"); do
                echo -n "\"${PREFIX}-zerofails.txt\" "
                echo -n "using 5:(column(${ALGO_COLS[0]})-column(${ALGO_COLS[${IDX}]})) "
                echo -n "with points title \"${ALGOS[${IDX}]}\","
                echo -n "\"${PREFIX}-zerofails.txt\" using "
                echo -n "5:(column(${ALGO_COLS[0]})-column(${ALGO_COLS[${IDX}]})) "
                echo -n "with linespoints smooth unique title \"${ALGOS[${IDX}]}\""
                if [ ${IDX} -lt ${MAX_ALGO_IDX} ]; then
                    echo ", \\"
                fi
            done
            echo
        fi
    ) | gnuplot
done
