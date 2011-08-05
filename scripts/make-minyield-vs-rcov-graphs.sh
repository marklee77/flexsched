#!/bin/bash

DIR=$1
ALGO_LIST="METAHVP RRNZ METAGREEDY METAGREEDYLIGHT METAVP VP_PP_DSUM_W1 HVP_PP_DSUM_ASUM_W1"
ALGO_LIST="METAHVP METAVP METAGREEDY RRNZ"
LINETYPE_LIST="1 3 2 4"
POINTTYPE_LIST="6 4 10 12"
POINTCOLOR_LIST="6 5 10 11"


SCRIPTDIR=$(dirname $0)

for FILE in minyield-vs-cov-*-0.?.txt; do
    PREFIX=$(basename ${FILE} .txt)
    IFS='-' read JUNK1 JUNK2 JUNK3 NUMHOSTS NUMJOBS SLACK <<< "${PREFIX}"
    read -a HEADERS <<< "$(head -1 ${FILE})"
    ALGOS=(${ALGO_LIST})
    LINETYPES=(${LINETYPE_LIST})
    POINTTYPES=(${POINTTYPE_LIST})
    POINTCOLORS=(${POINTCOLOR_LIST})
    ALGO_COLS=()
    IDX=0
    while [ ${IDX} -lt ${#ALGOS[*]} ]; do
        ALGO=${ALGOS[${IDX}]}
        COL=0
        for HEADER in "${HEADERS[@]}"; do
            if [ "${ALGO}" = "${HEADER}" ]; then
                ALGO_COLS[${IDX}]=${COL}
                IDX=$((${IDX}+1))
                break
            fi
            COL=$((${COL}+1))
        done
        if [ ${COL} -ge ${#HEADERS[*]} ]; then
            unset ALGOS[${IDX}] 
            ALGOS=("${ALGOS[@]}")
            unset LINETYPES[${IDX}]
            LINETYPES=("${LINETYPES[@]}")
            unset POINTTYPES[${IDX}]
            POINTTYPES=("${POINTTYPES[@]}")
            unset POINTCOLORS[${IDX}]
            POINTCOLORS=("${POINTCOLORS[@]}")
        fi
    done

    BASE_ALGO="${ALGOS[0]}"
    ALGOS=("${ALGOS[@]}")
    MAX_ALGO_IDX=$((${#ALGOS[*]}-1))
    BASE_COL="${ALGO_COLS[0]}"

    # could also do in BASH with read -a...
    PERL_CODE="print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$F[3]\""
    ALLFAIL_CODE="(\$. == 1) ? \"allfail\" : ("
    for COL in "${ALGO_COLS[@]}"; do
        PERL_CODE="${PERL_CODE}, \"\\t\", \$F[${COL}], \"\\t\", "
        PERL_CODE="${PERL_CODE}(\$F[${BASE_COL}] < 0.0 or \$F[${COL}] < 0.0) "
        PERL_CODE="${PERL_CODE}?  \"?\" : \"\", (\$. == 1) ? \$F[${COL}] : "
        PERL_CODE="${PERL_CODE}\$F[${BASE_COL}] - \$F[${COL}]"
        ALLFAIL_CODE="${ALLFAIL_CODE}\$F[${COL}] < 0.0 and "
    done
    ALLFAIL_CODE="${ALLFAIL_CODE}1) ? 1 : 0"
    PERL_CODE="${PERL_CODE}, \"\\t\", ${ALLFAIL_CODE}, \"\\n\";"

    PLOTFILE="${PREFIX}-gnuplot.txt"

    perl -ane "${PERL_CODE}" < "${FILE}" > "${PLOTFILE}"

    (
        echo "set term png size 1024,768 font \"Helvetica,16\""
        #echo "set size 0.68,0.68"
        #echo "set term postscript eps color font \"Times-Roman\" 16"

        echo "set key bottom outside horizontal width +1" 

        echo "set xlabel \"coefficient of variation\""
        echo "set xrange[0.0:0.9]"

        echo -n "set title \"Failure Rate vs. Specified Coefficient of "
        echo -n "Varation\\nfor ${NUMHOSTS} hosts, ${NUMJOBS} jobs, "
        echo "slack = ${SLACK}\" enhanced" 

        echo "set yrange [0.0:100.0]"
        echo "set ylabel \"failure rate (%)\""

        echo -n "set output \""
        echo -n "failrate-vs-cov-hosts-${NUMHOSTS}-jobs-${NUMJOBS}-slack-${SLACK/./}.png"
        echo "\""

        FAILCOL=$((7+2*${MAX_ALGO_IDX}))
        echo "set style fill transparent solid 0.25"
        echo -n "plot "
        for ((IDX=${MAX_ALGO_IDX}; IDX>=0; IDX--)); do
            PLOTCOL=$((5+2*${IDX}))
            echo -n "\"${PLOTFILE}\" using 2:"
            echo -n "((\$${PLOTCOL} < 0.0) ? 100.0 : 0.0) "
            echo -n "smooth unique with filledcurves x1 "
            echo "title \"${ALGOS[${IDX}]//_/ }\", \\"
        done
        echo -n "\"${PLOTFILE}\" using 2:((\$${FAILCOL} > 0) ? 100.0 : 0.0) "
        echo "smooth unique with filledcurves x1 title \"all algorithms\""

        echo -n "set title \"Minyield vs. Actual Coefficient of Variation\\n"
        echo -n "for ${NUMHOSTS} hosts, ${NUMJOBS} jobs, slack = ${SLACK}, "
        echo "failure = 0\""

        echo "set yrange [0.0:1.00]"
        echo "set ylabel \"minimum yield\""

        echo -n "set output \""
        echo -n "minyield-vs-rcov-hosts-${NUMHOSTS}-jobs-${NUMJOBS}-slack-${SLACK/./}-zerofails.png"
        echo "\""

        echo -n "plot "
        for ((IDX=0; IDX <= ${MAX_ALGO_IDX}; IDX++)); do
            PLOTCOL=$((5+2*${IDX}))
            PLOTVAL="((\$${PLOTCOL} < 0.0) ? 0.0 : \$${PLOTCOL})"
            echo -n "\"${PLOTFILE}\" using 4:${PLOTVAL} with points "
            echo -n "pt ${POINTTYPES[${IDX}]} lc ${POINTCOLORS[${IDX}]} "
            echo -n "title \"${ALGOS[${IDX}]//_/ }\", "
        done
        for ((IDX=0; IDX <= ${MAX_ALGO_IDX}; IDX++)); do
            PLOTCOL=$((5+2*${IDX}))
            PLOTVAL="((\$${PLOTCOL} < 0.0) ? 0.0 : \$${PLOTCOL})"
            echo -n "\"${PLOTFILE}\" using "
            echo -n "(floor(10.0*column(4)+0.5)/10.0):${PLOTVAL} "
            echo -n "with linespoints smooth unique "
            echo -n "lt ${LINETYPES[${IDX}]} lw 2 lc ${LINETYPES[${IDX}]} "
            echo -n "title \"${ALGOS[${IDX}]//_/ } (avg.)\""
            if [ ${IDX} -lt ${MAX_ALGO_IDX} ]; then
                echo ", \\"
            fi
        done
        echo

        echo -n "set title \"Minyield Difference from ${BASE_ALGO} vs. "
        echo -n "Actual Coefficient of Variation\\nfor ${NUMHOSTS} hosts, "
        echo "${NUMJOBS} jobs, slack = ${SLACK}, no failures\""

        echo "set yrange [-0.35:1.00]"
        echo "set ylabel \"minimum yield difference\""

        echo -n "set output \""
        echo -n "minyield-diff-vs-rcov-hosts-${NUMHOSTS}-jobs-${NUMJOBS}-slack-${SLACK/./}-nofails.png"
        echo "\""

        echo -n "plot "
        for ((IDX=0; IDX <= ${MAX_ALGO_IDX}; IDX++)); do
            PLOTCOL=$((6+2*${IDX}))
            echo -n "\"${PLOTFILE}\" using 4:${PLOTCOL} with points "
            echo -n "pt ${POINTTYPES[${IDX}]} lc ${POINTCOLORS[${IDX}]} "
            echo -n "title \"${ALGOS[${IDX}]//_/ }\", "
        done
        for ((IDX=0; IDX <= ${MAX_ALGO_IDX}; IDX++)); do
            PLOTCOL=$((6+2*${IDX}))
            echo -n "\"${PLOTFILE}\" using "
            echo -n "(floor(10.0*column(4)+0.5)/10.0):${PLOTCOL} "
            echo -n "with linespoints smooth unique "
            echo -n "lt ${LINETYPES[${IDX}]} lw 2 lc ${LINETYPES[${IDX}]} "
            echo -n "title \"${ALGOS[${IDX}]//_/ } (avg.)\""
            if [ ${IDX} -lt ${MAX_ALGO_IDX} ]; then
                echo ", \\"
            fi
        done
        echo

    ) | gnuplot
    rm ${PREFIX}-gnuplot.txt
done
