#/bin/bash
#$ -j y
#$ -o /home/mstillwe/work/flexsched_simulator/batch_output
#$ -v HOME,LD_LIBRARY_PATH
BASEDIR=${HOME}/work/flexsched_simulator
FLEXSCHED="${BASEDIR}/src/flexsched.glpk"
ALGOS="LPBOUND RRND RRNZ METAGREEDY METAGREEDYLIGHT METAHVP"

SORTS="ALEX AMAX ASUM AMAXRATIO AMAXDIFF DLEX DMAX DSUM DMAXRATIO DMAXDIFF"
HVPOPTS="E S S_E R R_E R_S R_S_E"

LOCKDIR=${BASEDIR}/locks

PROBDIR=$1

for PROBFILEPATH in ${PROBDIR}/*.prob; do
    PROBFILE=$(basename ${PROBFILEPATH})
    exec 200>${LOCKDIR}/${PROBFILE}.lck
    flock -xn 200
    if [ 0 = $? ]; then
        echo "Succeeded getting lock for ${PROBFILE}, solving..."
        RESULTFILEPATH=${PROBDIR}/$(basename ${PROBFILE} .prob).result
        #LD_LIBRARY_PATH=${HOME}/progs/base/lib 
        ${FLEXSCHED} "${ALGOS}" ${PROBFILEPATH} ${RESULTFILEPATH}
        for ITMSORT in ${SORTS}; do
            ${FLEXSCHED} "VP_PP_${ITMSORT}_W1 VP_FF_${ITMSORT} VP_BF_${ITMSORT}" ${PROBFILEPATH} ${RESULTFILEPATH}
            for BINSORT in ${SORTS}; do
                ${FLEXSCHED} "HVP_PP_${ITMSORT}_${BINSORT}_W1" ${PROBFILEPATH} ${RESULTFILEPATH}
                ${FLEXSCHED} "HVP_FF_${ITMSORT}_${BINSORT} HVP_BF_${ITMSORT}_${BINSORT}" ${PROBFILEPATH} ${RESULTFILEPATH}
                ${FLEXSCHED} "HVP_FF_${ITMSORT}_${BINSORT}_R HVP_BF_${ITMSORT}_${BINSORT}_R" ${PROBFILEPATH} ${RESULTFILEPATH}
                for HVPOPT in ${HVPOPTS}; do
                    ${FLEXSCHED} "HVP_PP_${ITMSORT}_${BINSORT}_W1_${HVPOPT}" ${PROBFILEPATH} ${RESULTFILEPATH}
                done
            done
        done
        flock -u 200
    else
        echo "Failed to get lock for ${PROBFILE}, trying next..."
    fi
    exec 200>&-
done
