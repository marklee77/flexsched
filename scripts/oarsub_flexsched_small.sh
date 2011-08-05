#/bin/bash
#OAR -l /core=1,walltime=30:0:0
#OAR -O /home/mstillwell/work/flexsched_simulator/batch_output/%jobid%.out
#OAR -E /home/mstillwell/work/flexsched_simulator/batch_output/%jobid%.err
HOME=/home/mstillwell
export LD_LIBRARY_PATH=${HOME}/progs/base/lib
BASEDIR=${HOME}/work/flexsched_simulator
FLEXSCHED="${BASEDIR}/src/flexsched.glpk"
ALGOS="LPBOUND RRNZ METAGREEDY METAGREEDYLIGHT METAHVP VP_PP_DSUM_W1 HVP_PP_DSUM_ASUM_W1 HVP_PP_DSUM_ASUM_W1_R_E"

LOCKDIR=${BASEDIR}/locks
PROBDIR=/home/mstillwell/work/flexsched_simulator/probsets/large_synth_probs-hmem

for T in {1..8}; do
(
for PROBFILEPATH in ${PROBDIR}/prob-64-250-*.prob; do
    PROBFILE=$(basename ${PROBFILEPATH})
    exec 200>${LOCKDIR}/${PROBFILE}.lck
    flock -xn 200
    if [ 0 = $? ]; then
        RESULTFILEPATH=${PROBDIR}/$(basename ${PROBFILE} .prob).result
        LD_LIBRARY_PATH=${HOME}/progs/base/lib ${FLEXSCHED} "${ALGOS}" ${PROBFILEPATH} ${RESULTFILEPATH}
        flock -u 200
    #else echo "Failed to get lock for ${PROBFILE}, trying next..."
    fi
    exec 200>&-
done
) &
done
wait
