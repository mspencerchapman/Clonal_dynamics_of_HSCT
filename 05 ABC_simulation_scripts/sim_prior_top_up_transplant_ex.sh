#!/bin/bash
while getopts a:b:c:d:e:f:g:h:i:j: flag
do
    case "${flag}" in
        a) SIM_RUN_INDEX=${OPTARG};;
        b) SOURCE_PATH=${OPTARG};;
        c) LOGS_PATH=${OPTARG};;
        d) PAIR_ID=${OPTARG};;
        e) MODEL_NAME=${OPTARG};;
        f) SIM_PATH=${OPTARG};;
        g) SIM_TREES_PATH=${OPTARG};;
        h) PRIOR_SAMPLE_DRIFT_FILE=${OPTARG};;
        i) PRIOR_SAMPLE_DRIVERS_FILE=${OPTARG};;
        j) N_SIMS_PER_RUN=${OPTARG};;
    esac
done
echo "SIM_RUN_INDEX: ${SIM_RUN_INDEX}";
echo "SOURCE_PATH: ${SOURCE_PATH}";
echo "LOGS_PATH: ${LOGS_PATH}";
echo "PAIR_ID: ${PAIR_ID}";
echo "MODEL_NAME: ${MODEL_NAME}";
echo "SIM_PATH: ${SIM_PATH}";
echo "SIM_TREES_PATH: ${SIM_TREES_PATH}";
echo "PRIOR_SAMPLE_DRIFT_FILE: ${PRIOR_SAMPLE_DRIFT_FILE}";
echo "PRIOR_SAMPLE_DRIVERS_FILE: ${PRIOR_SAMPLE_DRIVERS_FILE}";
echo "N_SIMS_PER_RUN: ${N_SIMS_PER_RUN}";

N_RUNS_ACCUMULATED=0

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

# SIM_RUN_INDEX=${JOB_INDEX}
echo "SIM_RUN_INDEX = ${SIM_RUN_INDEX}"

JOB_INDEX=${SIM_RUN_INDEX}

# RAND_SEED=${JOB_INDEX} # Not used!

R_PARAMS="${SOURCE_PATH} ${SIM_PATH} ${SIM_TREES_PATH} ${PRIOR_SAMPLE_DRIFT_FILE} ${PRIOR_SAMPLE_DRIVERS_FILE} ${PAIR_ID} ${MODEL_NAME} ${N_SIMS_PER_RUN} ${SIM_RUN_INDEX} ${N_RUNS_ACCUMULATED}"

echo "R_PARAMS = ${R_PARAMS}"

R_script_name="sim_prior_transplant.R"
Rout_filename="sim_prior_transplant"

module load R/4.1.0
# export R_LIBS="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1"
# export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.2.7"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.2.7:/nfs/users/nfs_m/ms56/R/x86_64-pc-linux-gnu-library/4.1"


CMD="R CMD BATCH"

cd ${SOURCE_PATH}

${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_PAIR_${PAIR_ID}_JOB_${JOB_ID}_RUN_INDEX_${SIM_RUN_INDEX}.Rout"


