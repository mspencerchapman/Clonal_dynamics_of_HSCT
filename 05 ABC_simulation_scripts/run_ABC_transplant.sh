#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
SIM_PATH=${3}
ABC_OUTPUT_PATH=${4}
POST_SAMPLE_FILE_NAME=${5}
OBS_STATS_PATH=${6}
OBS_STATS_FILE_NAME=${7}
POOLED_SIM_TREES_PATH=${8}
N_SIMS_ACCEPT=${9}
N_SIMS_MAX=${10}
MODEL_NAME=${11}
STAT_SET_NAME=${12}
ABC_METHOD=${13}
PAIR_NAME=${14}

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

R_PARAMS="${SOURCE_PATH} ${SIM_PATH} ${ABC_OUTPUT_PATH} ${POST_SAMPLE_FILE_NAME} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${POOLED_SIM_TREES_PATH} ${N_SIMS_ACCEPT} ${N_SIMS_MAX} ${MODEL_NAME} ${STAT_SET_NAME} ${ABC_METHOD} ${PAIR_NAME}"

echo "R_PARAMS = ${R_PARAMS}"


R_script_name="run_ABC_transplant.R"
Rout_filename="run_ABC_transplant"


module load R/4.1.0
# export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop"
# export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.3.0"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.3.0:/nfs/users/nfs_m/ms56/R/x86_64-pc-linux-gnu-library/4.1"


CMD="R CMD BATCH"

cd ${SOURCE_PATH}

${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_${ABC_METHOD}_PATIENT_${PAIR_NAME}_JOB_${JOB_ID}.Rout"

