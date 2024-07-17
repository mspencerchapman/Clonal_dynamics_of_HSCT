#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
CONDITIONAL_RESIM_PATH=${3}
PPC_OUTPUT_PATH=${4}
POST_SAMPLE_FILE=${5}
OBS_STATS_PATH=${6}
OBS_STATS_FILE_NAME=${7}
MIN_RESIMS_PER_POST_OBS=${8}
MODEL_NAME=${9}
STAT_SET_NAME=${10}
ABC_METHOD=${11}
PAIR_NAME=${12}

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

R_PARAMS="${SOURCE_PATH} ${CONDITIONAL_RESIM_PATH} ${PPC_OUTPUT_PATH} ${POST_SAMPLE_FILE} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${MIN_RESIMS_PER_POST_OBS} ${MODEL_NAME} ${STAT_SET_NAME} ${ABC_METHOD} ${PAIR_NAME}"


R_script_name="compute_ppc_results_seq_transplant.R"
Rout_filename="compute_ppc_results_seq_transplant"

module load R/4.1.0
# export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop"
# export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.3.0"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.3.0:/nfs/users/nfs_m/ms56/R/x86_64-pc-linux-gnu-library/4.1"

CMD="R CMD BATCH"

cd ${SOURCE_PATH}

${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_${ABC_METHOD}_PATIENT_${PAIR_NAME}_JOB_${JOB_ID}.Rout"


