#!/bin/bash
while getopts a:b:c:d:e:f:g:h:i:j:k:l:m: flag
do
    case "${flag}" in
        a) JOB_INDEX=${OPTARG};;
        b) SOURCE_PATH=${OPTARG};;
        c) LOGS_PATH=${OPTARG};;
        d) PAIR_ID=${OPTARG};;
        e) MODEL_NAME=${OPTARG};;
        f) RESIM_PATH=${OPTARG};;
        g) RESIM_TREES_PATH=${OPTARG};;
        h) POST_SAMPLE_FILE=${OPTARG};;
        i) BEGIN_POST_OBS=${OPTARG};;
        j) END_POST_OBS=${OPTARG};;
        k) N_RUNS_PER_POST_OBS=${OPTARG};;
        l) N_RUNS_ACCUMULATED_PER_POST_OBS=${OPTARG};;
        m) N_RESIMS_PER_RUN=${OPTARG};;
    esac
done
echo "JOB_INDEX: ${JOB_INDEX}";
echo "SOURCE_PATH: ${SOURCE_PATH}";
echo "LOGS_PATH: ${LOGS_PATH}";
echo "PAIR_ID: ${PAIR_ID}";
echo "MODEL_NAME: ${MODEL_NAME}";
echo "RESIM_PATH: ${RESIM_PATH}";
echo "RESIM_TREES_PATH: ${RESIM_TREES_PATH}";
echo "POST_SAMPLE_FILE: ${POST_SAMPLE_FILE}";
echo "BEGIN_POST_OBS: ${BEGIN_POST_OBS}";
echo "END_POST_OBS: ${END_POST_OBS}";
echo "N_RUNS_PER_POST_OBS: ${N_RUNS_PER_POST_OBS}";
echo "N_RUNS_ACCUMULATED_PER_POST_OBS: ${N_RUNS_ACCUMULATED_PER_POST_OBS}";
echo "N_RESIMS_PER_RUN: ${N_RESIMS_PER_RUN}";


# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}
echo "JOB_INDEX = ${JOB_INDEX}"
echo "JOB_ID = ${JOB_ID}"


TEMP_JOB_INDEX=$((${JOB_INDEX}-1))
# TEMP_POST_OBS_INDEX=$((${TEMP_JOB_INDEX}/${N_RUNS_PER_POST_OBS}))
TEMP_OBS_COUNT=$((${TEMP_JOB_INDEX}/${N_RUNS_PER_POST_OBS}))
TEMP_RESIM_RUN_INDEX=$((${TEMP_JOB_INDEX}%${N_RUNS_PER_POST_OBS}))

echo "TEMP_JOB_INDEX = ${TEMP_JOB_INDEX}"
echo "TEMP_OBS_COUNT = ${TEMP_OBS_COUNT}"
echo "TEMP_RESIM_RUN_INDEX = ${TEMP_RESIM_RUN_INDEX}"

# POST_OBS_INDEX=$((${TEMP_POST_OBS_INDEX}+1))
POST_OBS_INDEX=$((${TEMP_OBS_COUNT}+${BEGIN_POST_OBS}))
RESIM_RUN_INDEX=$((${TEMP_RESIM_RUN_INDEX}+1))
echo "POST_OBS_INDEX = ${POST_OBS_INDEX}"
echo "RESIM_RUN_INDEX = ${RESIM_RUN_INDEX}"

# RAND_SEED=${JOB_INDEX} # Not used!

R_PARAMS="${SOURCE_PATH} ${RESIM_PATH} ${RESIM_TREES_PATH} ${POST_SAMPLE_FILE} ${POST_OBS_INDEX} ${RESIM_RUN_INDEX} ${N_RESIMS_PER_RUN} ${N_RUNS_ACCUMULATED_PER_POST_OBS} ${PAIR_ID} ${MODEL_NAME}"

echo "R_PARAMS = ${R_PARAMS}"

R_script_name="resim_conditionals_transplant.R"
Rout_filename="resim_conditionals"

module load R/4.1.0
# export R_LIBS="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1"
# export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.2.7"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/4.1_rsimpop_2.2.7:/nfs/users/nfs_m/ms56/R/x86_64-pc-linux-gnu-library/4.1"


CMD="R CMD BATCH"

cd ${SOURCE_PATH}

${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_PAIR_${PAIR_ID}_JOB_${JOB_ID}_POST_OBS_${POST_OBS_INDEX}_RUN_INDEX_${RESIM_RUN_INDEX}.Rout"

