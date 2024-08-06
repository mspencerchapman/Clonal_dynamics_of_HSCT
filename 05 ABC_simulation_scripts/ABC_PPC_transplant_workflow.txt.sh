
####################################################################################################
####################################################################################################
#
# "Round 1" simulations are the simulations used to perform ABC (Approximate Bayesian Computation) analysis.
# Run script: sim_prior_transplant.nf;
#
####################################################################################################
# 
# Do ABC (Approximate Bayesian Computation) analysis;
# Run script: run_ABC_transplant.sh;
#
####################################################################################################
#
# "Round 2" simulations ("resimulations") are the simulations used to perform PPC analysis.
# Run script: resim_conditionals_transplant.nf;
#
####################################################################################################
# 
# Do PPC (Posterior Predictive Model Checks) analysis;
# Run script: compute_ppc_results_seq_transplant.sh;
#
####################################################################################################
####################################################################################################

####################################################################################################
#
# "Round 1" simulations are the simulations used to perform ABC (Approximate Bayesian Computation) analysis.
# Run script: sim_prior_transplant.nf;
#
####################################################################################################

# load software modules;
LSB_DEFAULTGROUP=team154-grp
module load nextflow

# Specify model;

MODEL_NAME="m1" 
# MODEL_NAME="m2"
# MODEL_NAME="m3"

# Set-up work directory;
farm_path="/lustre/scratch125/casm/team273jn/kd7/farm_work/"
cd ${farm_path}

project_path=${farm_path}"ABC/PPC_for_transplants/rerun_ABC_PPC/ABC_PPC/array_jobs_3/"
cd ${farm_path}
mkdir -p ABC/PPC_for_transplants/rerun_ABC_PPC/ABC_PPC/array_jobs_3/

SOURCE_PATH=${project_path}"source/"

project_out_dir_name="output_MODEL_"${MODEL_NAME}
project_out_path=${project_path}${project_out_dir_name}"/"
cd ${project_path}
mkdir ${project_out_dir_name}

####################################################################################################
# Specify patient (donor-recipient) pair array;

PAIR_ID_ARRAY=(11 13 21 24 25 28 31 38 40 41)
PAIR_NAME_ARRAY=("Pair11" "Pair13" "Pair21" "Pair24" "Pair25" "Pair28" "Pair31" "Pair38" "Pair40" "Pair41")

N_PAIRS=${#PAIR_NAME_ARRAY[@]}

# Specify priors;

# Specify prior on HSC population size. 
# Use posterior sample generated from two young individuals (KX001, KX002) from E. Mitchell et al, 2022, using sequential ABC.
# This is the file: data/reference_files/HSC_population_posterior_sample.txt within the GitHub repo.
PRIOR_SAMPLE_DRIFT_FILE="/lustre/scratch126/casm/team273jn/kd7/farm_work/ABC/PPC_for_transplants/rerun_ABC_PPC/seq_ABC_KX001_KX002/seq_ABC_1_2/array_jobs_6/output/output_KX002/abc_output_round_1/abc_ridge_output/HSC_population_posterior_sample.txt"

# Specify prior on "driver" model parameters;
# Use posterior sample generated from 4 oldest individuals from E. Mitchell et al, 2022, using sequential ABC.
# This is the file: data/reference_files/driver_parameter_posterior_sample.txt within the GitHub repo.
PRIOR_SAMPLE_DRIVERS_FILE="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/driver_parameter_posterior_sample.txt"

####################################################################################################

N_SIMS_MAX_1=100000 # N_SIMS_MAX_1=100,000
N_SIMS_ACCEPT_REJECTION_1=1000
N_SIMS_ACCEPT_REGRESSION_1=2000

N_SIMS_PER_RUN=100
N_RUNS_1=1000

N_SIMS_1=$((${N_SIMS_PER_RUN}*${N_RUNS_1}))

N_RUNS_ACCUMULATED=0

####################################################################################################
# Run for-loop:

POST_SAMPLE_FILE_NAME="posterior_sample.txt"

for i in "${!PAIR_NAME_ARRAY[@]}"
do
 echo "i = $i"
 
 PAIR_NAME=${PAIR_NAME_ARRAY[$i]}
 PAIR_ID=${PAIR_ID_ARRAY[$i]}
 echo "PAIR_NAME = ${PAIR_NAME}"
 echo "PAIR_ID = ${PAIR_ID}"
 
 ################################################################################################################################
 # create output directories;
 
 PATIENT_OUT_PATH=${project_out_path}"output_"${PAIR_NAME}"/"
 cd ${project_out_path}
 mkdir "output_"${PAIR_NAME}
 
 LOGS_PATH=${PATIENT_OUT_PATH}"logs/"
 LOGS_SIM_PATH=${PATIENT_OUT_PATH}"logs_sim/"
 PRIOR_SIM_TREES_PATH=${PATIENT_OUT_PATH}"sim_trees_round_1/"
 POOLED_SIM_TREES_PATH=${PATIENT_OUT_PATH}"pooled_sim_trees_round_1/"
 PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"sims_round_1/"
 POOLED_PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"pooled_sims_round_1/"
 
 cd ${PATIENT_OUT_PATH}
 mkdir logs
 mkdir logs_sim
 mkdir sim_trees_round_1
 mkdir pooled_sim_trees_round_1
 mkdir sims_round_1
 mkdir pooled_sims_round_1
 
 cd ${PRIOR_SIM_PATH}
 mkdir original
 mkdir pre_interval_divided
 mkdir peri_interval_divided
 mkdir peri_interval_narrower
 mkdir peri_interval_wider
 
 ################################################################################################################################
 # Run script: sim_prior_transplant.nf;
 
 JOB_NAME_SIM="sim_prior_1_"${PAIR_NAME}
 
 # NFLOW_SCRIPT_SIM="sim_prior_transplant.nf"
 NFLOW_SCRIPT_SIM="sim_prior_top_up_transplant.nf"
 STDOUT_NAME_SIM="sim_prior_transplant"
 MEM_SIM=5000 # 5,000 MB
 QUEUE_SIM=week
 
 cd ${SOURCE_PATH}
 chmod +x *.sh
 
 bsub -G team154-grp -J "${JOB_NAME_SIM}" -q "${QUEUE_SIM}" -n 1 -M ${MEM_SIM} -R "select[mem>${MEM_SIM}] rusage[mem=${MEM_SIM}]" -o ${LOGS_PATH}${STDOUT_NAME_SIM}_%J.stdout -e ${LOGS_PATH}${STDOUT_NAME_SIM}_%J.stderr "${SOURCE_PATH}${NFLOW_SCRIPT_SIM}" --SOURCE_PATH "${SOURCE_PATH}" --project_out_path "${project_out_path}" --PAIR_ID "${PAIR_ID}" --N_SIMS_PER_RUN "${N_SIMS_PER_RUN}" --N_RUNS "${N_RUNS_1}" --MODEL_NAME "${MODEL_NAME}" --PRIOR_SAMPLE_DRIFT_FILE "${PRIOR_SAMPLE_DRIFT_FILE}" --PRIOR_SAMPLE_DRIVERS_FILE "${PRIOR_SAMPLE_DRIVERS_FILE}" --N_RUNS_ACCUMULATED "${N_RUNS_ACCUMULATED}"
 
 ################################################################################################################################
 
done

################################################################################################################################

################################################################################################################################
# 
# Do ABC (Approximate Bayesian Computation) analysis;
# Run script: run_ABC_transplant.sh;
# Run the following script ONLY AFTER (round 1) sims (for ALL donors) have completed.
#
################################################################################################################################

# load software modules;
LSB_DEFAULTGROUP=team154-grp
# module load nextflow

# Specify model;

MODEL_NAME="m1" 
# MODEL_NAME="m2"
# MODEL_NAME="m3"

# Specify summary stats sets;
STAT_SET_NAME_ARRAY=("original" "pre_interval_divided" "peri_interval_divided" "peri_interval_narrower" "peri_interval_wider")

# Specify work directory;
farm_path="/lustre/scratch125/casm/team273jn/kd7/farm_work/"
project_path=${farm_path}"ABC/PPC_for_transplants/rerun_ABC_PPC/ABC_PPC/array_jobs_3/"
SOURCE_PATH=${project_path}"source/"
project_out_dir_name="output_MODEL_"${MODEL_NAME}
project_out_path=${project_path}${project_out_dir_name}"/"

####################################################################################################
# Specify patient (donor-recipient) pair array;

PAIR_ID_ARRAY=(11 13 21 24 25 28 31 38 40 41)
PAIR_NAME_ARRAY=("Pair11" "Pair13" "Pair21" "Pair24" "Pair25" "Pair28" "Pair31" "Pair38" "Pair40" "Pair41")

N_PAIRS=${#PAIR_NAME_ARRAY[@]}

####################################################################################################

N_SIMS_MAX_1=100000 # N_SIMS_MAX_1=100,000
N_SIMS_ACCEPT_REJECTION_1=1000
N_SIMS_ACCEPT_REGRESSION_1=2000

####################################################################################################
# Run for-loop:

POST_SAMPLE_FILE_NAME="posterior_sample.txt"

for i in "${!PAIR_NAME_ARRAY[@]}"
do
 echo "i = $i"
 
 PAIR_NAME=${PAIR_NAME_ARRAY[$i]}
 PAIR_ID=${PAIR_ID_ARRAY[$i]}
 echo "PAIR_NAME = ${PAIR_NAME}"
 echo "PAIR_ID = ${PAIR_ID}"
 
 ################################################################################################################################
 # create output directories;
 
 PATIENT_OUT_PATH=${project_out_path}"output_"${PAIR_NAME}"/"
 cd ${project_out_path}
 mkdir "output_"${PAIR_NAME}
 
 LOGS_PATH=${PATIENT_OUT_PATH}"logs/"
 LOGS_SIM_PATH=${PATIENT_OUT_PATH}"logs_sim/"
 PRIOR_SIM_TREES_PATH=${PATIENT_OUT_PATH}"sim_trees_round_1/"
 POOLED_SIM_TREES_PATH=${PATIENT_OUT_PATH}"pooled_sim_trees_round_1/"
 PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"sims_round_1/"
 POOLED_PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"pooled_sims_round_1/"
 
 ################################################################################################################################
 
 for STAT_SET_NAME in ${STAT_SET_NAME_ARRAY[*]}
 do
  echo "STAT_SET_NAME = ${STAT_SET_NAME}"
  
  ################################################################################################################################
  # create output directories;
  
  PRIOR_SIM_STAT_SET_PATH=${PRIOR_SIM_PATH}${STAT_SET_NAME}"/"
  
  PATIENT_STAT_SET_DIR_NAME="post_output_STATS_"${STAT_SET_NAME}
  PATIENT_STAT_SET_OUT_PATH=${PATIENT_OUT_PATH}${PATIENT_STAT_SET_DIR_NAME}"/"
  cd ${PATIENT_OUT_PATH}
  mkdir ${PATIENT_STAT_SET_DIR_NAME}
  
  LOGS_ABC_PATH=${PATIENT_STAT_SET_OUT_PATH}"logs_abc/"
  OBS_STATS_PATH=${PATIENT_STAT_SET_OUT_PATH}"obs_stats_round_1/"
  ABC_OUTPUT_PATH=${PATIENT_STAT_SET_OUT_PATH}"abc_output_round_1/"
  
  cd ${PATIENT_STAT_SET_OUT_PATH}
  mkdir logs_abc
  mkdir obs_stats_round_1
  mkdir abc_output_round_1
  
  ABC_REJECTION_OUTPUT_PATH=${ABC_OUTPUT_PATH}"abc_rejection_output/"
  ABC_RIDGE_OUTPUT_PATH=${ABC_OUTPUT_PATH}"abc_ridge_output/"
  ABC_LOCAL_OUTPUT_PATH=${ABC_OUTPUT_PATH}"abc_local_output/"
  ABC_NN_OUTPUT_PATH=${ABC_OUTPUT_PATH}"abc_nn_output/"
  
  cd ${ABC_OUTPUT_PATH}
  mkdir abc_rejection_output
  mkdir abc_ridge_output
  mkdir abc_local_output
  mkdir abc_nn_output
  
  ################################################################################################################################
  # Apply ABC "rejection" method to generate sample from (approximate) posterior distribution;
  # Run script: run_ABC_transplant.sh;
  
  ABC_METHOD="rejection" # ABC rejection method;
  
  N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REJECTION_1}
  # N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REGRESSION_1}
  
  PARAMS_ABC="${SOURCE_PATH} ${LOGS_ABC_PATH} ${PRIOR_SIM_STAT_SET_PATH} ${ABC_REJECTION_OUTPUT_PATH} ${POST_SAMPLE_FILE_NAME} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${POOLED_SIM_TREES_PATH} ${N_SIMS_ACCEPT_1} ${N_SIMS_MAX_1} ${MODEL_NAME} ${STAT_SET_NAME} ${ABC_METHOD} ${PAIR_NAME}"
  
  JOB_NAME_ABC="ABC_1_rejection_"${PAIR_NAME}"_"${STAT_SET_NAME}
  
  SH_SCRIPT_ABC="run_ABC_transplant.sh"
  STDOUT_NAME_ABC="run_ABC_transplant_1_rejection"
  MEM_ABC=5000
  QUEUE_ABC=normal
  
  cd ${SOURCE_PATH}
  chmod +x *.sh
  
  bsub -G team154-grp -J ${JOB_NAME_ABC} -q "${QUEUE_ABC}" -M ${MEM_ABC} -R "select[mem>${MEM_ABC}] rusage[mem=${MEM_ABC}] span[hosts=1]" -o ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stdout -e ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stderr "${SOURCE_PATH}${SH_SCRIPT_ABC} ${PARAMS_ABC}"
  
  ################################################################################################################################
  # Apply ABC "ridge regression" method to generate sample from (approximate) posterior distribution;
  # Run script: run_ABC_transplant.sh;
  
  ABC_METHOD="ridge" # ABC ridge regression method;
  
  # N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REJECTION_1}
  N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REGRESSION_1}
  
  PARAMS_ABC="${SOURCE_PATH} ${LOGS_ABC_PATH} ${PRIOR_SIM_STAT_SET_PATH} ${ABC_REJECTION_OUTPUT_PATH} ${POST_SAMPLE_FILE_NAME} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${POOLED_SIM_TREES_PATH} ${N_SIMS_ACCEPT_1} ${N_SIMS_MAX_1} ${MODEL_NAME} ${STAT_SET_NAME} ${ABC_METHOD} ${PAIR_NAME}"
  
  JOB_NAME_ABC="ABC_1_ridge_"${PAIR_NAME}"_"${STAT_SET_NAME}
  
  SH_SCRIPT_ABC="run_ABC_transplant.sh"
  STDOUT_NAME_ABC="run_ABC_transplant_1_ridge"
  MEM_ABC=5000
  QUEUE_ABC=normal
  
  cd ${SOURCE_PATH}
  chmod +x *.sh
  
  bsub -G team154-grp -J ${JOB_NAME_ABC} -q "${QUEUE_ABC}" -M ${MEM_ABC} -R "select[mem>${MEM_ABC}] rusage[mem=${MEM_ABC}] span[hosts=1]" -o ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stdout -e ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stderr "${SOURCE_PATH}${SH_SCRIPT_ABC} ${PARAMS_ABC}"
  
  ################################################################################################################################
  # Apply ABC "local linear regression" method to generate sample from (approximate) posterior distribution;
  # Run script: run_ABC_transplant.sh;
  
  ABC_METHOD="loclinear" # ABC local linear regression method;
  
  # N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REJECTION_1}
  N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REGRESSION_1}
  
  PARAMS_ABC="${SOURCE_PATH} ${LOGS_ABC_PATH} ${PRIOR_SIM_STAT_SET_PATH} ${ABC_REJECTION_OUTPUT_PATH} ${POST_SAMPLE_FILE_NAME} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${POOLED_SIM_TREES_PATH} ${N_SIMS_ACCEPT_1} ${N_SIMS_MAX_1} ${MODEL_NAME} ${STAT_SET_NAME} ${ABC_METHOD} ${PAIR_NAME}"
  
  JOB_NAME_ABC="ABC_1_local_"${PAIR_NAME}"_"${STAT_SET_NAME}
  
  SH_SCRIPT_ABC="run_ABC_transplant.sh"
  STDOUT_NAME_ABC="run_ABC_transplant_1_local"
  MEM_ABC=5000
  QUEUE_ABC=normal
  
  cd ${SOURCE_PATH}
  chmod +x *.sh
  
  bsub -G team154-grp -J ${JOB_NAME_ABC} -q "${QUEUE_ABC}" -M ${MEM_ABC} -R "select[mem>${MEM_ABC}] rusage[mem=${MEM_ABC}] span[hosts=1]" -o ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stdout -e ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stderr "${SOURCE_PATH}${SH_SCRIPT_ABC} ${PARAMS_ABC}"
  
  ################################################################################################################################
  # Apply ABC "neural network regression" method to generate sample from (approximate) posterior distribution;
  # Run script: run_ABC_transplant.sh;
  
  ABC_METHOD="neuralnet" # ABC neural network regression method;
  
  # N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REJECTION_1}
  N_SIMS_ACCEPT_1=${N_SIMS_ACCEPT_REGRESSION_1}
  
  PARAMS_ABC="${SOURCE_PATH} ${LOGS_ABC_PATH} ${PRIOR_SIM_STAT_SET_PATH} ${ABC_REJECTION_OUTPUT_PATH} ${POST_SAMPLE_FILE_NAME} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${POOLED_SIM_TREES_PATH} ${N_SIMS_ACCEPT_1} ${N_SIMS_MAX_1} ${MODEL_NAME} ${STAT_SET_NAME} ${ABC_METHOD} ${PAIR_NAME}"
  
  JOB_NAME_ABC="ABC_1_nn_"${PAIR_NAME}"_"${STAT_SET_NAME}
  
  SH_SCRIPT_ABC="run_ABC_transplant.sh"
  STDOUT_NAME_ABC="run_ABC_transplant_1_nn"
  MEM_ABC=5000
  QUEUE_ABC=normal
  
  cd ${SOURCE_PATH}
  chmod +x *.sh
  
  bsub -G team154-grp -J ${JOB_NAME_ABC} -q "${QUEUE_ABC}" -M ${MEM_ABC} -R "select[mem>${MEM_ABC}] rusage[mem=${MEM_ABC}] span[hosts=1]" -o ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stdout -e ${LOGS_PATH}${STDOUT_NAME_ABC}_%J.stderr "${SOURCE_PATH}${SH_SCRIPT_ABC} ${PARAMS_ABC}"
  
  ################################################################################################################################
  
  done # for STAT_SET_NAME in ${STAT_SET_NAME_ARRAY[*]}
  
 ################################################################################################################################
 
done # for i in "${!PAIR_NAME_ARRAY[@]}"

####################################################################################################

####################################################################################################
#
# "Round 2" simulations ("resimulations") are the simulations used to perform PPC analysis.
# Run script: resim_conditionals_transplant.nf;
# Run the following script ONLY AFTER the ABC jobs have completed.
#
####################################################################################################

# load software modules;
LSB_DEFAULTGROUP=team154-grp
module load nextflow

# Specify model;

MODEL_NAME="m1" 
# MODEL_NAME="m2"
# MODEL_NAME="m3"

# Specify ABC summary stats set;
ABC_STAT_SET_NAME="original"

# Specify ABC method;

ABC_METHOD_NAME_ARRAY=("rejection" "ridge" "loclinear" "neuralnet")
ABC_METHOD_DIR_ARRAY=("abc_rejection_output" "abc_ridge_output" "abc_local_output" "abc_nn_output")

ABC_METHOD_INDEX=0 # "rejection"
# ABC_METHOD_INDEX=1 # "ridge"
# ABC_METHOD_INDEX=2 # "loclinear"
# ABC_METHOD_INDEX=3 # "neuralnet"

ABC_METHOD=${ABC_METHOD_NAME_ARRAY[${ABC_METHOD_INDEX}]}
POST_SAMPLE_DIR_NAME=${ABC_METHOD_DIR_ARRAY[${ABC_METHOD_INDEX}]}
POST_SAMPLE_FILE_NAME="posterior_sample.txt"


# Specify work directory;
farm_path="/lustre/scratch125/casm/team273jn/kd7/farm_work/"
project_path=${farm_path}"ABC/PPC_for_transplants/rerun_ABC_PPC/ABC_PPC/array_jobs_3/"
SOURCE_PATH=${project_path}"source/"
project_out_dir_name="output_MODEL_"${MODEL_NAME}
project_out_path=${project_path}${project_out_dir_name}"/"

####################################################################################################
# Specify patient (donor-recipient) pair array;

PAIR_ID_ARRAY=(11 13 21 24 25 28 31 38 40 41)
PAIR_NAME_ARRAY=("Pair11" "Pair13" "Pair21" "Pair24" "Pair25" "Pair28" "Pair31" "Pair38" "Pair40" "Pair41")

N_PAIRS=${#PAIR_NAME_ARRAY[@]}

################################################################################################################################

BEGIN_POST_OBS=1
END_POST_OBS=1000

N_POST_OBS=$((${END_POST_OBS}-${BEGIN_POST_OBS}))
N_POST_OBS=$((${N_POST_OBS}+1))

# N_SIMS_PER_RUN=100
N_RESIMS_PER_RUN=100
N_RUNS_PER_POST_OBS=10
N_RUNS_ACCUMULATED_PER_POST_OBS=0

N_RESIM_RUNS=$((${N_RUNS_PER_POST_OBS}*${N_POST_OBS}))
N_RESIMS_TOTAL=$((${N_RESIMS_PER_RUN}*${N_RESIM_RUNS}))
N_RESIMS_PER_POST_OBS=$((${N_RESIMS_PER_RUN}*${N_RUNS_PER_POST_OBS}))

echo "N_POST_OBS = ${N_POST_OBS}"
echo "N_RESIMS_PER_RUN = ${N_RESIMS_PER_RUN}"
echo "N_RUNS_PER_POST_OBS = ${N_RUNS_PER_POST_OBS}"
echo "N_RESIMS_PER_POST_OBS = ${N_RESIMS_PER_POST_OBS}"

####################################################################################################
# Run for-loop:

POST_SAMPLE_FILE_NAME="posterior_sample.txt"

for i in "${!PAIR_NAME_ARRAY[@]}"
do
 echo "i = $i"
 
 PAIR_NAME=${PAIR_NAME_ARRAY[$i]}
 PAIR_ID=${PAIR_ID_ARRAY[$i]}
 echo "PAIR_NAME = ${PAIR_NAME}"
 echo "PAIR_ID = ${PAIR_ID}"
 
 ################################################################################################################################
 # create output directories;
 
 PATIENT_OUT_PATH=${project_out_path}"output_"${PAIR_NAME}"/"
 
 LOGS_PATH=${PATIENT_OUT_PATH}"logs/"
 LOGS_SIM_PATH=${PATIENT_OUT_PATH}"logs_sim/"
 LOGS_CONDITIONAL_RESIM_PATH=${PATIENT_OUT_PATH}"logs_resim/"
 
 PRIOR_SIM_TREES_PATH=${PATIENT_OUT_PATH}"sim_trees_round_1/"
 POOLED_SIM_TREES_PATH=${PATIENT_OUT_PATH}"pooled_sim_trees_round_1/"
 PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"sims_round_1/"
 POOLED_PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"pooled_sims_round_1/"
 
 CONDITIONAL_RESIM_TREES_PATH=${PATIENT_OUT_PATH}"conditional_resims_trees/"
 POOLED_CONDITIONAL_RESIM_TREES_PATH=${PATIENT_OUT_PATH}"pooled_conditional_resims_trees/"
 CONDITIONAL_RESIM_PATH=${PATIENT_OUT_PATH}"conditional_resims/"
 POOLED_CONDITIONAL_RESIM_PATH=${PATIENT_OUT_PATH}"pooled_conditional_resims/"
 PPC_OUTPUT_PATH=${PATIENT_OUT_PATH}"ppc_results/"
 
 cd ${PATIENT_OUT_PATH}
 mkdir logs_resim
 mkdir conditional_resims_trees
 mkdir pooled_conditional_resims_trees
 mkdir conditional_resims
 mkdir pooled_conditional_resims
 mkdir ppc_results
 
 cd ${CONDITIONAL_RESIM_PATH}
 mkdir original
 mkdir pre_interval_divided
 mkdir peri_interval_divided
 mkdir peri_interval_narrower
 mkdir peri_interval_wider
 
 PATIENT_ABC_STAT_SET_DIR_NAME="post_output_STATS_"${ABC_STAT_SET_NAME}
 PATIENT_ABC_STAT_SET_OUT_PATH=${PATIENT_OUT_PATH}${PATIENT_ABC_STAT_SET_DIR_NAME}"/"
 cd ${PATIENT_OUT_PATH}
 mkdir ${PATIENT_ABC_STAT_SET_DIR_NAME}
 
 OBS_STATS_PATH=${PATIENT_ABC_STAT_SET_OUT_PATH}"obs_stats_round_1/"
 OBS_STATS_FILE_NAME="summary_stats.txt"
 
 ABC_OUTPUT_PATH=${PATIENT_ABC_STAT_SET_OUT_PATH}"abc_output_round_1/"
 
 POST_SAMPLE_PATH=${ABC_OUTPUT_PATH}${POST_SAMPLE_DIR_NAME}"/"
 POST_SAMPLE_FILE=${POST_SAMPLE_PATH}${POST_SAMPLE_FILE_NAME}
 
 ################################################################################################################################
 # Run script: resim_conditionals_transplant.nf;
 
 JOB_NAME_RESIM_CONDITIONALS="resim_conditionals_"${MODEL_NAME}"_"${PAIR_NAME}
 
 NFLOW_SCRIPT_RESIM_CONDITIONALS="resim_conditionals_transplant.nf"
 STDOUT_NAME_RESIM_CONDITIONALS="resim_conditionals"
 MEM_RESIM_CONDITIONALS=5000 # 5,000 MB
 QUEUE_RESIM_CONDITIONALS=week
 
 cd ${SOURCE_PATH}
 chmod +x *.sh
 
 bsub -G team154-grp -J "${JOB_NAME_RESIM_CONDITIONALS}" -q "${QUEUE_RESIM_CONDITIONALS}" -M ${MEM_RESIM_CONDITIONALS} -R "select[mem>${MEM_RESIM_CONDITIONALS}] rusage[mem=${MEM_RESIM_CONDITIONALS}]" -o ${LOGS_PATH}${STDOUT_NAME_RESIM_CONDITIONALS}_%J.stdout -e ${LOGS_PATH}${STDOUT_NAME_RESIM_CONDITIONALS}_%J.stderr "${SOURCE_PATH}${NFLOW_SCRIPT_RESIM_CONDITIONALS}" --SOURCE_PATH "${SOURCE_PATH}" --project_out_path "${project_out_path}" --PAIR_ID "${PAIR_ID}" --MODEL_NAME "${MODEL_NAME}" --POST_SAMPLE_FILE "${POST_SAMPLE_FILE}" --BEGIN_POST_OBS "${BEGIN_POST_OBS}" --END_POST_OBS "${END_POST_OBS}" --N_RUNS_PER_POST_OBS "${N_RUNS_PER_POST_OBS}" --N_RUNS_ACCUMULATED_PER_POST_OBS "${N_RUNS_ACCUMULATED_PER_POST_OBS}" --N_RESIMS_PER_RUN "${N_RESIMS_PER_RUN}" --N_RESIM_RUNS "${N_RESIM_RUNS}"
 
 ################################################################################################################################
 
done

################################################################################################################################

################################################################################################################################
# 
# Do PPC (Posterior Predictive Model Checks) analysis;
# Run script: compute_conditionals_drivers.sh;
# Run the following loop ONLY AFTER (round 2) resims (for ALL donors) have completed.
#
################################################################################################################################

# load software modules;
LSB_DEFAULTGROUP=team154-grp
# module load nextflow

# Specify model;

MODEL_NAME="m1" 
# MODEL_NAME="m2"
# MODEL_NAME="m3"

# Specify ABC summary stats set;
ABC_STAT_SET_NAME="original"

# Specify ABC method;

ABC_METHOD_NAME_ARRAY=("rejection" "ridge" "loclinear" "neuralnet")
ABC_METHOD_DIR_ARRAY=("abc_rejection_output" "abc_ridge_output" "abc_local_output" "abc_nn_output")

ABC_METHOD_INDEX=0 # "rejection"
# ABC_METHOD_INDEX=1 # "ridge"
# ABC_METHOD_INDEX=2 # "loclinear"
# ABC_METHOD_INDEX=3 # "neuralnet"

ABC_METHOD=${ABC_METHOD_NAME_ARRAY[${ABC_METHOD_INDEX}]}
POST_SAMPLE_DIR_NAME=${ABC_METHOD_DIR_ARRAY[${ABC_METHOD_INDEX}]}
POST_SAMPLE_FILE_NAME="posterior_sample.txt"


# Specify work directory;
farm_path="/lustre/scratch125/casm/team273jn/kd7/farm_work/"
project_path=${farm_path}"ABC/PPC_for_transplants/rerun_ABC_PPC/ABC_PPC/array_jobs_3/"
SOURCE_PATH=${project_path}"source/"
project_out_dir_name="output_MODEL_"${MODEL_NAME}
project_out_path=${project_path}${project_out_dir_name}"/"

####################################################################################################
# Specify patient (donor-recipient) pair array;

PAIR_ID_ARRAY=(11 13 21 24 25 28 31 38 40 41)
PAIR_NAME_ARRAY=("Pair11" "Pair13" "Pair21" "Pair24" "Pair25" "Pair28" "Pair31" "Pair38" "Pair40" "Pair41")

N_PAIRS=${#PAIR_NAME_ARRAY[@]}

################################################################################################################################

BEGIN_POST_OBS=1
END_POST_OBS=1000

N_POST_OBS=$((${END_POST_OBS}-${BEGIN_POST_OBS}))
N_POST_OBS=$((${N_POST_OBS}+1))

# N_SIMS_PER_RUN=100
N_RESIMS_PER_RUN=100
N_RUNS_PER_POST_OBS=10
N_RUNS_ACCUMULATED_PER_POST_OBS=0

N_RESIM_RUNS=$((${N_RUNS_PER_POST_OBS}*${N_POST_OBS}))
N_RESIMS_TOTAL=$((${N_RESIMS_PER_RUN}*${N_RESIM_RUNS}))
N_RESIMS_PER_POST_OBS=$((${N_RESIMS_PER_RUN}*${N_RUNS_PER_POST_OBS}))

echo "N_POST_OBS = ${N_POST_OBS}"
echo "N_RESIMS_PER_RUN = ${N_RESIMS_PER_RUN}"
echo "N_RUNS_PER_POST_OBS = ${N_RUNS_PER_POST_OBS}"
echo "N_RESIMS_PER_POST_OBS = ${N_RESIMS_PER_POST_OBS}"

####################################################################################################

POST_SAMPLE_FILE_NAME="posterior_sample.txt"

OBS_STATS_DIR_NAME="obs_stats_round_1"
OBS_STATS_FILE_NAME="summary_stats.txt"

POOLED_CONDITIONAL_RESIM_FILE_NAME="pooled_conditional_resims.txt"

####################################################################################################
# Specify minimum number of resims required for each "observation"
# (from the approximate posterior distribution)
# to be used in the posterior p-value calculation;
#

# MIN_RESIMS_PER_POST_OBS=1000
MIN_RESIMS_PER_POST_OBS=500

####################################################################################################
# Run for-loop:

POST_SAMPLE_FILE_NAME="posterior_sample.txt"

for i in "${!PAIR_NAME_ARRAY[@]}"
do
 echo "i = $i"
 
 PAIR_NAME=${PAIR_NAME_ARRAY[$i]}
 PAIR_ID=${PAIR_ID_ARRAY[$i]}
 echo "PAIR_NAME = ${PAIR_NAME}"
 echo "PAIR_ID = ${PAIR_ID}"
 
 ################################################################################################################################
 # create output directories;
 
 PATIENT_OUT_PATH=${project_out_path}"output_"${PAIR_NAME}"/"
 
 LOGS_PATH=${PATIENT_OUT_PATH}"logs/"
 LOGS_SIM_PATH=${PATIENT_OUT_PATH}"logs_sim/"
 LOGS_CONDITIONAL_RESIM_PATH=${PATIENT_OUT_PATH}"logs_resim/"
 
 PRIOR_SIM_TREES_PATH=${PATIENT_OUT_PATH}"sim_trees_round_1/"
 POOLED_SIM_TREES_PATH=${PATIENT_OUT_PATH}"pooled_sim_trees_round_1/"
 PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"sims_round_1/"
 POOLED_PRIOR_SIM_PATH=${PATIENT_OUT_PATH}"pooled_sims_round_1/"
 
 CONDITIONAL_RESIM_TREES_PATH=${PATIENT_OUT_PATH}"conditional_resims_trees/"
 POOLED_CONDITIONAL_RESIM_TREES_PATH=${PATIENT_OUT_PATH}"pooled_conditional_resims_trees/"
 CONDITIONAL_RESIM_PATH=${PATIENT_OUT_PATH}"conditional_resims/"
 POOLED_CONDITIONAL_RESIM_PATH=${PATIENT_OUT_PATH}"pooled_conditional_resims/"
 PPC_OUTPUT_PATH=${PATIENT_OUT_PATH}"ppc_results/"
 
 PATIENT_ABC_STAT_SET_DIR_NAME="post_output_STATS_"${ABC_STAT_SET_NAME}
 PATIENT_ABC_STAT_SET_OUT_PATH=${PATIENT_OUT_PATH}${PATIENT_ABC_STAT_SET_DIR_NAME}"/"
 cd ${PATIENT_OUT_PATH}
 mkdir ${PATIENT_ABC_STAT_SET_DIR_NAME}
 
 OBS_STATS_PATH=${PATIENT_ABC_STAT_SET_OUT_PATH}"obs_stats_round_1/"
 OBS_STATS_FILE_NAME="summary_stats.txt"
 
 ABC_OUTPUT_PATH=${PATIENT_ABC_STAT_SET_OUT_PATH}"abc_output_round_1/"
 
 POST_SAMPLE_PATH=${ABC_OUTPUT_PATH}${POST_SAMPLE_DIR_NAME}"/"
 POST_SAMPLE_FILE=${POST_SAMPLE_PATH}${POST_SAMPLE_FILE_NAME}
 
 ################################################################################################################################
 # Run script: compute_ppc_results_seq_transplant.sh;
 
 PARAMS_COMPUTE_PPC_RESULTS="${SOURCE_PATH} ${LOGS_PATH} ${CONDITIONAL_RESIM_ABC_STAT_SET_PATH} ${PPC_RESULTS_PATH} ${POST_SAMPLE_FILE} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${MIN_RESIMS_PER_POST_OBS} ${MODEL_NAME} ${ABC_STAT_SET_NAME} ${ABC_METHOD} ${PAIR_NAME}"
 
 JOB_NAME_COMPUTE_PPC_RESULTS="compute_ppc_results_seq_"${PAIR_NAME}
 
 SH_SCRIPT_COMPUTE_PPC_RESULTS="compute_ppc_results_seq_transplant.sh"
 STDOUT_NAME_COMPUTE_PPC_RESULTS="compute_ppc_results_seq"
 MEM_COMPUTE_PPC_RESULTS=5000 # 5,000 MB
 QUEUE_COMPUTE_PPC_RESULTS=normal
 
 cd ${SOURCE_PATH}
 chmod +x *.sh
 
 bsub -G team154-grp -J ${JOB_NAME_COMPUTE_PPC_RESULTS} -q "${QUEUE_COMPUTE_PPC_RESULTS}" -M ${MEM_COMPUTE_PPC_RESULTS} -R "select[mem>${MEM_COMPUTE_PPC_RESULTS}] rusage[mem=${MEM_COMPUTE_PPC_RESULTS}] span[hosts=1]" -o ${LOGS_PATH}${STDOUT_NAME_COMPUTE_PPC_RESULTS}_%J.stdout -e ${LOGS_PATH}${STDOUT_NAME_COMPUTE_PPC_RESULTS}_%J.stderr "${SOURCE_PATH}${SH_SCRIPT_COMPUTE_PPC_RESULTS} ${PARAMS_COMPUTE_PPC_RESULTS}"
 
 
 ################################################################################################################################
 
done # for i in "${!PAIR_NAME_ARRAY[@]}"

####################################################################################################
