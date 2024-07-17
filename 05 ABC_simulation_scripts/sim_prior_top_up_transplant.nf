#!/usr/bin/env nextflow

nextflow.enable.dsl=2 // what does this do?

// inputs
params.SOURCE_PATH = "/source/" // default value (debugging only)
params.project_out_path = "/output/" // default value (debugging only)
params.PAIR_ID = 11 // default value (debugging only)
params.N_SIMS_PER_RUN = 2 // default value (debugging only)
params.N_RUNS_BEGIN = 10 // default value (debugging only)
params.N_RUNS_END = 20 // default value (debugging only)
params.MODEL_NAME = "m1" // default value (debugging only)
params.PRIOR_SAMPLE_DRIFT_FILE ="/lustre/scratch126/casm/team273jn/kd7/farm_work/ABC/PPC_for_transplants/rerun_ABC_PPC/seq_ABC_KX001_KX002/seq_ABC_1_2/array_jobs_6/output/output_KX002/abc_output_round_1/abc_ridge_output/posterior_sample.txt"
params.PRIOR_SAMPLE_DRIVERS_FILE = "/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/posterior_sample.txt"


process run_sims {

	input:
	val SIM_RUN_INDEX
	val SOURCE_PATH
	val LOGS_SIM_PATH
	val PAIR_ID
	val MODEL_NAME
	val SIM_PATH
	val SIM_TREES_PATH
	val PRIOR_SAMPLE_DRIFT_FILE
    val PRIOR_SAMPLE_DRIVERS_FILE
	val N_SIMS_PER_RUN
	
	output:
    stdout
    
    executor = 'lsf'
    // memory { 2.GB * task.attempt }
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 4 
    

    cpus 1
    // memory '4 GB'
    // clusterOptions = '-G team78-grp'
    clusterOptions = '-G team154-grp'
    
    // time 10.h
    
    // queue 'normal'
    queue 'long'
    // queue 'basement'
    
    // queueSize = 50 // on casm3
    queueSize = 1000 // on farm5
    
    
    // SH_SCRIPT_SIM = "sim_prior_transplant_ex.sh"
    // SH_SCRIPT_SIM = "sim_prior_top_up_transplant_ex.sh"
    
    script:
    """
    echo "projectDir = ${workflow.projectDir}"
    ${workflow.projectDir}/sim_prior_top_up_transplant_ex.sh \\
        -a "${SIM_RUN_INDEX}" \\
        -b "${SOURCE_PATH}" \\
        -c "${LOGS_SIM_PATH}" \\
        -d "${PAIR_ID}" \\
        -e "${MODEL_NAME}" \\
        -f "${SIM_PATH}" \\
        -g "${SIM_TREES_PATH}" \\
        -h "${PRIOR_SAMPLE_DRIFT_FILE}" \\
        -i "${PRIOR_SAMPLE_DRIVERS_FILE}" \\
        -j "${N_SIMS_PER_RUN}"
    """
    
}


workflow {

    
	// create new paths (strings)
	
	PATIENT_OUT_PATH = params.project_out_path + "output_Pair" + params.PAIR_ID + "/"
	
	LOGS_PATH = PATIENT_OUT_PATH + "logs/"
	LOGS_SIM_PATH = PATIENT_OUT_PATH + "logs_sim/"
	SIM_PATH = PATIENT_OUT_PATH + "sims_round_1/"
	SIM_TREES_PATH = PATIENT_OUT_PATH + "sim_trees_round_1/"
	
	println "PATIENT_OUT_PATH = ${PATIENT_OUT_PATH}"
	println "LOGS_PATH = ${LOGS_PATH}"
	println "LOGS_SIM_PATH = ${LOGS_SIM_PATH}"
	println "SIM_PATH = ${SIM_PATH}"
	println "SIM_TREES_PATH = ${SIM_TREES_PATH}"
	
	
	// create directories at these paths
	
	LOGS_DIR = file("${LOGS_PATH}")
	result = LOGS_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${LOGS_DIR}"
	
	LOGS_SIM_DIR = file("${LOGS_SIM_PATH}")
	result = LOGS_SIM_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${LOGS_SIM_DIR}"
	
	SIM_DIR = file("${SIM_PATH}")
	result = SIM_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${SIM_PATH}"
	
	SIM_TREES_DIR = file("${SIM_TREES_PATH}")
	result = SIM_TREES_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${SIM_TREES_PATH}"
	
	
	// run pipeline
	
    // def run_index_list = 1..params.N_RUNS
    // def run_index_list = params.N_RUNS..1000
    def run_index_list = params.N_RUNS_BEGIN..params.N_RUNS_END
    
    def run_ch = Channel.fromList( run_index_list )
	
    run_sims( run_ch, params.SOURCE_PATH, LOGS_SIM_PATH, params.PAIR_ID, params.MODEL_NAME, SIM_PATH, SIM_TREES_PATH, params.PRIOR_SAMPLE_DRIFT_FILE, params.PRIOR_SAMPLE_DRIVERS_FILE, params.N_SIMS_PER_RUN )
	
}




























