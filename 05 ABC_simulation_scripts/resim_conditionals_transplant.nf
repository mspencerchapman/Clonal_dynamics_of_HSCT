#!/usr/bin/env nextflow

nextflow.enable.dsl=2 // what does this do?

// inputs
params.SOURCE_PATH = "/source/" // default value (debugging only)
params.project_out_path = "/output/" // default value (debugging only)
params.PAIR_ID = 11 // default value (debugging only)
params.MODEL_NAME = "m1" // default value (debugging only)
params.POST_SAMPLE_FILE = "/abc_output_round_1/abc_rejection_output/test_posterior_sample.txt" // default value (debugging only)
params.BEGIN_POST_OBS = 1 // default value (debugging only)
params.END_POST_OBS = 10 // default value (debugging only)
params.N_RUNS_PER_POST_OBS = 5 // default value (debugging only)
params.N_RESIMS_PER_RUN = 2 // default value (debugging only)
params.N_RUNS_ACCUMULATED_PER_POST_OBS = 0 // default value (debugging only)
params.N_RESIM_RUNS = 8 // default value (debugging only)


process run_sims {

	input:
	val JOB_INDEX
	val SOURCE_PATH
	val LOGS_RESIM_PATH
	val PAIR_ID
	val MODEL_NAME
	val RESIM_PATH
	val RESIM_TREES_PATH
	val POST_SAMPLE_FILE
    val BEGIN_POST_OBS
    val END_POST_OBS
	val N_RUNS_PER_POST_OBS
	val N_RUNS_ACCUMULATED_PER_POST_OBS
	val N_RESIMS_PER_RUN
	
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
    
    
    // SH_SCRIPT_SIM = "resim_conditionals_transplant_ex.sh"
    
    script:
    """
    echo "projectDir = ${workflow.projectDir}"
    ${workflow.projectDir}/resim_conditionals_transplant_ex.sh \\
        -a "${JOB_INDEX}" \\
        -b "${SOURCE_PATH}" \\
        -c "${LOGS_RESIM_PATH}" \\
        -d "${PAIR_ID}" \\
        -e "${MODEL_NAME}" \\
        -f "${RESIM_PATH}" \\
        -g "${RESIM_TREES_PATH}" \\
        -h "${POST_SAMPLE_FILE}" \\
        -i "${BEGIN_POST_OBS}" \\
        -j "${END_POST_OBS}" \\
        -k "${N_RUNS_PER_POST_OBS}" \\
        -l "${N_RUNS_ACCUMULATED_PER_POST_OBS}" \\
        -m "${N_RESIMS_PER_RUN}"
    """
    
}


workflow {

    
	// create new paths (strings)
	
	PATIENT_OUT_PATH = params.project_out_path + "output_Pair" + params.PAIR_ID + "/"
	
	LOGS_PATH = PATIENT_OUT_PATH + "logs/"
	LOGS_SIM_PATH = PATIENT_OUT_PATH + "logs_sim/"
	LOGS_RESIM_PATH = PATIENT_OUT_PATH + "logs_resim/"
	RESIM_PATH = PATIENT_OUT_PATH + "conditional_resims/"
	RESIM_TREES_PATH = PATIENT_OUT_PATH + "conditional_resims_trees/"
	
	println "PATIENT_OUT_PATH = ${PATIENT_OUT_PATH}"
	println "LOGS_PATH = ${LOGS_PATH}"
	println "LOGS_SIM_PATH = ${LOGS_SIM_PATH}"
	println "LOGS_RESIM_PATH = ${LOGS_RESIM_PATH}"
	println "RESIM_PATH = ${RESIM_PATH}"
	println "RESIM_TREES_PATH = ${RESIM_TREES_PATH}"
	
	
	// create directories at these paths
	
	LOGS_DIR = file("${LOGS_PATH}")
	result = LOGS_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${LOGS_DIR}"
	
	LOGS_SIM_DIR = file("${LOGS_SIM_PATH}")
	result = LOGS_SIM_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${LOGS_SIM_DIR}"
	
	LOGS_RESIM_DIR = file("${LOGS_SIM_PATH}")
	result = LOGS_RESIM_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${LOGS_RESIM_DIR}"
	
	RESIM_DIR = file("${RESIM_PATH}")
	result = RESIM_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${RESIM_PATH}"
	
	RESIM_TREES_DIR = file("${RESIM_TREES_PATH}")
	result = RESIM_TREES_DIR.mkdir()
	println result ? "OK" : "Cannot create directory: ${RESIM_TREES_PATH}"
	
	N_POST_OBS = params.END_POST_OBS - params.BEGIN_POST_OBS
	N_POST_OBS = N_POST_OBS + 1
	// N_RESIM_RUNS = params.N_RUNS_PER_POST_OBS * N_POST_OBS
	
	
	// run pipeline
	
    // def run_index_list = 1..N_RESIM_RUNS
    def run_index_list = 1..params.N_RESIM_RUNS
    
    def run_ch = Channel.fromList( run_index_list )
	
    run_sims( run_ch, params.SOURCE_PATH, LOGS_RESIM_PATH, params.PAIR_ID, params.MODEL_NAME, RESIM_PATH, RESIM_TREES_PATH, params.POST_SAMPLE_FILE, params.BEGIN_POST_OBS, params.END_POST_OBS, params.N_RUNS_PER_POST_OBS, params.N_RUNS_ACCUMULATED_PER_POST_OBS, params.N_RESIMS_PER_RUN )
	
}


























