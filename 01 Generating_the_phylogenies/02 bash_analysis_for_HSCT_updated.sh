#Load all the necessary modules
module load cgpVAFcommand
module load mpboot
module load perl
module load vagrent
module load tabix
module load pcap-core
module load dataImportExport
module load picard-tools
module load canPipe/live
module load canned-queries-client
module load bwa

#Reference file paths
GENOME_FILE=/lustre/scratch124/casm/team78pipelines/reference/human/GRCh37d5/genome.fa
HIGH_DEPTH_REGIONS=/lustre/scratch124/casm/team78pipelines/reference/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz
GENES_TO_ANNOTATE=/lustre/scratch126/casm/team154pc/ms56/reference_files/chip_drivers.txt

#Program file paths
PROGRAMS_DIR=/lustre/scratch126/casm/team154pc/ms56/my_programs
IMPORT_NEW_SAMPLES_SCRIPT=${PROGRAMS_DIR}/import_new_samples_only.R
CREATE_SPLIT_CONFIG_SCRIPT=${PROGRAMS_DIR}/create_split_config_ini.R
FILTERING_LCMOUTPUT_SCRIPT=${PROGRAMS_DIR}/filter_for_bed_file.pl
MUTATION_PARAMETER_SCRIPT=${PROGRAMS_DIR}/Mutation_filtering_get_parameters.R
SENSITIVITY_ANALYSIS_SCRIPT=${PROGRAMS_DIR}/Sensitivity_analysis_from_SNPs.R
TREE_BUILDING_SCRIPT=${PROGRAMS_DIR}/filtering_from_table_mix_remove.R

#Set the main study directory
STUDY_DIR=/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT

#Secondary variables
PD_NUMBERS=${DONOR_ID},${RECIP_ID}
MATS_AND_PARAMS_DIR=${STUDY_DIR}/filtering_runs/mats_and_params
ROOT_DIR=${STUDY_DIR}/Pair_${EXP_NO}
EXP_ID=Pair${EXP_NO}
SNV_BEDFILE_NAME=Zur_HSCT_${EXP_ID}_caveman.bed
INDEL_BEDFILE_NAME=Zur_HSCT_${EXP_ID}_pindel.bed
MS_FILTERED_BEDFILE_NAME=Zur_HSCT_${EXP_ID}_postMS_SNVs.bed
IFS=',' read -r -a PROJECT_ARRAY <<< "$ALL_PROJECT_NUMBERS"
RUN_ID=$EXP_ID
RUN_ID_M=${EXP_ID}_m40
RUN_ID_TB=${RUN_ID_M}_postMS_reduced

#--------------------Submit LCM filtering jobs-----------------------------
mkdir -p $ROOT_DIR/MS_filters/output_files
cd $ROOT_DIR/MS_filters
/lustre/scratch126/casm/team154pc/ms56/my_programs/Submitting_Mathijs_filters_jobs.R -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -o $ROOT_DIR/MS_filters/output_files -q long

#----------------SNV analysis-----------------
mkdir -p $ROOT_DIR/caveman_raw/caveman_pileup/output; cd $ROOT_DIR/caveman_raw
Rscript $IMPORT_NEW_SAMPLES_SCRIPT -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -t SNVs
cut -f 1,2,4,5 *.caveman_c.annot.vcf_pass_flags|sort|uniq>caveman_pileup/$SNV_BEDFILE_NAME
cd $ROOT_DIR/caveman_raw/caveman_pileup

#Run createVafCmd for each project id containing samples
for PROJECT_NUMBER in "${PROJECT_ARRAY[@]}"; do
    echo "3"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $SNV_BEDFILE_NAME
done

#Then run the "create_split_config_ini.R" script in the output folder
PROJECT_NUMBER=${PROJECT_ARRAY[0]}
cd $ROOT_DIR/caveman_raw/caveman_pileup/output
Rscript $CREATE_SPLIT_CONFIG_SCRIPT -p $PROJECT_NUMBER
cd $ROOT_DIR/caveman_raw/caveman_pileup

#Then re-run the createVafCmd.pl script with the new config file
echo "3"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -i output/${PROJECT_NUMBER}_cgpVafConfig_split.ini -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $SNV_BEDFILE_NAME

#Update the run_bsub.sh command to allow more jobs in the array to run together and get more memory
sed -e 's/\%5/\%200/g;s/2000/4000/g;s/500/1000/g;' run_bsub.sh >run_bsub_updated.sh

#Run this if need to switch to the long queue (not normally necessary)
#sed -i -e 's/normal/long/g' run_bsub_updated.sh

bash run_bsub_updated.sh

#----------------INDEL analysis-----------------
mkdir -p $ROOT_DIR/pindel_raw/pindel_pileup/output; cd $ROOT_DIR/pindel_raw
mkdir -p $ROOT_DIR/pindel_raw/pindel_pileup_FF016
Rscript $IMPORT_NEW_SAMPLES_SCRIPT -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -t indels
cut -f 1,2,4,5 *.pindel.annot.vcf_pass_flags|sort|uniq>pindel_pileup/$INDEL_BEDFILE_NAME

cd pindel_pileup
#Run createVafCmd for each project id containing samples
for PROJECT_NUMBER in "${PROJECT_ARRAY[@]}"; do
    echo "1"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $SNV_BEDFILE_NAME
done

#Then run the "create_split_config_ini.R" script in the output folder
PROJECT_NUMBER=${PROJECT_ARRAY[0]}
cd $ROOT_DIR/pindel_raw/pindel_pileup/output
Rscript $CREATE_SPLIT_CONFIG_SCRIPT -p $PROJECT_NUMBER
cd $ROOT_DIR/pindel_raw/pindel_pileup

#Then re-run the createVafCmd.pl script with the new config file
echo "1"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -i output/${PROJECT_NUMBER}_cgpVafConfig_split.ini -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $INDEL_BEDFILE_NAME

#Update the run_bsub.sh command to allow more jobs in the array to run together and get more memory
sed -e 's/\%5/\%200/g;s/2000/4000/g;s/500/1000/g;' run_bsub.sh >run_bsub_updated.sh

#Run this if need to switch to the long queue (not normally necessary)
#sed -i -e 's/normal/long/g' run_bsub_updated.sh

bash run_bsub_updated.sh

#--------------------ONCE cgpVAF HAS COMPLETED-----------------------------
#-----------------------------SNV merge-----------------------------
cd $ROOT_DIR/caveman_raw/caveman_pileup/output/output/PDv37is/snp
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,24,26,39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 $(sed -n '1p' files) > temp.1   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do
    cut -f 39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 $FILE >  temp.$FILE
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_SNVs_${EXP_ID}.tsv && rm output.*

mv $ROOT_DIR/caveman_raw/caveman_pileup/output/output/PDv37is/snp/merged_SNVs_${EXP_ID}.tsv $ROOT_DIR/

#-----------------------------INDEL merge-----------------------------
cd $ROOT_DIR/pindel_raw/pindel_pileup/output/output/PDv37is/indel
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,16,18,25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117 $(sed -n '1p' files) > temp.1   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do
    cut -f 25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117 $FILE >  temp.$FILE
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_indels_${EXP_ID}.tsv && rm output.*
mv $ROOT_DIR/pindel_raw/pindel_pileup/output/output/PDv37is/indel/merged_indels_${EXP_ID}.tsv $ROOT_DIR/

#--------------------AFTER ALL CGPVAF MATRICES ARE GENERATED-----------------------------
cd $ROOT_DIR
mkdir -p log_files
mkdir -p err_files

#Run the "Mutation_filtering_get_params" script. Run twice: (1) not excluding samples, (2) excluding colonies with peak VAF <0.4
bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q long -G team78-grp -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n6 -J GET_PARAMS \
    $MUTATION_PARAMETER_SCRIPT \
    -r $RUN_ID \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR
bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q long -G team78-grp -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n6 -J GET_PARAMS \
    $MUTATION_PARAMETER_SCRIPT \
    -r $RUN_ID_M \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR \
    -m \
    -v 0.4

#ONCE MS FILTERS JOBS HAVE FINISHED
cd $ROOT_DIR/MS_filters/output_files
perl $FILTERING_LCMOUTPUT_SCRIPT
cut -f 1,2,4,5 *complete_final_retained_3.vcf_for_bed_file|sort|uniq>$ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME

#ONCE Mutation_filtering_get_paramaters.R SCRIPT HAS COMPLETED
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce /lustre/scratch119/casm/team154pc/ms56/my_programs/Reducing_mutset_from_MSfilters.R -r $RUN_ID -b $ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce /lustre/scratch119/casm/team154pc/ms56/my_programs/Reducing_mutset_from_MSfilters.R -r $RUN_ID_M -b $ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR

#Run the sensitivity analysis
bsub -o $ROOT_DIR/log_files/sensitivity.log.%J -e $ROOT_DIR/err_files/sensitivity.err.%J \
    -q normal -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
    -M4000 -n1 -J sensitivity \
    $SENSITIVITY_ANALYSIS_SCRIPT \
    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID}_postMS \
    -o $ROOT_DIR -n sensitivity_analysis_${EXP_ID} \
    -i $ROOT_DIR/pindel_raw \
    -s $ROOT_DIR/MS_filters/output_files \
    -x '_complete_final_retained_3.vcf_for_bed_file'

#SETTING OFF THE TREE-BUILDING SCRIPT
RUN_ID=${EXP_ID}_m40
RUN_ID_TB=${RUN_ID}_postMS_reduced

#-----------------------Run the tree-building script - this has lots of options

# -i The id for the tree-building run - will be included in the output file names.
# -m Path to mutation filtering output file
# -f Option to do p-value based filtering, or vaf based filtering (pval or vaf)
# -c Cut off to exclude samples with low coverage
# -o Folder for storing script output files. Folder & subfolders will be created if don't already exist.
# -s Path to the sensitivity dataframe
# -p Save trees as polytomous trees (not bifurcating trees)
# -t The age ('time') of the individual - for building the age-adjusted ultrametric tree
# -j Option to do initial tree-building with just the SNVs (i.e. don't' include indels)
# -a Option to keep an ancestral branch

###IF THERE IS ANY RECIPIENT CHIMERISM THIS SHOULD BE SET AND REMOVED
bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q basement -R 'select[type==X86_64 && mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n1 -J tree_build \
    $TREE_BUILDING_SCRIPT \
    -i ${RUN_ID_TB} \
    -m ${STUDY_DIR}/mats_and_params/mats_and_params_${RUN_ID_TB} \
    -f pval \
    -d $DONOR_ID \
    -r $RECIP_ID \
    -c 4 \
    -o ${STUDY_DIR}/filtering_runs \
    -s ${STUDY_DIR}/Pair_${EXP_NO}/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -a \
    -x $RECIP_CHIMERISM

bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q basement -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n1 -J tree_build \
    $TREE_BUILDING_SCRIPT \
    -i ${RUN_ID_TB} \
    -m ${STUDY_DIR}/filtering_runs/mats_and_params/mats_and_params_${RUN_ID_TB} \
    -f vaf \
    -d $DONOR_ID \
    -r $RECIP_ID \
    -c 4 \
    -o ${STUDY_DIR}/filtering_runs \
    -s ${STUDY_DIR}/Pair_${EXP_NO}/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -a \
    -x $RECIP_CHIMERISM


#IF THERE IS RECIPIENT CHIMERISM AND WANT TO INCLUDE THIS

bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q long -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n1 -J tree_build \
    $TREE_BUILDING_SCRIPT \
    -i ${RUN_ID_TB}_wr \
    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID_TB} \
    -f pval \
    -d $DONOR_ID \
    -r $RECIP_ID \
    -c 4 \
    -o ${STUDY_DIR}/filtering_runs \
    -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -a

bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q long -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n1 -J tree_build \
    $TREE_BUILDING_SCRIPT \
    -i ${RUN_ID_TB}_wr \
    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID_TB} \
    -f vaf \
    -d $DONOR_ID \
    -r $RECIP_ID \
    -c 4 \
    -o ${STUDY_DIR}/filtering_runs \
    -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -a