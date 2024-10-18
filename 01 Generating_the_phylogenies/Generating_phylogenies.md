# Generating phylogenies from WGS data

This markdown script runs through all the steps used to generate the phylogenies.
This version is not readily generalizable to environments outside the Wellcome Sanger Institute, as a lot of the algorithms are in-house generated. However, all the software is available.
The generic pipeline is discussed at length within our Nature Protocols paper and this should be your first port of call if trying to adapt this pipeline to your own data and setup: \
Coorens, T.H.H., Spencer Chapman, M., Williams, N. et al. Reconstructing phylogenetic trees from genome-wide somatic mutations in clonal samples. Nat Protoc 19, 1866–1886 (2024). https://doi.org/10.1038/s41596-024-00962-8


## Set the variables that are specific for this analysis
For each individual need to define: \
1. The experiment number (used for naming files throughout the analysis) \
2. The project number in which the bams/ variant annotation files are stored - this is specific to the Sanger insitute, so will need to be adapted if reproducing analysis. \
3. The IDs of samples from the donor, and those from the recipient - used to separate out the donor & recipient phylogenies in the final script. \
4. The age of the donor at the time of sampling - used for scaling the trees to  'actual time' \
5. Any samples that should be exclude as being from a different germline background. Here they are labelled as being recipient chimerism, but in fact are more likely from contamination from other individuals. \

Note that other sets of variables are in the script HSCT_set_variables.sh & can be swapped in.
```bash
EXP_NO=11
ALL_PROJECT_NUMBERS=2189,2256
DONOR_ID=PD45792b
RECIP_ID=PD45793b
DONOR_AGE=74.8
RECIP_CHIMERISM=PD45793b_lo0094,PD45793b_lo0095,PD45793b_lo0096
```

## Load up all necessary modules (this is done in Linux Ubuntu 22)
Some of these are specific to the CASM setup at the Sanger, so may need to be modified. However, most are available from github.com/cancerit

```bash
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
```

#Set reference file paths
```bash
#Reference file paths
GENOME_FILE=/lustre/scratch124/casm/team78pipelines/reference/human/GRCh37d5/genome.fa
HIGH_DEPTH_REGIONS=/lustre/scratch124/casm/team78pipelines/reference/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz
GENES_TO_ANNOTATE=/lustre/scratch126/casm/team154pc/ms56/reference_files/chip_drivers.txt
```

## Program file paths
Define paths to all the different custom scripts used during the analysis.

```bash
PROGRAMS_DIR=/lustre/scratch126/casm/team154pc/ms56/my_programs
IMPORT_NEW_SAMPLES_SCRIPT=${PROGRAMS_DIR}/import_new_samples_only.R
CREATE_SPLIT_CONFIG_SCRIPT=${PROGRAMS_DIR}/create_split_config_ini.R
FILTERING_LCMOUTPUT_SCRIPT=${PROGRAMS_DIR}/filter_for_bed_file.pl
MUTATION_PARAMETER_SCRIPT=${PROGRAMS_DIR}/Mutation_filtering_get_parameters.R
SENSITIVITY_ANALYSIS_SCRIPT=${PROGRAMS_DIR}/Sensitivity_analysis_from_SNPs.R
TREE_BUILDING_SCRIPT=${PROGRAMS_DIR}/filtering_from_table_mix_remove.R
REDUCE_MUTSET_SCRIPT=${PROGRAMS_DIR}/Reducing_mutset_from_MSfilters.R
```

## Set the main study directory, and define relative file paths
Set a main 'study directory' from which all other file paths are set relative to.

```bash
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
```

## Filtering artefacts from the low-input pipeline
The low-input DNA pipeline used for sequencing the colonies is known to be associated with a very specific set of artefacts related to cruciform DNA formation.
A specific algorithm has been designed to filter these and is available from https://github.com/cancerit/hairpin-wrapper

```bash
#--------------------Submit LCM filtering jobs-----------------------------
mkdir -p $ROOT_DIR/MS_filters/output_files
cd $ROOT_DIR/MS_filters
/lustre/scratch126/casm/team154pc/ms56/my_programs/Submitting_Mathijs_filters_jobs.R -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -o $ROOT_DIR/MS_filters/output_files -q long
```


## SNV analysis
This code is to create a bedfile of all SNVs called in at least one sample, and then used cgpVaf to extract wild type and mutation allele read counts from all individual samples.

```bash
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
```

## INDEL analysis
This code is to create a bedfile of all indels called in at least one sample, and then used cgpVaf to extract wild type and mutation allele read counts from all individual samples.

```bash
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
```

## SNV merge
To run once cgpVAF has completed for SNVs
This extracts the relevant columns from the cgpVAF output to create matrices of all samples that only includes the mutant allele & total depth columns.

```bash
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
```
## INDEL merge
To run once cgpVAF has completed for indels
This extracts the relevant columns from the cgpVAF output to create matrices of all samples that only includes the mutant allele & total depth columns.

```bash
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
```

## Generate filtering parameters for subsequent use in the tree-building script
This script iterates across all mutation loci from the original bed file and uses the read count matrices to generate several parameters used for subsequent filtering.
1. Exact binomial test \
2. Beta-binomial test \
3. P-value test \
4. 
5. Maximum depth in positive samples \
6. Maximum VAF in positive samples \

Note that this is run twice: \
(1) not excluding samples, \
(2) excluding colonies with peak VAF <0.4

```bash

#--------------------AFTER ALL CGPVAF MATRICES ARE GENERATED-----------------------------
cd $ROOT_DIR
mkdir -p log_files
mkdir -p err_files

#Run the "Mutation_filtering_get_params" script.

#1. Without excluding any samples
bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q long -G team78-grp -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n6 -J GET_PARAMS \
    Rscript $MUTATION_PARAMETER_SCRIPT \
    -r $RUN_ID \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR

#2. Excluding samples with a peak vaf < 0.4
bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q long -G team78-grp -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
    -M24000 -n6 -J GET_PARAMS \
    Rscript $MUTATION_PARAMETER_SCRIPT \
    -r $RUN_ID_M \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR \
    -m \
    -v 0.4

```

## Reducing the matrix
The matrices and filter parameters from the previous step included all mutations called by Caveman & Pindel, without using the hairpin filtering.
Therefore, we can now exclude any mutation that was filtered by the hairpin filters, and reduce the mutation set accordingly.

```bash
#ONCE MS FILTERS JOBS HAVE FINISHED
cd $ROOT_DIR/MS_filters/output_files
perl $FILTERING_LCMOUTPUT_SCRIPT
cut -f 1,2,4,5 *complete_final_retained_3.vcf_for_bed_file|sort|uniq>$ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME

#ONCE Mutation_filtering_get_paramaters.R SCRIPT HAS COMPLETED
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce Rscript $REDUCE_MUTSET_SCRIPT -r $RUN_ID -b $ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce Rscript $REDUCE_MUTSET_SCRIPT -r $RUN_ID_M -b $ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR
```

## Sensitivity analysis
This script uses the sensitivity for calling germline variants (that are known to be present in all samples) as a surrogate for their sensitivity for calling somatic variants.
This file is then used in the subsequent tree-building analysis for correcting branch lengths.

```bash
#Run the sensitivity analysis
bsub -o $ROOT_DIR/log_files/sensitivity.log.%J -e $ROOT_DIR/err_files/sensitivity.err.%J \
    -q normal -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
    -M4000 -n1 -J sensitivity \
    Rscript $SENSITIVITY_ANALYSIS_SCRIPT \
    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID}_postMS \
    -o $ROOT_DIR -n sensitivity_analysis_${EXP_ID} \
    -i $ROOT_DIR/pindel_raw \
    -s $ROOT_DIR/MS_filters/output_files \
    -x '_complete_final_retained_3.vcf_for_bed_file'
```

## Run the tree-building script

This is now the business end of the analysis & is the main script that generates the phylogenies and the final set of somatic mutations.
This script does quite a lot of work and has many options.

|Option (short)|Option (long)|Description|
|---|-----|------------------------|
|-w|--working_dir|Working directory for analysis, if not set will default to the current directory|
|-i|--id_run|Run ID for this filtering run|
|-m|--mats_and_params|path for the mats and params file, if not set will default to the usual filename within the output directory|
|-f|--filtering_type|Set as pval or vaf to choose filtering type|
|-c|--covcut|Remove samples with mean coverage below specified cut-off|
|-d|--donor_id|ID for the transplant donor|
|-r|--recip_id|ID for the transplant recipient|
|-o|--output_dir|/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs|
|-s|--sensitivity_path|Path for the saved sensitivity dataframe|
|-x|--xclude_samples|option to manually exclude certain samples from the analysis,separate with a comma|
|-e|--exclude_muts|option to manually exclude mutations,even if pass filtering|
|-k|--keep_muts|option to manually retain mutation,even if fail filtering|
|-p|--polytomous_tree|option to make the tree polytomous i.e. multi-furcating|
|-b|--bbcutoff|cut-off value for beta-binomial filter|
|-y|--y_filter|do not use Y mutations, or the X chromosome PARs for tree-building. Useful if recurrent loss of Y is likely.|
|-a|--ancestral|option to keep the ancestral branch for mutation allocation|
|-n|--nonclonal|option to switch off the removal of non clonal samples by testing against tree|
|-g|--germline_addback|option to switch off the 'check_for_false_germline_calls' function|
|-j|--just_snvs|option to only use the SNVs for tree building, useful if high indel error rate|
|-l|--loh_and_dels|table of copy number losses to incorporate for tree building|
|-q|--qgenes|path to list of genes to annotate on the trees (e.g. driver genes)|
|-t|--time|donor age at time of sampling for scaling of ultrametric tree|

Now run the script on as a job on LSF using the following commands.
These are the options used to generate the data in the manuscript.

```bash
#SETTING OFF THE TREE-BUILDING SCRIPT
RUN_ID=${EXP_ID}_m40
RUN_ID_TB=${RUN_ID}_postMS_reduced

###IF THERE IS ARE SAMPLES FROM A DIFFERENT GERMLINE BACKGROUND, THESE SHOULD BE REMOVED USING the -x flag
MEM=24000
bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q basement -R "select[type==X86_64 && mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" \
    -M${MEM} -n1 -J tree_build \
    Rscript $TREE_BUILDING_SCRIPT \
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
    -q basement -R "select[type==X86_64 && mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" \
    -M${MEM} -n1 -J tree_build \
    Rscript $TREE_BUILDING_SCRIPT \
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
```

## Alternatively, if you want to include samples with a different germline
Can rerun the tree without excluding samples that are from a different germline background.
This is only useful if you want to get a list of the germline SNPs in the samples that have a different germline e.g. to see if they come from another individual from the study (due to cross-contamination) or are likely to be residual recipient chimerism.

```bash
MEM=24000
bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" \
    -M${MEM} -n1 -J tree_build \
    Rscript $TREE_BUILDING_SCRIPT \
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
    -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" \
    -M${MEM} -n1 -J tree_build \
    Rscript $TREE_BUILDING_SCRIPT \
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
    
```