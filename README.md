# Clonal_dynamics_of_HSCT
Code accompanying the manuscript 'Clonal dynamics after allogeneic haematopoietic cell transplantation using genome-wide somatic mutations'

# General notes
The code is provided for all analyses from mutation calls from Caveman and Pindel on the original bam files (run unmatched) up to the generation of figures used in the manuscript. The numbered folders go through the different steps in the analysis, roughly in the order in which they were performed. However, (most) intermediate data is available such that each stage of the analysis can be readily performed without having to re-run earlier stages. \
Smaller data objects are provided within the github repository. Larger data objects need to be downloaded from Mendeley Data at doi: 10.17632/m7nz2jk8wb.2 \
Raw whole-genome sequencing data is available via the EGA (https://ega-archive.org/), accession number: EGAD00001010872 \
Raw targeted sequencing data is available via the EGA (https://ega-archive.org/), accession number: EGAD00001010874
If there are any queries, feel free to contact me via email: ms56@sanger.ac.uk

## Objects that need to be downloaded from Mendeley Data
data/Targeted_sequencing_data \
data/annot_files \
data/annot_files_no_dups \
data/HDP/HDP_multi_chain.Rds \
data/genomic_loci_reference_files

# Notes on specific stages of data generation

## 01 Generating the phylogenies and mutation lists
This is potentially the most challenging to reproduce exactly. You will first need to download the original bam files from European Genome/Phenome Archive (EGA), accession number: EGAD00001010872.
You will then need to rerun Caveman & Pindel against an unmatched synthetic reference genome to produce mutation calls from each clone.
You will then need to run the hairpin filtering program that we have developed to remove artefacts from our low-input pipeline data.
You will next need to run cgpVaf, a programme to extract read counts across a large set of bams at genomic positions specified in a bed file.
If you are trying to replicate this pipeline for your own data, it is likely going to be easier to adapt this stage to your local set up and the mutation-calling algorithm that you are most familiar with e.g. Mutect2, Strelka, VarScan. As long as these are run in 'unmatched' mode, all produce output that can be used in subsequent analyses. Please see the following reference for further details:
Coorens, T.H.H., Spencer Chapman, M., Williams, N. et al. Reconstructing phylogenetic trees from genome-wide somatic mutations in clonal samples. Nat Protoc 19, 1866–1886 (2024). https://doi.org/10.1038/s41596-024-00962-8 

## 02 Running_HDP_mutation_signature_extraction
The scripts provided here allow this to be re-run locally, or better, on a compute farm with parallel computing capability. All input data & scripts are provided.

## 03 Compiling_data_for_downstream_analysis
This is a single script that takes the data from the previous two stages, as well as a comprehensive list of sample metadata to generate the final set of data objects for all downstream analyses.

## 04 Gibbs_sampler_for_targ_seq
These scripts are written in julia. They need to be run on data using the three separate baitsets individually. This is because for each baitset, the samples coming from different individuals to the one in which the mutation was called are used to estimate the locus-specific sequencing error rates.

## 05 ABC_simulation_scripts
These scripts must be run on a compute farm, and use nextflow for job submission.
While they are set up for submission on LSF, they could readily be adapted to other systems.

## 06 Generating_figures
Within this folder there are subfolders with scripts for generating each figure from the manuscript.
Each figure uses downstream data files that are saved within the data/ folder (or from Mendeley data)
All generated plots are saved within the plots/ folder.

# OTHER FOLDERS
## Original analysis scripts
These are the original scripts used to analyse the data. They are not quite as 'tidy' as those within the '06 Generating_figures/' folder, and much of this is duplications of the same analyses. However, some scripts contain code for some additional plots and analyses that are not included in the manuscript and are therefore included for reference.

## Other simulation scripts
Simulation scripts used for analyses other than the main 'engrafting cell number' ABC.
This includes:
1. estimating the phylogenetic age \
2. estimating the effect that increased T-cell clone longevity may have on clonal composition relative to the myeloid fraction

## data
data/reference_files/ - includes reference files used in various analyses \
data/metadata_files/ - includes individual-level metadata, and sample-level metadata \
data/tree_and_mutation_files/ - includes all saved objects relating to tree structures or mutation information from the WGS \
data/SV_and_CNA_data/ - includes summaries of structural variants and copy number alterrations from GRIDSS and ASCAT. Loss-of-Y information is derived from the mean coverage data. \
data/HDP/ - data files relating mutational signature extraction with HDP \
data/APOBEC_VCFs - vcf files containing only the likely APOBEC mutations from branches affected by APOBEC/ \
data/targeted_sequencing_data - raw and inferred data related to the targeted sequencing. The 'data_by_baitset' folder includes \
data/ABC_simulation_results - the posterior results from the models with different ABCs

