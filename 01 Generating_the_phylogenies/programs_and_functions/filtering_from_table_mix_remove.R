#!/software/R-4.1.0/bin/Rscript

#This script takes the output from the "Mutation_filtering_get_parameters.R" scripts and performs several steps
#(1) Filters the mutations, can use two primary strategies - either based around the mutation VAF or a binomial test of the read counts to test for consistency with a somatic mutation
#(2) Runs the tree building, using mpboot
#(3) Assign mutations to each tree branch using treemut. An "ancestral" branch can optionally be retained for this stage (to which unfiltered germline or artefact mutations may be assigned)
#(4) An additional step to filter samples that are mixed colonies, by testing the vafs against the initial tree
#(5) Runs the tree building a second time without the mixed colonies and again reassigns the mutations to each branch (mutations that are "private" to the mixed colonies are removed)

#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("optparse","ggplot2","dplyr","ape","dichromat","seqinr","stringr","readr","phytools","tidyr","data.table")
bioconductor_packages=c("ggtree")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

#========================================#
# Define the options list from the command line ####
#========================================#

option_list = list(
  make_option(c("-w", "--working_dir"), action="store", default=NULL, type='character', help="Working directory for analysis, if not set will default to the current directory"),
  make_option(c("-i", "--id_run"), action="store", default='HSPC_filter', type='character', help="Run ID for this filtering run"),
  make_option(c("-m", "--mats_and_params"), action="store", default=NULL, type='character', help="path for the mats and params file, if not set will default to the usual filename within the output directory"),
  make_option(c("-f", "--filtering_type"), action="store", type='character', help="Set as pval or vaf to choose filtering type"),
  make_option(c("-c", "--covcut"), action="store", default=0, type='numeric', help="Remove samples with mean coverage below specified cut-off"),
  make_option(c("-d", "--donor_id"), action="store", default=NULL, type='character', help="ID for the transplant donor"),
  make_option(c("-r", "--recip_id"), action="store", default='none', type='character', help="ID for the transplant recipient"),
  make_option(c("-o", "--output_dir"), action="store", default="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs", type='character', help="output directory for files"),
  make_option(c("-s", "--sensitivity_path"), action="store", default=FALSE, type='character', help="Path for the saved sensitivity dataframe"),
  make_option(c("-x", "--xclude_samples"), action="store", default=NULL, type='character', help="option to manually exclude certain samples from the analysis,separate with a comma"),
  make_option(c("-e", "--exclude_muts"), action="store", default=NULL, type='character', help="option to manually exclude mutations,even if pass filtering"),
  make_option(c("-k", "--keep_muts"), action="store", default=NULL, type='character', help="option to manually retain mutation,even if fail filtering"),
  make_option(c("-p", "--polytomous_tree"), action="store_true", default=FALSE, type='logical', help="option to make the tree polytomous i.e. multi-furcating"),
  make_option(c("-b", "--bbcutoff"), action="store", default=0.1, type='numeric', help="cut-off value for beta-binomial filter"),
  make_option(c("-y", "--y_filter"), action="store_true", default=FALSE, type='logical', help="do not use Y mutations, or the X chromosome PARs for tree-building. Useful if recurrent loss of Y is likely."),
  make_option(c("-a", "--ancestral"), action="store_true", default=FALSE, type='logical', help="option to keep the ancestral branch for mutation allocation"),
  make_option(c("-n", "--nonclonal"), action="store_false", default=TRUE, type='logical', help="option to switch off the removal of non clonal samples by testing against tree"),
  make_option(c("-g", "--germline_addback"), action="store_false", default=TRUE, type='logical', help="option to switch off the 'check_for_false_germline_calls' function"),
  make_option(c("-j", "--just_snvs"), action="store_true", default=FALSE, type='logical', help="option to only use the SNVs for tree building, useful if high indel error rate"),
  make_option(c("-l", "--loh_and_dels"), action="store", default=NULL, type='character', help="table of copy number losses to incorporate for tree building"),
  make_option(c("-q", "--qgenes"), action="store", default="/lustre/scratch126/casm/team154pc/ms56/reference_files/chip_drivers.txt", type='character', help="file of genes to annotate, one per line"),
  make_option(c("-t", "--time"), action="store", default=NA, type='numeric', help="donor age at time of sampling for scaling of ultrametric tree")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

print(opt)

##SET RUN_ID AND FILEPATHS
Run_ID <- opt$i
filtering_ID <- opt$f
donor_ID<-opt$d
recip_ID<-opt$r
output_dir<-opt$o
sensitivity_df_path<-opt$s
min_sample_mean_cov<-opt$c
if(is.null(opt$x)) {other_samples_to_remove<-NULL} else {other_samples_to_remove<-unlist(strsplit(x=opt$x,split = ","))}
keep_ancestral<-opt$a
create_multi_tree<-opt$p
filter_y<-opt$y
bbcutoff<-opt$b
donor_age<-opt$t
use_SNV_only<-opt$j
if(is.null(opt$r)) {retain_muts<-NULL} else {retain_muts<-unlist(strsplit(x=opt$r,split = ","))}
if(is.null(opt$e)) {exclude_muts<-NULL} else {exclude_muts<-unlist(strsplit(x=opt$e,split = ","))}
if(is.null(opt$m)) {mats_and_params_file<-paste0(output_dir,"/mats_and_params/mats_and_params_", Run_ID)} else {mats_and_params_file<-opt$m}
if(is.null(opt$w)) {my_working_directory<-getwd()} else {my_working_directory<-opt$w}
if(is.null(opt$l)) {CN_table<-NULL} else {CN_table<-read.csv(opt$l,stringsAsFactors = F)}
if(is.null(opt$q)) {genes_to_annotate<-NULL} else {genes_to_annotate = readLines(opt$q)}

print(str(CN_table))

#Other parameters used in script
treefit_pval_cutoff<-1e-3
XY_low_depth_cutoff = 3; XY_high_depth_cutoff = 10; AUTO_low_depth_cutoff = 6; AUTO_high_depth_cutoff = 23


#========================================#
# SET RUN_ID AND FILEPATHS --------------
#========================================#
#Check through these to make sure correct for your system

R_scripts_dir="/lustre/scratch126/casm/team154pc/ms56/my_functions"
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"

#Import functions/ files
setwd(my_working_directory)
R_function_files = list.files(R_scripts_dir,pattern=".R",full.names=TRUE)
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Create the required sub-directories if they are not already present
output_subdirectories=c("mats_and_params","filtered_muts","mpboot_files","tree_files","annotated_muts","mutation_vcfs","plots")
sapply(output_subdirectories,function(folder){dir.create(paste0(output_dir,"/",folder),recursive = T)})

#Set file paths for saved files
filtered_muts_file = paste0(output_dir,"/filtered_muts/filtered_muts_",Run_ID,"_",filtering_ID)
dna_string_file = paste0(output_dir,"/mpboot_files/Filter_", Run_ID,"_", filtering_ID,".fa")
mpboot_tree_file = paste0(dna_string_file,".treefile")
tree_v1_file_path = paste0(output_dir,"/tree_files/tree_", Run_ID,"_",filtering_ID, "_v1.tree")
tree_file_path = paste0(output_dir,"/tree_files/tree_", Run_ID,"_",filtering_ID, ".tree")
file_annot = paste0(output_dir,"/annotated_muts/annotated_mut_set_", Run_ID,"_",filtering_ID)
dna_string_file_post_mix = paste0(output_dir,"/mpboot_files/Filter_", Run_ID,"_", filtering_ID,"_post_mix.fa")
mpboot_tree_file_post_mix = paste0(dna_string_file_post_mix,".treefile")
tree_file_path_post_mix = paste0(output_dir,"/tree_files/tree_", Run_ID,"_",filtering_ID, "_post_mix.tree")
file_annot_post_mix = paste0(output_dir,"/annotated_muts/annotated_mut_set_", Run_ID,"_",filtering_ID,"_post_mix")
tree_file_path_post_mix_post_dup = paste0(output_dir,"/tree_files/tree_", Run_ID,"_",filtering_ID, "_post_mix_post_dup.tree")
file_annot_post_mix_post_dup = paste0(output_dir,"/annotated_muts/annotated_mut_set_", Run_ID,"_",filtering_ID,"_post_mix_post_dup")

#Paths for running vagrent
vcf_header_path = "/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt"
vcf_path = paste0(output_dir,"/mutation_vcfs/mutations_all_", Run_ID,"_",filtering_ID,".vcf")
vagrent_input_path = paste0(output_dir,"/mutation_vcfs/mutations_all_", Run_ID,"_",filtering_ID,"_header.vcf")
vagrent_output_path = paste0(vagrent_input_path,".annot")

#Paths for running vagrent
vcf_header_path = "/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt"
vagrent_cache_path="/lustre/scratch124/casm/team78pipelines/reference/human/GRCh37d5/vagrent/e75/vagrent.cache.gz"
species="Human"
genome_build="GRCm38"
vcf_path = paste0(output_dir,"/mutation_vcfs/mutations_all_", Run_ID,"_",filtering_ID,".vcf")
vagrent_input_path = paste0(output_dir,"/mutation_vcfs/mutations_all_", Run_ID,"_",filtering_ID,"_header.vcf")
vagrent_output_path = paste0(vagrent_input_path,".annot")

#========================================#
# START ANALYSIS: Load mats and params files ----
#========================================#

load(mats_and_params_file)

#Remove the low coverage samples and their private mutations
if(min_sample_mean_cov > 0|!is.null(other_samples_to_remove)) {
  output = remove_low_coverage_samples(COMB_mats = COMB_mats,
                                       filter_params = filter_params,
                                       min_sample_mean_cov = min_sample_mean_cov,
                                       other_samples_to_remove = other_samples_to_remove,
                                       min_variant_reads_auto = 3, #these parameters are to remove mutations from the matrix that are no longer positive in any samples
                                       min_variant_reads_xy = 2)
  COMB_mats= output$COMB_mats
  filter_params = output$filter_params
}

if(filtering_ID=="pval") {
  min_pval_for_true_somatic=0.1
  min_pval_for_true_somatic_SHARED=0.05
  min_vaf=NA
  min_vaf_SHARED=NA
} else if (filtering_ID=="vaf") {
  min_pval_for_true_somatic=NA
  min_pval_for_true_somatic_SHARED=NA
  min_vaf=c(0.2,0.4)
  min_vaf_SHARED=c(0.15,0.30)
} else {
  stop(print("Unrecognised filtering type"))
}

#========================================#
# MUTATION FILTERING 1 ----
#========================================#

#REVIEW MEAN DEPTH HISTOGRAMS TO DECIDE MANUAL CUT-OFFS FOR EXCLUDING OUTLIERS
#hist(filter_params$mean_depth, breaks = 100, xlim = c(0,60))
#hist(filter_params$mean_depth[COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 400, xlim = c(0,30))
#hist(filter_params$mean_depth[!COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 400, xlim = c(0,30))

#This is the main function - applies the set cut-offs to the mutation set, filtering the mutations, assigning genotypes for each to all the samples, and building the dummy dna strings for tree building.
filtered_muts = get_filtered_mut_set(input_set_ID = Run_ID,  #the Run_ID of the unfiltered mutation set used as input - though won't account for any removal of samples from set
                                     COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                     filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                     gender = COMB_mats$gender, #patient's gender
                                     
                                     #These parameters decide whether a mutation is retained in the "true somatic mutation" set
                                     retain_muts = retain_muts,  #any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                     exclude_muts = exclude_muts,
                                     germline_pval = -10,  #the log10 p-value cutoff for mutations coming from an expected germline distribution
                                     rho = bbcutoff,  #rho cutoff for the beta-binomial filter, a measure of how "over-dispersed" the counts are compared to a binomial distribution
                                     mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order 1. lower threshold for autosomes, 2. upper threshold for autosomes, 3. lower threshold for XY, 4. upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                     pval_dp2=NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads. Suitable for low coverage samples only.
                                     pval_dp3=0.001,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping)
                                     min_depth = c(6,4), #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                     min_pval_for_true_somatic = min_pval_for_true_somatic,   #Default: 0.1. the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic
                                     min_vaf = min_vaf, #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY)
                                     
                                     #These parameters decide the genotype for each sample for each "true somatic mutation".  These may be less stringent than the initial parameters.
                                     min_variant_reads_SHARED = 2,  #the minimum number of reads for samples to be assigned a positive genotype
                                     min_pval_for_true_somatic_SHARED = min_pval_for_true_somatic_SHARED,  #the p-value for coming from "true somatic mutation" read distribution to be assigned a positive genotype
                                     min_vaf_SHARED = min_vaf_SHARED) #Numeric vector of length 2, defining minimum vaf to be assigned a positive genotype

if(filter_y) {
  print("Removing variants on the Y chromosome and X chromosome PARs from the dummy dna strings used for tree building")
  #Extract the chromosome and position info from the rows of the "Genotype_shared_bin" object used for tree building
  shared_muts_df=data.frame(mut_ref=rownames(filtered_muts$Genotype_shared_bin))
  shared_muts_df$Chrom=str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,1]
  shared_muts_df$Pos=as.integer(str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,2])
  
  #Remove the variant from the genotype_shared_bin if it is on the Y chromosome or the PAR of the X chromosome
  retain_mut=ifelse(shared_muts_df$Chrom=="Y"|(shared_muts_df$Chrom=="X" & shared_muts_df$Pos %in% c(60001:2699520,154931044:155260560)),FALSE,TRUE)
  filtered_muts$Genotype_shared_bin<-filtered_muts$Genotype_shared_bin[retain_mut,]
  filtered_muts$dna_strings=dna_strings_from_genotype(filtered_muts$Genotype_shared_bin)
}

#Update the genotype table and DNA strings to take account of LOH and del events that make the genotype unclear
if(!is.null(CN_table)) {
  #Extract the chromosome and position info from the rows of the "Genotype_shared_bin" object used for tree building
  shared_muts_df=data.frame(mut_ref=rownames(filtered_muts$Genotype_shared_bin))
  shared_muts_df$Chrom=str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,1]
  shared_muts_df$Pos=as.integer(str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,2])
  
  for(k in 1:nrow(CN_table)){
    this_sample<-CN_table$Sample[k]
    if(CN_table$Type[k]%in%c("LOH","DEL")) {
      current_genotypes<-filtered_muts$Genotype_shared_bin[shared_muts_df$Chrom==CN_table$Chrom[k] &
                                                             shared_muts_df$Pos>CN_table$Pos_min[k]&
                                                             shared_muts_df$Pos<CN_table$Pos_max[k],this_sample]
      new_genotypes=pmax(0.5,current_genotypes)
      cat(paste("Altering genotypes at",length(new_genotypes),"locations in",this_sample),sep = "\n")
      filtered_muts$Genotype_shared_bin[shared_muts_df$Chrom==CN_table$Chrom[k] &
                                          shared_muts_df$Pos>CN_table$Pos_min[k]&
                                          shared_muts_df$Pos<CN_table$Pos_max[k],this_sample] <-new_genotypes
    }
  }
  filtered_muts$dna_strings=dna_strings_from_genotype(filtered_muts$Genotype_shared_bin)
}

if(use_SNV_only) {
  filtered_muts$dna_strings=dna_strings_from_genotype(filtered_muts$Genotype_shared_bin[sapply(rownames(filtered_muts$Genotype_shared_bin),is.snv),])
}

#Decide an ID for this filtered set, depending on approach taken, and save
save(filtered_muts, file = filtered_muts_file)

#Write a fasta file of the dummy DNA strings - this can then be used by mpboot for treebuilding
write.fasta(filtered_muts$dna_strings, names=names(filtered_muts$dna_strings), dna_string_file)

#========================================#
# BUILD INITIAL TREE with MPBoot ----
#========================================#

system(paste0("mpboot -s ", dna_string_file," -bb 1000"))

## Import the tree into R using ape ----
tree <- read.tree(mpboot_tree_file)

## Assign mutations to branches using treemut ----
res=assign_mutations_to_branches(tree,filtered_muts,keep_ancestral=keep_ancestral,create_multi_tree=create_multi_tree,treefit_pval_cutoff=treefit_pval_cutoff)

#Add node and pval information to the filtered_muts object
print("Storing node and p-value information to the info dataframe")
filtered_muts$COMB_mats.tree.build$mat$node <- res$tree$edge[res$summary$edge_ml,2]
filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval

tree<-res$tree #Make the output of the res object the new 'current tree'

#========================================#
# REMOVE POORLY FITTING INDELS ----
#========================================#
#Remove mutations that:
#(1) Are indels, AND
#(2) fit the tree poorly AND have borderline rho values in the beta-binomial filter, OR
#(3) have been assigned to the ancestral branch - if this has been included
low_pval_indels<-Reduce(intersect,list(filtered_muts$COMB_mats.tree.build$mat$mut_ref[filtered_muts$COMB_mats.tree.build$mat$pval<1e-5 & filtered_muts$COMB_mats.tree.build$mat$Mut_type=="INDEL"],rownames(filter_params[filter_params[,"bb_rhoval"]<0.15,])))
if("Ancestral"%in%tree$tip.label) {
  ROOT=tree$edge[1,1]
  ancestral_branch=tree$edge[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral"),2]
  ancestral_branch_indels=filtered_muts$COMB_mats.tree.build$mat$mut_ref[filtered_muts$COMB_mats.tree.build$mat$node==ancestral_branch & filtered_muts$COMB_mats.tree.build$mat$Mut_type=="INDEL"]
  remove_muts=unique(c(ancestral_branch_indels,low_pval_indels))
} else {
  remove_muts<-low_pval_indels
}

if(length(remove_muts)!=0) {
  exclude_muts<-c(exclude_muts,remove_muts)
  print(paste("Removing",length(remove_muts),"mutations due to poor fit with tree and borderline rho values < 0.15"))
  print(remove_muts)
  pdf(file = paste0(output_dir,"/plots/",Run_ID,"_",filtering_ID,"_muts_removed_as_artefact.pdf"),width=15,height=7)
  temp=sapply(remove_muts,function(mut) {
    tree=plot_tree(tree,cex.label = 0)
    add_annotation(tree=tree,
                   details=filtered_muts$COMB_mats.tree.build$mat,
                   matrices=list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR),
                   annot_function = plot_MAV_mut,
                   mut1=mut)
  })
  dev.off()
}

#========================================#
# CHECK FOR POTENTIAL FALSE GERMLINE CALLS ----
#========================================#

#Check for false germline calls, if this option is selected
if(opt$g) {
  false_germline_muts=unlist(check_for_false_germline_calls(tree, COMB_mats,filter_params,max_clade_prop = 0.1,CN_table=CN_table))
} else {
  false_germline_muts<-NULL
}

if(!is.null(false_germline_muts)|length(remove_muts)!=0){
  print("Saving the initial version of the tree - made including the artefactual mutations and without false-germline calls")
  write.tree(tree,file=tree_v1_file_path)
  
  print("Rerunning mutation filtering and tree building including the false germline calls and removing likely artefactual mutations")
  filtered_muts = get_filtered_mut_set(input_set_ID = Run_ID,  #the Run_ID of the unfiltered mutation set used as input - though won't account for any removal of samples from set
                                       COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       gender = COMB_mats$gender, #patient's gender
                                       
                                       #These parameters decide whether a mutation is retained in the "true somatic mutation" set
                                       retain_muts = c(retain_muts,false_germline_muts),  #any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                       exclude_muts = exclude_muts,
                                       germline_pval = -10,  #the log10 p-value cutoff for mutations coming from an expected germline distribution
                                       rho = bbcutoff,  #rho cutoff for the beta-binomial filter, a measure of how "over-dispersed" the counts are compared to a binomial distribution
                                       mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order 1. lower threshold for autosomes, 2. upper threshold for autosomes, 3. lower threshold for XY, 4. upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                       pval_dp2=NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads. Suitable for low coverage samples only.
                                       pval_dp3=0.001,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping)
                                       min_depth = c(6,4), #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                       min_pval_for_true_somatic = min_pval_for_true_somatic,   #Default: 0.1. the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic
                                       min_vaf = min_vaf, #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY)
                                       
                                       #These parameters decide the genotype for each sample for each "true somatic mutation".  These may be less stringent than the initial parameters.
                                       min_variant_reads_SHARED = 2,  #the minimum number of reads for samples to be assigned a positive genotype
                                       min_pval_for_true_somatic_SHARED = min_pval_for_true_somatic_SHARED,  #the p-value for coming from "true somatic mutation" read distribution to be assigned a positive genotype
                                       min_vaf_SHARED = min_vaf_SHARED) #Numeric vector of length 2, defining minimum vaf to be assigned a positive genotype
  
  if(filter_y) {
    print("Removing variants on the Y chromosome and X chromosome PARs from the dummy dna strings used for tree building")
    #Extract the chromosome and position info from the rows of the "Genotype_shared_bin" object used for tree building
    shared_muts_df=data.frame(mut_ref=rownames(filtered_muts$Genotype_shared_bin))
    shared_muts_df$Chrom=str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,1]
    shared_muts_df$Pos=as.integer(str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,2])
    
    #Remove the variant from the genotype_shared_bin if it is on the Y chromosome or the PAR of the X chromosome
    retain_mut=ifelse(shared_muts_df$Chrom=="Y"|(shared_muts_df$Chrom=="X" & shared_muts_df$Pos %in% c(60001:2699520,154931044:155260560)),FALSE,TRUE)
    filtered_muts$Genotype_shared_bin<-filtered_muts$Genotype_shared_bin[retain_mut,]
    filtered_muts$dna_strings=dna_strings_from_genotype(filtered_muts$Genotype_shared_bin)
  }
  
  #Update the genotype table and DNA strings to take account of LOH and del events that make the genotype unclear
  if(!is.null(CN_table)) {
    #Extract the chromosome and position info from the rows of the "Genotype_shared_bin" object used for tree building
    shared_muts_df=data.frame(mut_ref=rownames(filtered_muts$Genotype_shared_bin))
    shared_muts_df$Chrom=str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,1]
    shared_muts_df$Pos=as.integer(str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,2])
    
    for(k in 1:nrow(CN_table)){
      this_sample<-CN_table$Sample[k]
      if(CN_table$Type[k]%in%c("LOH","DEL")) {
        current_genotypes<-filtered_muts$Genotype_shared_bin[shared_muts_df$Chrom==CN_table$Chrom[k] &
                                                               shared_muts_df$Pos>CN_table$Pos_min[k]&
                                                               shared_muts_df$Pos<CN_table$Pos_max[k],this_sample]
        new_genotypes=pmax(0.5,current_genotypes)
        filtered_muts$Genotype_shared_bin[shared_muts_df$Chrom==CN_table$Chrom[k] &
                                            shared_muts_df$Pos>CN_table$Pos_min[k]&
                                            shared_muts_df$Pos<CN_table$Pos_max[k],this_sample] <-new_genotypes
      }
    }
    filtered_muts$dna_strings=dna_strings_from_genotype(filtered_muts$Genotype_shared_bin)
  }
  
  if(use_SNV_only) {
    filtered_muts$dna_strings=dna_strings_from_genotype(filtered_muts$Genotype_shared_bin[sapply(rownames(filtered_muts$Genotype_shared_bin),is.snv),])
  }
  
  #Decide an ID for this filtered set, depending on approach taken, and save
  save(filtered_muts, file = filtered_muts_file)
  
  #Write a fasta file of the dummy DNA strings - this can then be used by mpboot for treebuilding
  write.fasta(filtered_muts$dna_strings, names=names(filtered_muts$dna_strings), dna_string_file)
  
  ## RE-BUILD TREE with MPBoot ----
  system(paste0("mpboot -s ", dna_string_file," -bb 1000"))
  
  #Import the tree into R using ape
  tree <- read.tree(mpboot_tree_file)
  
  #Re-assign mutations to branches using treemut
  res=assign_mutations_to_branches(tree,filtered_muts,keep_ancestral=keep_ancestral,create_multi_tree=create_multi_tree,treefit_pval_cutoff=treefit_pval_cutoff)
  
  #Add node and pval information to the filtered_muts object
  print("Storing node and p-value information to the info dataframe")
  filtered_muts$COMB_mats.tree.build$mat$node <- res$tree$edge[res$summary$edge_ml,2]
  filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval
  
  #Make the res$tree output the new 'current working tree'
  tree<-res$tree
  
  #Assign edge lengths to the tree - safest to use this approach to assignment if have filtered any poor fit mutations
  print("Adjust tree edge lengths according to the new node information")
  tree$edge.length <- sapply(tree$edge[,2],function(node) {sum(filtered_muts$COMB_mats.tree.build$mat$node==node)})
}

print("Save the tree")
write.tree(tree, file = tree_file_path)

#========================================#
#ANNOTATION OF MUTATION SET USING VAGRENT ----
#========================================#

#Run VAGRENT on mutation set to annotate the filtered mutations
vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat)
write.table(vcf_file, sep = "\t", quote = FALSE, file = vcf_path, row.names = FALSE)

#1. paste vcf file to a dummy header file
system(paste0("cat ",vcf_header_path," ",vcf_path," > ", vagrent_input_path))

#2. commands to run vagrent
system(paste0("AnnotateVcf.pl -i ",vagrent_input_path," -o ",vagrent_output_path," -sp ",species," -as ",genome_build," -c ",vagrent_cache_path))

#3. import vagrent output and add into the filtered_muts object
vagrent_output = fread(vagrent_output_path,skip = "#CHROM")
annot_info = as.data.frame(str_split(vagrent_output$INFO, pattern = ";",simplify = TRUE), stringsAsFactors = FALSE)
colnames(annot_info) <- c("VT","VD","VC","VW")

annot_info$VC <- gsub(x=annot_info$VC, pattern = "VC=", replacement = "")
annot_info$VT <- gsub(x=annot_info$VT, pattern = "VT=", replacement = "")
annot_info$VW <- gsub(x=annot_info$VW, pattern = "VW=", replacement = "")
annot_info$VD <- gsub(x=annot_info$VD, pattern = "VD=", replacement = "")

filtered_muts$COMB_mats.tree.build$mat <- cbind(filtered_muts$COMB_mats.tree.build$mat,split_vagrent_output(df = annot_info,split_col = "VD"))

#Save the annotated filtered_muts files (post tree filtering)
save(filtered_muts, file = file_annot)

#========================================#
# SECOND STEP FOR FILTERING MIXED COLONIES----------
#========================================#

details=filtered_muts$COMB_mats.tree.build$mat
matrices=list(mtr=filtered_muts$COMB_mats.tree.build$NV,dep=filtered_muts$COMB_mats.tree.build$NR)

if(opt$n) {
  #Test branch VAFs for each sample to see if there are any that are inconsistent with a clonal sample
  sample_test=lapply(tree$tip.label,function(samples,min.depth=1) {
    print(samples)
    node_test=sapply(tree$edge[,2],function(node) {
      info=get_edge_info(tree,details,node)
      if(length(info$idx.in.details)>0) {
        df=data.frame(mtr=sum(matrices$mtr[info$idx,samples],na.rm = TRUE),
                      dep=sum(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
        df=df[which(df$dep>=min.depth),]
        df$vaf=df$mtr/df$dep
        df=df[which(!is.na(df$vaf)),]
        N=dim(df)[1]
        MTR=sum(df$mtr)
        DEP=sum(df$dep)
        min.mean.vaf=0.425
        if(DEP>0) {
          z1=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
          z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
          comb.p.value=max(z1$p.value,z2$p.value)
          if(comb.p.value<0.05){
            if(comb.p.value<0.05/dim(tree$edge)[1]){
              return(2)
            }else{
              return(1)
            }
          } else if(z1$p.value>0.05&z2$p.value<0.05){
            return(3)
          } else {
            return(0)
          }
        } else {
          return(NA)
        }
      } else {
        return(NA)
      }
    })
    names(node_test)<-tree$edge[,2]
    node_test=node_test[!is.na(node_test)]
    
    #Apply the 'consistent branch' test - checks that all the "POSITIVE" branches are those expected to be positive
    consistent_branches=as.numeric(names(node_test[node_test==3]))
    consistent_branch_lengths=sapply(consistent_branches,function(node) {sum(details$node==node)})
    consistent_branch_test=all(consistent_branches[consistent_branch_lengths>3]%in%get_ancestral_nodes(which(tree$tip.label==samples),tree$edge))
    
    #Combine into a single row data frame - at this stage, remove the terminal branch from the testing, this is allowed the have a few subclonal (potentially invitro) mutations that would result in the sample being inappropriately filtered as 'mixed'
    if(keep_ancestral){
      #If have an ancestral branch - ignore this branch for classifying mixed colonies. Mutations on this branch likely to be artefactual.
      root_branches<-tree$edge[which(tree$edge[,1]==tree$edge[1,1]),2]
      branches_to_ignore=c(root_branches,which(tree$tip.label==samples))
    } else {
      branches_to_ignore=which(tree$tip.label==samples)
    }
    df=data.frame(sample=samples,high_prob=sum(node_test==2 & !as.numeric(names(node_test))%in%branches_to_ignore),low_prob=sum(node_test==1& !as.numeric(names(node_test))%in%branches_to_ignore),consistent_pos_branches=consistent_branch_test)
    return(df)
  })
  sample_test_df=Reduce(rbind,sample_test)
  
  #Define which samples have "passed" this test of clonality
  sample_test_df$result=ifelse(sample_test_df$high_prob>0|sample_test_df$low_prob>2|!sample_test_df$consistent_pos_branches,"FAIL","PASS")
  #sample_test_df$result=ifelse(sample_test_df$high_prob>0|sample_test_df$low_prob>2,"FAIL","PASS")
  mixed_samples=as.character(sample_test_df$sample[sample_test_df$result=="FAIL"])
  if(length(mixed_samples)==0) {
    print("No samples failed the test for clonality")
  } else {
    print(paste(mixed_samples,"will be removed as failed the test for clonality"))
    #Plot the vaf trees of the failing samples
    pdf(paste0(output_dir,"/plots/",Run_ID,"_",filtering_ID,"_mixed_vaf_trees.pdf"),width=30,height = 10)
    temp=lapply(mixed_samples,function(sample) {
      tree=plot_tree(tree,cex.label=0)
      temp=plot_tree_vaf(tree,
                         details=filtered_muts$COMB_mats.tree.build$mat,
                         matrices=list(mtr=filtered_muts$COMB_mats.tree.build$NV,dep=filtered_muts$COMB_mats.tree.build$NR),
                         samples=sample)
    })
    dev.off()
    
    #Remove the mixed colony samples and their private mutations
    filtered_muts$COMB_mats.tree.build$Genotype_bin<-filtered_muts$COMB_mats.tree.build$Genotype_bin[,!colnames(filtered_muts$COMB_mats.tree.build$Genotype_bin)%in%mixed_samples]
    filtered_muts$COMB_mats.tree.build$NV<-filtered_muts$COMB_mats.tree.build$NV[,!colnames(filtered_muts$COMB_mats.tree.build$NV)%in%mixed_samples]
    filtered_muts$COMB_mats.tree.build$NR<-filtered_muts$COMB_mats.tree.build$NR[,!colnames(filtered_muts$COMB_mats.tree.build$NR)%in%mixed_samples]
    filtered_muts$COMB_mats.tree.build$PVal<-filtered_muts$COMB_mats.tree.build$PVal[,!colnames(filtered_muts$COMB_mats.tree.build$PVal)%in%mixed_samples]
    still_positive=rowSums(filtered_muts$COMB_mats.tree.build$Genotype_bin==1) > 0
    filtered_muts$COMB_mats.tree.build=list_subset(filtered_muts$COMB_mats.tree.build,select_vector=still_positive)
    
    #Decide an ID for this filtered set, depending on approach taken, and save
    save(filtered_muts, file = file_annot_post_mix)
    
    #Get the new "genotype_shared_bin" for tree building
    Genotype_shared_bin=filtered_muts$COMB_mats.tree.build$Genotype_bin[rowSums(filtered_muts$COMB_mats.tree.build$Genotype_bin == 1) > 1,]
    
    if(filter_y) {
      print("Removing variants on the Y chromosome and X chromosome PARs from the dummy dna strings used for tree building")
      #Extract the chromosome and position info from the rows of the "Genotype_shared_bin" object used for tree building
      shared_muts_df=data.frame(mut_ref=rownames(Genotype_shared_bin))
      shared_muts_df$Chrom=str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,1]
      shared_muts_df$Pos=str_split(shared_muts_df$mut_ref,pattern = "-",simplify = T)[,2]
      
      #Remove the variant from the genotype_shared_bin if it is on the Y chromosome or the PAR of the X chromosome
      retain_mut=ifelse(shared_muts_df$Chrom=="Y"|(shared_muts_df$Chrom=="X" & shared_muts_df$Pos %in% c(60001:2699520,154931044:155260560)),FALSE,TRUE)
      filtered_muts$dna_strings=dna_strings_from_genotype(Genotype_shared_bin[retain_mut,])
    } else {
      filtered_muts$dna_strings=dna_strings_from_genotype(Genotype_shared_bin)
    }
    
    #Write a fasta file of the dummy DNA strings
    write.fasta(filtered_muts$dna_strings, names=names(filtered_muts$dna_strings), dna_string_file_post_mix)
    
    ## RE-BUILD TREE with MPBoot ----
    print("Rerunning tree-building now that mixed samples have been removed")
    system(paste0("mpboot -s ", dna_string_file_post_mix," -bb 1000"))
    
    #Import the tree into R using ape
    tree <- read.tree(mpboot_tree_file_post_mix)
    
    ## Re-assign mutations using treemut ----
    res=assign_mutations_to_branches(tree,filtered_muts,keep_ancestral=keep_ancestral,create_multi_tree=create_multi_tree,treefit_pval_cutoff=treefit_pval_cutoff)
    
    #Add node and pval information to the filtered_muts object
    filtered_muts$COMB_mats.tree.build$mat$node <- res$tree$edge[res$summary$edge_ml,2]
    filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval
    tree<-res$tree
    
    #Assign edge lengths to the tree - safest to use this approach to assignment if have filtered any poor fit mutations
    tree$edge.length <- sapply(tree$edge[,2],function(node) {sum(filtered_muts$COMB_mats.tree.build$mat$node==node)})
  }
  
  #Save the res file and the tree file
  write.tree(tree, file = tree_file_path_post_mix)
  
  #Re-save the filtered_muts file with the new node numbers and pvals
  save(filtered_muts, file = file_annot_post_mix)
}

#========================================#
# CORRECT THE EDGE LENGTHS BASED ON SAMPLE SENSITIVITY ----
#========================================#

#NB - for function to work correctly, the sensitivity_df must be set up in exactly the correct format
#If there are no SNVs in the tree, and/or sensitivity analysis for these, set "include_indels" to FALSE
#NOTE ON APPROACH:
#CORRECTION is done BEFORE removing duplicates
#Duplicate samples are then removed, and mutations reassigned to the new tree without changing branch lengths

#Create trees of only SNVs or only INDELs (this function does not do any correction)
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")

#Import the sensitivity analysis file
sensitivity_analysis_df <- read.delim(sensitivity_df_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")
sensitivity_analysis_df<-rbind(sensitivity_analysis_df,data.frame(Sample="Ancestral",SNV_sensitivity=1,INDEL_sensitivity=1))
sensitivity_analysis_df$SNV_sensitivity[sensitivity_analysis_df$SNV_sensitivity==0]<-0.001
sensitivity_analysis_df$INDEL_sensitivity[sensitivity_analysis_df$INDEL_sensitivity==0]<-0.001

#Create corrected trees
tree_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = TRUE, sensitivity_df = sensitivity_analysis_df)
tree_SNV_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_analysis_df)
tree_INDEL_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_SNVs = FALSE, sensitivity_df = sensitivity_analysis_df)

#Visualize the different trees
pdf(paste0(output_dir,"/plots/",Run_ID,"_",filtering_ID,"_summary_trees_with_duplicates.pdf"),width=30,height = 10)
par(mfrow = c(2,2))
plot(tree, show.tip.label = FALSE, main = "Uncorrected tree: SNVs and INDELs",direction="downwards",edge.color=highlight_groups(tree,group1 = grep(donor_ID,tree$tip.label,value = T),group2 = grep(recip_ID,tree$tip.label,value = T)))
plot(tree_c, show.tip.label = FALSE, main = "Corrected tree: SNVs and INDELs",direction="downwards",edge.color=highlight_groups(tree,group1 = grep(donor_ID,tree$tip.label,value = T),group2 = grep(recip_ID,tree$tip.label,value = T)))
plot(tree_SNV_c, show.tip.label = FALSE, main = "Corrected tree: SNVs only",direction="downwards",edge.color=highlight_groups(tree,group1 = grep(donor_ID,tree$tip.label,value = T),group2 = grep(recip_ID,tree$tip.label,value = T)))
plot(tree_INDEL_c, show.tip.label = FALSE, main = "Corrected tree: INDELs only",direction="downwards",edge.color=highlight_groups(tree,group1 = grep(donor_ID,tree$tip.label,value = T),group2 = grep(recip_ID,tree$tip.label,value = T)))
par(mfrow = c(1,1))
dev.off()

#========================================#
#DETECT DUPLICATE SAMPLES AND REMOVE ----
#========================================#

mean_mut_burden=mean(get_mut_burden(tree))
private_mut_threshold=max(30,0.03*mean_mut_burden)
print(paste("Threshold for determining duplicate samples is",private_mut_threshold,"mutations."))
pseudo_terminal_nodes=sapply(tree$edge[,2][!tree$edge[,2]%in%1:length(tree$tip.label)],function(node) {
  node_height=nodeheight(tree = tree,node=node)
  samples=get_edge_info(tree,filtered_muts$COMB_mats.tree.build$mat,node)$samples
  sample_heights=nodeHeights(tree)[tree$edge[,2]%in%which(tree$tip.label%in%samples),2]
  
  if(all((sample_heights-node_height)<private_mut_threshold)){ #This is the method of determining the duplicates - may need to alter the "30" for low mutation burden samples
    return(node)
  }else{
    return(NA)
  }
})

pseudo_terminal_nodes=pseudo_terminal_nodes[!is.na(pseudo_terminal_nodes)]
duplicate_samples=lapply(pseudo_terminal_nodes,function(node) getTips(tree,node=node))
lapply(duplicate_samples,function(x) print(paste(paste(x,collapse=" & "),"are recognised as duplicate samples")))

#Choose which of the duplicate samples to drop
drop_samples_list=lapply(duplicate_samples,function(samples) {
  if(sum(samples%in%tree$tip.label)>1) {
    included_samples=samples[samples%in%tree$tip.label]
    sample_heights=sapply(included_samples,function(sample) {nodeheight(tree = tree,node=which(tree$tip.label==sample))})
    retain_sample=included_samples[sample_heights==min(sample_heights)][1]
    return(included_samples[!included_samples==retain_sample])
  }else{
    return(NULL)
  }
})
drop_samples=unlist(drop_samples_list)

if(length(drop_samples)==0) {
  print("No duplicate samples identified")
} else {
  print(paste("Dropping",drop_samples,"as one of a duplicate set"))
  
  #Now remove these samples and any mutations that are only called "positive" in that sample
  filtered_muts$COMB_mats.tree.build$Genotype_bin<-filtered_muts$COMB_mats.tree.build$Genotype_bin[,!colnames(filtered_muts$COMB_mats.tree.build$Genotype_bin)%in%drop_samples]
  filtered_muts$COMB_mats.tree.build$NV<-filtered_muts$COMB_mats.tree.build$NV[,!colnames(filtered_muts$COMB_mats.tree.build$NV)%in%drop_samples]
  filtered_muts$COMB_mats.tree.build$NR<-filtered_muts$COMB_mats.tree.build$NR[,!colnames(filtered_muts$COMB_mats.tree.build$NR)%in%drop_samples]
  filtered_muts$COMB_mats.tree.build$PVal<-filtered_muts$COMB_mats.tree.build$PVal[,!colnames(filtered_muts$COMB_mats.tree.build$PVal)%in%drop_samples]
  still_positive=rowSums(filtered_muts$COMB_mats.tree.build$Genotype_bin==1) > 0
  filtered_muts$COMB_mats.tree.build=list_subset(filtered_muts$COMB_mats.tree.build,select_vector=still_positive)
  
  #Now import the tree from mpboot and drop the duplicate samples - do this to the version that is already corrected
  tree_SNV_c <- drop.tip(tree_SNV_c,drop_samples)
  
  #Assign mutations back to the tree
  if(!keep_ancestral) {
    #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
    df = reconstruct_genotype_summary(tree_SNV_c) #Define df (data frame) for treeshape
    
    #Get matrices in order, and run the main assignment functions
    mtr = filtered_muts$COMB_mats.tree.build$NV; mtr = as.matrix(mtr)
    depth = filtered_muts$COMB_mats.tree.build$NR; depth = as.matrix(depth)
    p.error = sapply(df$samples,function(x) ifelse(x=="Ancestral",1e-6,0.01))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
    
  } else {
    #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
    df = reconstruct_genotype_summary(tree_SNV_c) #Define df (data frame) for treeshape
    
    #Get matrices in order, and run the main assignment functions
    mtr = filtered_muts$COMB_mats.tree.build$NV; mtr$Ancestral=0;mtr = as.matrix(mtr)
    depth = filtered_muts$COMB_mats.tree.build$NR; depth$Ancestral=10; depth = as.matrix(depth)
    p.error = sapply(df$samples,function(x) ifelse(x=="Ancestral",1e-6,0.01))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  }
  
  #Get res (results!) object
  filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval
  filtered_muts$COMB_mats.tree.build$mat$node <- tree_SNV_c$edge[res$summary$edge_ml,2]
}

write.tree(tree_SNV_c,file=tree_file_path_post_mix_post_dup)

#========================================#
# DEFINE & ANNOTATE MUTATIONS OF INTEREST ----
#========================================#

filtered_muts$COMB_mats.tree.build$mat$coding_change <- ifelse(filtered_muts$COMB_mats.tree.build$mat$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                                                                                  "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                                                                                  "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                                                                                  "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                                                                                  "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                                                                                  "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss")|
                                                                 grepl("splice_site_variant",filtered_muts$COMB_mats.tree.build$mat$Type),
                                                               "Coding change",
                                                               "no")
filtered_muts$COMB_mats.tree.build$mat$coding_change_chip = ifelse(filtered_muts$COMB_mats.tree.build$mat$Gene %in% genes_to_annotate & filtered_muts$COMB_mats.tree.build$mat$coding_change == "Coding change",
                                                                   "Coding change mutation in driver", "no")
filtered_muts$COMB_mats.tree.build$mat$ChromPos=paste(filtered_muts$COMB_mats.tree.build$mat$Chrom,filtered_muts$COMB_mats.tree.build$mat$Pos,sep="-")
filtered_muts$COMB_mats.tree.build$mat$variant_ID=paste(filtered_muts$COMB_mats.tree.build$mat$Gene, filtered_muts$COMB_mats.tree.build$mat$Protein, sep = " ")

save(filtered_muts, file = file_annot_post_mix_post_dup)

#========================================#
#FINAL SUMMARY PLOTS INCLUDING NON-SYNONYMOUS MUTATIONS IN CH GENES ----
#========================================#

if(is.na(donor_age)){donor_age<-mean(get_mut_burden(tree_SNV_c))} #If no donor age is provided, scale to the mean mutation burden across colonies

#Now get the donor & recipient trees
tree_SNV_c$D_or_R =ifelse(grepl(donor_ID,tree_SNV_c$tip.label),"D",ifelse(grepl("Ancestral",tree_SNV_c$tip.label),"Ancestral","R"))
tree_SNV_c$coords<-NULL #reset the coords so that the plot_tree function will rework them out

#Make ultrametric tree of the combined tree
tree_SNV_c_ultra <- make.ultrametric.tree(tree_SNV_c)
tree_SNV_c_ultra$edge.length[tree_SNV_c_ultra$edge[,2]==which(tree_SNV_c_ultra$tip.label=="Ancestral")]<-0
tree_SNV_c_ultra$edge.length<-tree_SNV_c_ultra$edge.length*donor_age

if(!all(tree_SNV_c$D_or_R%in%c("D","Ancestral"))) {
  tree_SNV_c.R=drop.tip(tree_SNV_c,tree_SNV_c$tip.label[tree_SNV_c$D_or_R=="D"])
  tree_SNV_c.D=drop.tip(tree_SNV_c,tree_SNV_c$tip.label[tree_SNV_c$D_or_R=="R"])
  
  #Make ultra-metric versions of the combined, donor and recipient trees
  tree_SNV_c_ultra.R <- make.ultrametric.tree(tree_SNV_c.R);tree_SNV_c_ultra.R$edge.length <- tree_SNV_c_ultra.R$edge.length*donor_age
  tree_SNV_c_ultra.R$edge.length[tree_SNV_c_ultra.R$edge[,2]==which(tree_SNV_c_ultra.R$tip.label=="Ancestral")]<-0
  tree_SNV_c_ultra.D <- make.ultrametric.tree(tree_SNV_c.D);tree_SNV_c_ultra.D$edge.length <- tree_SNV_c_ultra.D$edge.length*donor_age
  tree_SNV_c_ultra.D$edge.length[tree_SNV_c_ultra.D$edge[,2]==which(tree_SNV_c_ultra.D$tip.label=="Ancestral")]<-0
}

details=filtered_muts$COMB_mats.tree.build$mat

#PLOT THEM (the all in one)
pdf(file = paste0(output_dir,"/plots/",Run_ID,"_",filtering_ID,"_final_summary_trees.pdf"),width=15,height=7)

tree_SNV_c=plot_tree(tree_SNV_c,cex.label=0)
temp=add_annotation(tree_SNV_c,
                    details,
                    matrices,
                    annot_function=plot_sharing_info,
                    donor_ID=donor_ID,
                    recip_ID=recip_ID,
                    sharing_cols=c("black", "#11a0aa80", "#c8256580"))
temp=plot_tree_labels(tree_SNV_c,
                      details = details,
                      type="line",
                      query.field = "coding_change_chip",
                      data.frame(value="Coding change mutation in driver",col="red",pch = 17,stringsAsFactors = FALSE),
                      label.field = "variant_ID",
                      cex.label = 0.7,
                      lty=2)

tree_SNV_c_ultra=plot_tree(tree_SNV_c_ultra,cex.label=0)
temp=add_annotation(tree_SNV_c_ultra,
                    details,
                    matrices,
                    annot_function=plot_sharing_info,
                    donor_ID=donor_ID,
                    recip_ID=recip_ID,
                    sharing_cols=c("black", "#11a0aa80", "#c8256580"))
temp=plot_tree_labels(tree_SNV_c_ultra,
                      details = details,
                      type="line",
                      query.field = "coding_change_chip",
                      data.frame(value="Coding change mutation in driver",col="red",pch = 17,stringsAsFactors = FALSE),
                      label.field = "variant_ID",
                      cex.label = 0.7,
                      lty=2)

if(!all(tree_SNV_c$D_or_R%in%c("D","Ancestral"))) {
  tree_SNV_c.R=plot_tree(tree_SNV_c.R,cex.label=0)
  temp=add_annotation(tree_SNV_c.R,
                      details,
                      matrices,
                      annot_function=plot_sharing_info,
                      donor_ID=donor_ID,
                      recip_ID=recip_ID,
                      sharing_cols=c("black", "#11a0aa80", "#c8256580"))
  
  tree_SNV_c.D=plot_tree(tree_SNV_c.D,cex.label=0)
  temp=add_annotation(tree_SNV_c.D,
                      details,
                      matrices,
                      annot_function=plot_sharing_info,
                      donor_ID=donor_ID,
                      recip_ID=recip_ID,
                      sharing_cols=c("black", "#11a0aa80", "#c8256580"))
  
  tree_SNV_c_ultra.R=plot_tree(tree_SNV_c_ultra.R,cex.label=0)
  temp=add_annotation(tree_SNV_c_ultra.R,
                      details,
                      matrices,
                      annot_function=plot_sharing_info,
                      donor_ID=donor_ID,
                      recip_ID=recip_ID,
                      sharing_cols=c("black", "#11a0aa80", "#c8256580"))
  
  tree_SNV_c_ultra.D=plot_tree(tree_SNV_c_ultra.D,cex.label=0)
  temp=add_annotation(tree_SNV_c_ultra.D,
                      details,
                      matrices,
                      annot_function=plot_sharing_info,
                      donor_ID=donor_ID,
                      recip_ID=recip_ID,
                      sharing_cols=c("black", "#11a0aa80", "#c8256580"))
  
}

#Plot all coding change mutations shared by >1 sample
details$coding_change_shared=sapply(1:nrow(details),function(i) {
  if(length(getTips(tree_SNV_c,details$node[i]))>1 &
     details$coding_change[i] == "Coding change") {
    return("yes")
  } else {
    return("no")
  }
})

tree_SNV_c_ultra=plot_tree(tree_SNV_c_ultra, cex.label = 0)
temp=add_annotation(tree_SNV_c_ultra,
                    details,
                    matrices,
                    annot_function=plot_sharing_info,
                    donor_ID=donor_ID,
                    recip_ID=recip_ID,
                    sharing_cols=c("black", "#11a0aa80", "#c8256580"))
temp=plot_tree_labels(tree_SNV_c_ultra,
                      details = details,
                      type = "label",
                      query.field = "coding_change_shared",
                      data.frame(value="yes",col="red",lwd = 1,pch=1,stringsAsFactors = FALSE),
                      label.field = "variant_ID",
                      cex.label = 0.5)
dev.off()



