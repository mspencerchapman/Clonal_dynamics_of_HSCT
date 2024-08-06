library(optparse)

option_list = list(
  make_option(c("-r", "--run_ID"), action="store", type='character', help="Run ID for the filtered params output to be filtered"),
  make_option(c("-b", "--bed"), action="store", type='character', help="path to the MS filtered bedfile"),
  make_option(c("-d", "--dir"), action="store", default="SNVs", type='character', help="Directory for the filtered file")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

#Script for selecting out the SNVs that have passed the MS filters (Mathijs's filters) from the previously built COMB_mats sets.
#Removes the need to re-run cgpVAF on the reduced set
#Leaves in all INDELs
source("/lustre/scratch126/casm/team154pc/ms56/my_functions/foetal.filters.parallel.R")

Run_ID = opt$r
MS_bed_path = opt$b
output_directory = opt$d

#Load up the pre-MS COMB_mats file
mats_and_params_file = paste0(output_directory,"/mats_and_params_", Run_ID) #This script assumes this file naming format
load(mats_and_params_file)

#Read in the MS bed file, strip white space from the Pos column, and convert into a mut_ref column
MS <- read.delim(MS_bed_path, sep = "\t",header = FALSE, stringsAsFactors=FALSE, strip.white = TRUE)
colnames(MS) <- c("Chrom","Pos","Ref","Alt")
MS$Pos <- sapply(MS$Pos, function(x) gsub("^\\s+|\\s+$","",x))
MS$mut_ref <- apply(MS,1,paste,collapse = "-")

#Check that the (almost) all of the MS pass muts are in the COMB_mats mut_ref column
sum(MS$mut_ref %in% COMB_mats$mat$mut_ref)/ length(MS$mut_ref)
#MS$mut_ref[!MS$mut_ref %in% COMB_mats$mat$mut_ref]

#Filter the COMB_mats matrices to take out SNVs that didn't make it through the MS filters
MS_select = COMB_mats$mat$mut_ref%in%MS$mut_ref|COMB_mats$mat$Mut_type == "INDEL"
COMB_mats <- list_subset(list = COMB_mats,select_vector = MS_select)
filter_params <- filter_params[MS_select,]
save(COMB_mats, filter_params, file = paste0(output_directory,"/mats_and_params_", Run_ID,"_postMS"))

#Now save the "reduced" set without the definite germline mutations
reduced = filter_params$germline_pval < 0.01|filter_params$bb_rhoval >0.01
filter_params = filter_params[reduced,]
COMB_mats = list_subset(COMB_mats, reduced)
save(COMB_mats, filter_params, file = paste0(output_directory, "/mats_and_params_", Run_ID,"_postMS_reduced"))
