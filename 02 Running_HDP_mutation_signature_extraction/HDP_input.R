#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
ccran_packages=c("dplyr","stringr","seqinr","MASS")
bioconductor_packages=c("GenomicRanges","Rsamtools")

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
# Set file path for human genome ####
#========================================#
local_ref="~/R_work/Reference_files/genome.fa"
lustre_ref="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

#========================================#
# Define function for extracting trinucleotide counts ####
#========================================#
mutlist_to_96_contexts=function(mutlist,genomeFile=lustre_ref){
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  samples=unique(mutlist$SampleID)
  trinuc_mut_mat=matrix(0,ncol=96,nrow=length(samples))
  for (n in 1:length(samples)){
    s=samples[n]
    mutations=as.data.frame(mutlist[mutlist$SampleID==s,c("Chr","Pos","Ref","Alt")])
    colnames(mutations) = c("chr","pos","ref","mut")
    mutations$pos=as.numeric(mutations$pos)
    mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                       as.numeric(mutations$pos)+1))))
    ntcomp = c(T="A",G="C",C="G",A="T")
    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
    mutations$trinuc_ref_py = mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A","G")) { # Purine base
        mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
        mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
      }
    }
    freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    trinuc_mut_mat[n,]=freqs_full
    print(s)
  }
  colnames(trinuc_mut_mat)=full_vec
  rownames(trinuc_mut_mat)=samples
  return(trinuc_mut_mat)
}

#========================================#
# Set the root directory and read in the necessary files ####
#========================================#

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
source(paste0(root_dir,"/data/HSCT_functions.R"))
HDP_folder=paste0(root_dir,"/data/HDP")
annotated_muts_folder=paste0(root_dir,"/data/annot_files/")

#Import the annotated muts files
annotated_muts_paths=list.files(annotated_muts_folder,pattern="annotated_mut_set",full.names = T)
all.muts<-lapply(exp_nos,function(exp_no){load(grep(paste0("Pair",exp_no),annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
names(all.muts)<-paste0("Pair",exp_nos)

#========================================#
# Bind the mutations into a single dataframe ####
#========================================#

all_muts=Map(details=all.muts,patient=names(all.muts),function(details,patient) {
  details=details[details$Mut_type=="SNV",]
  details=details[,c("Chrom","Pos","Ref","Alt","node")]
  details$node=paste(details$node,patient,sep = "_")
  colnames(details)=c("Chr","Pos","Ref","Alt","SampleID")
  return(details)
})%>%dplyr::bind_rows()

dim(all_muts)

#========================================#
# Apply the function to extract trinuc contexts ####
#========================================#
trinuc_mut_mat=mutlist_to_96_contexts(all_muts,genomeFile = local_ref)

#========================================#
# Save output ####
#========================================#
samples=rownames(trinuc_mut_mat)
key_table=data.frame(Sample=samples,Patient=str_split(samples,pattern="_",simplify=TRUE)[,2])
write.table(trinuc_mut_mat,file = paste0(HDP_folder,"/trinuc_mut_mat.txt"))
write.table(key_table,file=paste0(HDP_folder,"/key_table.txt"))
