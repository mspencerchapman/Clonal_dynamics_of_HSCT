#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("dplyr","data.table")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

#========================================#
# Define fixed paths ####
#========================================#

root_dir="~/R_work/Clonal_dynamics_of_HSCT/"

metadata_file_path=paste0(root_dir,"/data/metadata_temp.tsv")
output_metadata_file_path=paste0(root_dir,"/data/metadata_temp_updated.tsv")
sensitivity_mat_path=paste0(root_dir,"/01 Generating_the_phylogenies/sensitivity_matrix.csv")
setwd(root_dir)

#========================================#
# Define custom functions ####
#========================================#

parse_MSfilter_output=function(file_path,for_bed_file=T){
  if(for_bed_file){
    MS=data.table::fread(file_path,header = F)[,c(1,2,4,5)]
  }else{
    MS=data.table::fread(file_path,skip='#CHROM',header = T)[,c(1,2,4,5)]
  }
  colnames(MS) <- c("Chrom","Pos","Ref","Alt")
  MS$Pos <- sapply(MS$Pos, function(x) gsub("^\\s+|\\s+$","",x))
  MS$mut_ref <- apply(MS,1,paste,collapse = "-")
  return(MS)
}

create_sample_mat=function(Sample,COMB_mats,germline_SNVs,MS_dir,MS_suffix="_complete_final_retained_3.vcf_for_bed_file"){
  sample_NV=COMB_mats$NV[germline_SNVs,Sample]
  sample_NR=COMB_mats$NR[germline_SNVs,Sample]
  called_muts<-parse_MSfilter_output(paste0(MS_dir,'/',Sample,MS_suffix))
  called<-germline_SNVs%in%called_muts$mut_ref
  sample_mat=data.frame(mut_ref=germline_SNVs,NV=sample_NV,NR=sample_NR,VAF=sample_NV/sample_NR,called=called)
  return(sample_mat)
}


# #---------CREATE SENSITIVITY MATRIX FOR ALL COMBINATIONS OF variant read/ depth NUMBERS----------------------
# #Define sensitivity for given NV and NR from the germline SNPs
# #Could use any samples from any Individual - but here uses first 20 samples from BCL002
# 
# ID="BCL002"
# mats_and_params_file=paste0("filtering_runs/mats_and_params/mats_and_params_",ID,"_m40")
# load(mats_and_params_file) #Path to full output file of the "Mutation_filtering_get_parameters.R" script
# 
# binom_cutoff=6
# rho_cutoff=0.05
# depth_cutoff=500
# 
# germline_SNVs = COMB_mats$mat$mut_ref[COMB_mats$mat$Mut_type == "SNV" & log10(filter_params$germline_pval) > -binom_cutoff & filter_params$bb_rhoval < rho_cutoff & rowSums(COMB_mats$NR) > depth_cutoff]
# test_samples=colnames(COMB_mats$NV)[1:20]
# 
# list_of_sample_mats<-lapply(test_samples,function(Sample){
#   print(Sample)
#   sample_mat<-create_sample_mat(Sample=Sample,COMB_mats=COMB_mats,germline_SNVs=germline_SNVs,MS_dir = 'BCL002/MS_filters/output_files')
# })
# names(list_of_sample_mats)<-test_samples
# 
# all_sample_mats<-dplyr::bind_rows(list_of_sample_mats)
# 
# sens_mat=matrix(NA,nrow = 34,ncol = 34,dimnames=list(NV=1:34,NR=1:34))
# n_samp_mat=matrix(NA,nrow = 34,ncol = 34,dimnames=list(NV=1:34,NR=1:34))
# 
# for(i in 1:nrow(sens_mat)){
#   print(i)
#   for(j in 1:ncol(sens_mat)){
#     n_samp=sum(all_sample_mats$NV==i&all_sample_mats$NR==j)
#     n_samp_mat[i,j]<-n_samp
#     if(n_samp>0) {
#       sens=sum(all_sample_mats$called[all_sample_mats$NV==i&all_sample_mats$NR==j])/sum(all_sample_mats$NV==i&all_sample_mats$NR==j)
#       sens_mat[i,j]<-sens
#     }
#   }
# }
# 
# sens_mat[n_samp_mat<10]<-NA
# 
# #If only 1-2 variant reads, sensitivity will always be 0
# sens_mat[1:2,]<-0
# 
# write.table(sens_mat,file = sensitivity_mat_path,quote = F,sep = ",",row.names = T,col.names = T)

#========================================#
# Use the SENSITIVITY MATRIX as new reference for defining 'SOMATIC SNV SENSITIVITY' ####
#========================================#

#THEN  USE THE 'SENSITIVITY MATRIX' REFERENCE FOR CALCULATING NEW SENSITIVITY INCLUDING THE CLONALITY DATA)
sample_metadata<-read.delim(metadata_file_path,sep="\t")
sens_mat=read.csv(sensitivity_mat_path)
rownames(sens_mat)<-colnames(sens_mat)<-1:ncol(sens_mat)

#Simulate sequencing data
sample_metadata=dplyr::bind_rows(lapply(unique(sample_metadata$ID),function(Individual) {
  cat(paste0(Individual,"\n"))
  
  #Get the subsetted metadata for that individual
  sample_metadata_Ind<-sample_metadata%>%dplyr::filter(ID==Individual)
  
  #Define the germline SNVs
  mats_and_params_file=paste0("filtering_runs/mats_and_params/mats_and_params_",Individual,"_m40")
  load(mats_and_params_file) #Path to full output file of the "Mutation_filtering_get_parameters.R" script
  binom_cutoff=6; rho_cutoff=0.05; depth_cutoff=500
  germline_SNVs = COMB_mats$mat$mut_ref[COMB_mats$mat$Mut_type == "SNV" & log10(filter_params$germline_pval) > -binom_cutoff & filter_params$bb_rhoval < rho_cutoff & rowSums(COMB_mats$NR) > depth_cutoff]
  
  #Import the baseline germline sensitivity data from pre-existing file
  sample_germline_sens<-read.delim(paste0(Individual,"/sensitivity_analysis_",Individual),sep=" ")
  colnames(sample_germline_sens)<-c("Sample","germline_SNV_sensitivity","germline_INDEL_sensitivity")
  sample_metadata_Ind<-left_join(sample_metadata_Ind,sample_germline_sens)
  
  #Now use the peak VAF and sensitivity measures to calculate a 'somatic SNV sensitivity'
  sample_metadata_Ind$somatic_SNV_sensitivity<-sapply(1:nrow(sample_metadata_Ind),function(i) {
    Sample=sample_metadata_Ind$Sample[i]
    peak_vaf=min(sample_metadata_Ind$peak_vaf[i],0.5)
    germline_sens<-sample_metadata_Ind$germline_SNV_sensitivity[i]
    
    cat(Sample,paste0("Peak vaf: ",peak_vaf),paste0("Sensitivity for germline SNPs: ",germline_sens),sep="\n")
    
    if(is.na(peak_vaf)){
      stop(return(NA))
    } else if(abs(peak_vaf-0.5)<0.02){ #if clonality is >0.48, return the germline sensitivity - otherwise proceed to calculate 'clonality corrected' sensitivity
        stop(return(germline_sens))
    }
    
    sim_NR=sample(x=COMB_mats$NR[germline_SNVs,Sample],size=1e4,replace=T)
    sim_NR[sim_NR>34]<-34 #sensitivity matrix only covers depths up to 34 -> therefore set max depth as 34
    sim_NR[sim_NR>ncol(sens_mat)]<-ncol(sens_mat)
    sim_NV=rbinom(n=length(sim_NR),size=sim_NR,prob=peak_vaf) #binomial draw from the depth data, with clonality
    sim_mat=data.frame(NV=sim_NV,NR=sim_NR,VAF=sim_NV/sim_NR)
    
    sim_mat$called=sapply(1:nrow(sim_mat),function(i) {
      if(sim_mat$NV[i]==0){stop(return(F))}
      sens=sens_mat[sim_mat$NV[i],sim_mat$NR[i]]
      if(is.na(sens)){
        return(NA)
      } else {
        return(runif(1)<sens)
      }
    })
    
    #Any mutation with a VAF < 0.2 will subsequently be filtered out
    sim_mat$called[sim_mat$VAF<0.2]<-F
    cat(paste0(sum(is.na(sim_mat$called))," results are out of range of the sensitivity matrix"),sep="\n")
    overall_sensitivity=sum(sim_mat$called[!is.na(sim_mat$called)])/sum(!is.na(sim_mat$called))
    cat(paste0("Inferred somatic SNV sensitivity: ",overall_sensitivity),sep="\n")
    return(overall_sensitivity)
  })
  
  return(sample_metadata_Ind)
}))

#========================================#
# SAVE THE UDPATED METADATA with SOMATIC SNV SENSITIVITY ####
#========================================#
cat("Saving updated metadata")
write.table(sample_metadata,file=output_metadata_file_path,quote=F,sep="\t",row.names=F)
