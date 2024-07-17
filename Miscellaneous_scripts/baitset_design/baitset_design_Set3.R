library(stringr)
library(ape)
library(seqinr)
library(data.table)
library(tidyr)
library(readr)
library(dplyr)

my_working_directory = "/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/"
my_functions_dir="/lustre/scratch119/casm/team154pc/ms56/my_functions"
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"

setwd(my_working_directory)
R_scripts=list.files(my_functions_dir,pattern = ".R",full.names = T)
sapply(R_scripts[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory) #R scripts have to be sourced within containing directory or outputs an error.

#Functions for this script
write_target_batches=function(details,sampleID,coord_buffer=60,batch_size=36275,output_dir) {
  print(paste("Writing mutation coordinates in",ceiling(length(details$mut_ref)/batch_size),"batches of",batch_size))
  for(i in 1:ceiling(length(details$mut_ref)/batch_size)) {
    start_idx =batch_size*(i-1) + 1; end_idx=batch_size*i
    mut_batch=details$mut_ref[start_idx:end_idx]
    mut_batch_coords <- as.data.frame(str_split(mut_batch, pattern = "-", simplify = TRUE)[,1:2], stringsAsFactors = FALSE)
    colnames(mut_batch_coords) = c("Chr", "Pos")
    mut_batch_coords$Chr <- paste0("chr", mut_batch_coords$Chr)
    mut_batch_coords$Pos <- as.numeric(mut_batch_coords$Pos)
    mut_batch_coords$Start <- as.character(mut_batch_coords$Pos - coord_buffer - 1)
    mut_batch_coords$End <- as.character(mut_batch_coords$Pos + coord_buffer)
    mut_batch_baitset <- mut_batch_coords[,-2] #Remove the original Pos column
    write_delim(mut_batch_baitset, file = paste0(output_dir,"/",sampleID,"_batch_muts_baits_",i,".txt"), delim = " ", col_names = FALSE)
  }
}

write_mut_coords_bedfile=function(details,all_mut_coords_path) {
  require(stringr)
  require(readr)
  mut_batch_coords <- as.data.frame(str_split(details$mut_ref, pattern = "-", simplify = TRUE)[,1:2], stringsAsFactors = FALSE)
  colnames(mut_batch_coords) = c("Chr", "Pos")
  mut_batch_coords$Chr <- paste0("chr", mut_batch_coords$Chr)
  mut_batch_coords$Pos <- as.numeric(mut_batch_coords$Pos)
  mut_batch_coords$Start <- mut_batch_coords$Pos
  mut_batch_coords$End <- mut_batch_coords$Pos
  mut_batch_coords <- mut_batch_coords[,-2] #Remove the original Pos column
  write_delim(mut_batch_coords, file = all_mut_coords_path, delim = "\t", col_names = FALSE)
}

read_in_bait_set_covered_regions=function(dir,stringency,file_suffix="_t2_1_Covered.bed") {
  files=list.files(dir,pattern = file_suffix,full.names = T)
  files<-grep(stringency,files,value=T)
  print("The following files have been found and will be read in:")
  print(files)
  batch_list=lapply(files,function(file) {return(read.delim(file,skip = 2,header=FALSE))})
  cat_bed=Reduce(rbind,batch_list)
  cat_bed<-cat_bed[,-4]
  colnames(cat_bed) <- c("Chr","Start","End")
  return(cat_bed)
}

read_in_covered_muts=function(included_muts_output_path) {
  covered_muts=read.delim(included_muts_output_path)
  covered_muts<-covered_muts[,1:2]
  covered_muts[,2] <- sapply(covered_muts[,2], function(x) gsub("^\\s+|\\s+$","",x)) #Strip the white space - the bedtools output has white space
  covered_coords=apply(covered_muts,1,paste,collapse="-")
  covered_coords<-gsub("chr","",covered_coords)
  return(covered_coords)
}


#---------------------------WRITE THE MUTATIONS FOR RUNNING THROUGH AGILENT---------------------------
# tree_folder="~/R_work/Zurich_HSCT/Data/tree_files"
# annotated_muts_folder="~/R_work/Zurich_HSCT/Data/annot_files"

tree_folder="filtering_runs2/tree_files"
annotated_muts_folder="filtering_runs2/annotated_muts"
output_dir="Baitset3/"
exp_nos=c(21,24,25)
#exp_nos=c(11,13,21,24,25,28,31,38,40,41)

#Set up the all.muts and all.trees objects - if previously run, read in the RDS object, otherwise import the details matrices from the annotated muts objects
if(file.exists(paste0(output_dir,"all.muts.RDS"))) {
  cat("Importing pre-existing mutation file",sep="\n")
  all.muts<-readRDS(paste0(output_dir,"all.muts.RDS"))
} else {
  annotated_muts_paths=list.files(annotated_muts_folder,pattern="vaf_post_mix_post_dup",full.names = T)
  all.muts<-lapply(exp_nos,function(exp_no){load(grep(paste0("Pair",exp_no),annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
  names(all.muts)<-paste0("Pair",exp_nos)
  all.muts<-all.muts[c("Pair21","Pair24","Pair25")]
}

tree_paths=list.files(tree_folder,pattern="vaf_post_mix_post_dup.tree",full.names = T)
all.trees<-lapply(exp_nos,function(exp_no){read.tree(grep(paste0("Pair",exp_no),tree_paths,value = T))})
names(all.trees)<-paste0("Pair",exp_nos)

all.trees<-all.trees[names(all.muts)]

#Write the mutation bed coordinates in batches of 80000 to run through the Agilent system at different stringencies (80,000 = maximum allowed for the agilent SureDesign system)
dir.create(paste0(output_dir,"All_mut_coords"))
temp=Map(details=all.muts, pair=names(all.muts), function(details,pair) {
  write_target_batches(details=details,sampleID=pair,batch_size = 80000,coord_buffer = 0,output_dir=output_dir)
  write_mut_coords_bedfile(details,all_mut_coords_path = paste0(output_dir,"All_mut_coords/",pair,"_mut_coords.bed")) #Then write a bedfile of all the actual mutation coordinates to intersect with the covered regions of the high/ mod stringency sets
})


write_target_batches(details=dplyr::bind_rows(all.muts),sampleID="HSCT_Set3_combined",batch_size = 80000,coord_buffer = 0,output_dir=output_dir)
write_mut_coords_bedfile(dplyr::bind_rows(all.muts),all_mut_coords_path = paste0(output_dir,"All_mut_coords/HSCT_Set3_combined_mut_coords.bed")) #Then write a bedfile of all the actual mutation coordinates to intersect with the covered regions of the high/ mod stringency sets


#---------------------------THEN HAVE TO RUN THESE BATCHES THROUGH THE AGILENT SUREDESIGN WEBSITE FOR EACH STRINGENCY (manual step)---------------

#---------------------------READ IN THE RESULTS FROM AGILENT AND DECIDE MUTATIONS TO INCLUDE-------------------------------------------------

if(!"bait_set"%in%colnames(all.muts[[1]])) {
  all.muts<-lapply(all.muts,function(details) mutate(details,bait_set="none"))
  
  #Run these functions for each stringency, starting from the lowest to the highest
  for(stringency in c("least","mod","most")) {
    cat_bed=read_in_bait_set_covered_regions(dir=paste0(output_dir,"/Covered_regions_beds"),stringency=stringency,file_suffix="_Covered.bed")
    cat_bed$Chr<-gsub("chr","",cat_bed$Chr)
    covered_regions_GR<-GenomicRanges::makeGRangesFromDataFrame(cat_bed,seqnames.field = "Chr",start.field = "Start",end.field = "End")
    
    all.muts<-lapply(all.muts,function(details) {
      details_GR<-GenomicRanges::makeGRangesFromDataFrame(details,seqnames.field = "Chrom",start.field="Pos",end.field = "Pos")
      overlaps<-GenomicRanges::findOverlaps(covered_regions_GR,details_GR)
      details$bait_set[overlaps@to] <- stringency #Record that each mutation at that coordinate is covered by the baitset of that stringency
      return(details)
    })
  }
}

lapply(all.muts,function(details) table(details$bait_set))

#Add the shearwater error rate information
shearwater_np_path="/lustre/scratch116/casm/cgp/users/nw14/resources/swaterNp_v1.bed.gz"
all.muts<-Map(details=all.muts,pair=names(all.muts),f=function(details,pair) {
  print(pair)
  if("sw_err_rate"%in%colnames(details)){
    return(details) #if the error rate is already included in the df, return it unchanged
  } else {
    mut_coords_path=paste0(output_dir,"DR",pair,"_loci.txt")
    shearwater_filtered_path=paste0(output_dir,"DR",pair,"_sw_output.txt")
    if(!file.exists(shearwater_filtered_path)) { #If the filtered shearwater file already exists, no need to RERUN
      temp=details[order(details$Chrom,details$Pos),c("Chrom","Pos")]
      colnames(temp)<-c("CHROM","POS")
      write.table(temp,file = mut_coords_path,sep="\t",col.names = F,quote = F,row.names = F)
      print("Running tabix to filter the shearwater reference file for the relevant locations")
      system(paste0("tabix -R ",mut_coords_path," ",shearwater_np_path,">",shearwater_filtered_path))
    }
    zz=read.table(shearwater_filtered_path,stringsAsFactors=FALSE) #Read in the shearwater filem and name the columns
    colnames(zz)<-c("chr","pos","ref","AF","TF","GF","CF","DELF","INSF","AB","TB","GB","CB","DELB","INSB","AR","TR","GR","CR","DELR","INSR")
    all_count_cols=c("AF","TF","GF","CF","DELF","INSF","AB","TB","GB","CB","DELB","INSB") #Define which columns refer to 'counts'
    
    #Now iterate through the df, and define the relevant base-specific error rate for each mutation
    details$sw_err_rate=NA
    for(i in 1:nrow(details)) {
      if(i%%10000==0) {print(i)}
      Chrom=details$Chrom[i]
      Pos=details$Pos[i]
      if(details$Mut_type[i]=="SNV") {
        Alt=details$Alt[i]
      } else if (details$Mut_type[i]=="INDEL") {
        if(nchar(details$Ref[i])>1) {
          Alt="DEL"
        } else {
          Alt="INS"
        }
      }
      mut_cols=all_count_cols[grepl(Alt,all_count_cols)]
      dep<-zz%>%filter(chr==Chrom & pos==Pos)%>%dplyr::select(all_of(all_count_cols))%>%sum()
      mtr<-zz%>%filter(chr==Chrom & pos==Pos)%>%dplyr::select(all_of(mut_cols))%>%sum()
      details$sw_err_rate[i]<-(mtr+1)/(dep+1)
    }
    return(details) #Return the updated details df that includes the base-specific error rate
  }
})

#Add the replication information
if(!file.exists(paste0(output_dir,"All_probegroups_HSCT_Set3_by_target"))) {
  replication_info_df=read.delim(paste0(output_dir,"HSCT_Set3_least_all_Probes_forMike.txt"))
  replication_info_df<-replication_info_df%>%filter(Replication!="Replication")%>%mutate(Replication=as.numeric(Replication))
  dim(replication_info_df)
  
  #Create a version of the df where each target is unique, and the replication is the sum of the replication number for each probe for that target
  unique_targets<-unique(replication_info_df$TargetID)
  length(unique_targets)
  library(parallel)
  replication_info_df=dplyr::bind_rows(mclapply(unique_targets,function(target) {
    if(which(unique_targets==target)%%100==0){print(which(unique_targets==target))}
    target_df=replication_info_df%>%filter(grepl(target,TargetID))
    return(data.frame(TargetID=target,ProbeIDs=paste0(target_df$ProbeID,collapse=","),mean_rep=mean(target_df$Replication)))
  },mc.cores=1))
  replication_info_df$Chrom=gsub("chr","",stringr::str_split(replication_info_df$TargetID,pattern=":",simplify=T)[,1])
  replication_info_df$Pos=as.numeric(stringr::str_split(stringr::str_split(replication_info_df$TargetID,pattern=":",simplify=T)[,2],pattern="-",simplify = T)[,1])
  write.table(replication_info_df,file=paste0(output_dir,"All_probegroups_HSCT_Set3_by_target"),sep = "\t",quote=F,row.names=F)
 } else {
 replication_info_df<-read.delim(paste0(output_dir,"All_probegroups_HSCT_Set3_by_target"),stringsAsFactors = F)
}

all.muts<-lapply(all.muts, function(details) {
  if("mean_rep"%in%colnames(details)) {
    return(details)
  } else {
    details<-left_join(details,replication_info_df%>%select(Chrom,Pos,ProbeIDs,mean_rep),by=c("Chrom","Pos"))
    return(details)
  }
})

# if(!"Replication"%in%colnames(details)|TRUE) {
#   details$Replication<-NULL
#   reps_path=paste0(output_dir,"DR",PairID,"_by_pos.tsv")
#   reps_df=read.table(reps_path,stringsAsFactors=FALSE,header=T)
#   reps_df$Pos=reps_df$Pos-1
#   reps_df$ChromPos=paste(reps_df$Chrom,reps_df$Pos,sep="-")
#   details=dplyr::left_join(details,reps_df[,c("ChromPos","Replication")],by="ChromPos")
# }

#SELECT THE MOST STRINGENT OF THE PRIVATE MUTS (rho > 0.3, vaf > 0.35), THEN RANDOMLY SELECT FROM EACH PRIVATE BRANCH TO MAKE A TOTAL

for(i in 1:length(all.muts)) {
  pair<-names(all.muts)[i]
  print(pair)
  details<-all.muts[[i]]
  if("stringent"%in%colnames(details)) {next}
  mats_and_params_file<-paste0("filtering_runs/mats_and_params/mats_and_params_",pair,"_m40_postMS_reduced")
  print("Loading up the mats and params file")
  load(mats_and_params_file)
  XY_low_depth_cutoff = 3; XY_high_depth_cutoff = 11; AUTO_low_depth_cutoff = 5; AUTO_high_depth_cutoff = 23
  print("Running filtering with more stringent parameters")
  stringent_muts=get_filtered_mut_set(input_set_ID = paste(pair,"_stringent"),  #the Run_ID of the unfiltered mutation set used as input - though won't account for any removal of samples from set
                                      COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                      filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                      gender = COMB_mats$gender, #patient's gender
                                      
                                      #These parameters decide whether a mutation is retained in the "true somatic mutation" set
                                      retain_muts = NULL,  #any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                      germline_pval = -10,  #the log10 p-value cutoff for mutations coming from an expected germline distribution
                                      rho = 0.3,  #rho cutoff for the beta-binomial filter, a measure of how "over-dispersed" the counts are compared to a binomial distribution
                                      mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order 1. lower threshold for autosomes, 2. upper threshold for autosomes, 3. lower threshold for XY, 4. upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                      pval_dp2=NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads
                                      pval_dp3=0.001,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping)
                                      min_depth = c(10,5), #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                      min_pval_for_true_somatic = 0.1,   #Default: 0.1. the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic
                                      min_vaf = c(0.35,0.8), #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY)
                                      
                                      #These parameters decide the genotype for each sample for each "true somatic mutation".  These may be less stringent than the initial parameters.
                                      min_variant_reads_SHARED = 2,  #the minimum number of reads for samples to be assigned a positive genotype
                                      min_pval_for_true_somatic_SHARED = 0.05,  #the p-value for coming from "true somatic mutation" read distribution to be assigned a positive genotype
                                      min_vaf_SHARED = NA) #Numeric vector of length 2, defining minimum vaf to be assigned a positive genotype
  
  #Select only the private muts that are in the more stringent mutation set
  details$stringent <- ifelse(details$mut_ref %in% stringent_muts$COMB_mats.tree.build$mat$mut_ref,1,0)
  all.muts[[i]]<-details
  print("Removing the mats and params files")
  rm(filter_params)
  rm(details)
  rm(stringent_muts)
  rm(COMB_mats)
}

lapply(all.muts,function(df) table(df$stringent))

saveRDS(all.muts,file=paste0(output_dir,"all.muts.RDS"))

#---------------------------NOW SELECT THE MUTATIONS TO INCLUDE BASED ON COLLECTED CRITERIA-------------------------------------------------

#SELECT THE SHARED & PRIVATE MUTS FROM DETAILS
#For shared mutations, select out the SNVs, and the ones that pass at least "mod stringency" masking
#Also include any coding change

samples_for_specificity_analysis=c()

all.muts<-Map(details=all.muts,pair=names(all.muts),tree=all.trees,f=function(details,pair,tree) {
  print(pair)
  details<-details[!duplicated(details$mut_ref),]
  shared_muts = details$mut_ref[!details$node %in% 1:length(tree$tip.label) &
                                  ((details$Mut_type == "SNV" &
                                      (details$bait_set == "most"|details$bait_set=="mod") &
                                      details$mean_rep <10)|
                                     (details$coding_change=="Coding change" & details$bait_set!="none"))]
  #Any shared mutation that is mod or most string
  
  #Loop to add in lower stringency mutations ("least string") if <0.5 of shared branches are included (e.g. short early branches)
  for(i in unique(details$node)) {
    if(i %in% 1:length(tree$tip.label)) {
      next
    } else if (length(details$mut_ref[details$node==i & details$mut_ref %in% shared_muts])/length(details$mut_ref[details$node==i]) > 0.5) {
      next
    } else {
      add_muts <- details$mut_ref[details$node==i & details$bait_set == "least" & details$sw_err_rate<2e-4]
      shared_muts <- c(shared_muts,add_muts)
    }
  }
  
  #For private mutations, select out SNVs, "most stringent" masking, and a more stringent filter parameter (as above)
  private_muts = details$mut_ref[details$node %in% 1:length(tree$tip.label) & details$Mut_type == "SNV" & details$bait_set == "most" &details$mean_rep<=2] #Only private mutations that pass stringent masking, stringent filtering parameters AND have low error rates on shearwater
  
  #Select all the "important mutations" (SNVs and indels) - but good to check if pass masking
  important_muts=details$mut_ref[details$coding_change_chip=="Coding change mutation in driver"]
  details[details$coding_change_chip=="Coding change mutation in driver",]
  
  #Need to decide how many private mutations to include for each pair
  total_private_muts_for_bait_set=10000
  total_private_muts=sum(details$node%in%1:length(tree$tip.label))
  to_include_proportion=total_private_muts_for_bait_set/total_private_muts
  
  private_muts_to_include_by_sample=lapply(1:length(tree$tip.label),function(i) {
    edge_muts<-details[get_edge_info(tree,details,node=i)$idx,]
    n_muts_to_include_for_branch=floor(nrow(edge_muts)*to_include_proportion)
    
    #Include minimum of 5 mutations for a branch
    n_muts_to_include_for_branch=max(5,n_muts_to_include_for_branch)
    
    reduced_muts<-edge_muts%>%
      filter(mut_ref%in%private_muts)%>%
      arrange(desc(stringent),sw_err_rate)%>%
      slice_head(n=n_muts_to_include_for_branch)%>%
      pull(mut_ref)
    
    return(reduced_muts)
  }
  )
  names(private_muts_to_include_by_sample) <- tree$tip.label
  unlist(lapply(private_muts_to_include_by_sample,length))
  print(sum(unlist(lapply(private_muts_to_include_by_sample,length)))) #This should be approximately the "total_private_muts_for_bait_set" figure
  private_muts_to_include=unlist(private_muts_to_include_by_sample)
  
  #Collect the more complete set of mutations for the 5 selected samples for specificity analysis
  selected_sample_nodes=which(tree$tip.label%in%samples_for_specificity_analysis)
  if(length(selected_sample_nodes)>0) {
    sensitivity_analysis_muts=details$mut_ref[details$node%in%selected_sample_nodes &
                                                details$bait_set!="none" &
                                                details$mean_rep<10]
  } else {
    sensitivity_analysis_muts=NULL
  }
  
  details$include_in_baitset <-ifelse(details$mut_ref %in% c(shared_muts,private_muts_to_include,important_muts,sensitivity_analysis_muts),1,0)
  return(details)
})
lapply(all.muts,function(details) table(details$include_in_baitset))

saveRDS(all.muts,file=paste0(output_dir,"all.muts.RDS"))

#---------------------------SUMMARISE BAITSET DETAILS & VISUALIZE---------------------------

pdf(paste0(output_dir,"Baitset_trees.pdf"),width = 15,height = 8)
temp=Map(details=all.muts,pair=names(all.muts),tree=all.trees,f=function(details,pair,tree) {
  print("Below is a summary of the masking stringency of included mutations:")
  print(table(details$bait_set[details$include_in_baitset==1]))
  print(paste("There are",sum(details$include_in_baitset[details$node %in% 1:length(tree$tip.label)]),"private mutations in the baitset"))
  print(paste("There are",sum(details$include_in_baitset[!details$node %in% 1:length(tree$tip.label)]),"shared mutations in the baitset"))
  print("Below is a summary of the coding mutations in driver genes with their masking stringencies:")
  print(details[details$coding_change_chip=="Coding change mutation in driver",c("mut_ref","Mut_type","variant_ID","sw_err_rate","bait_set","stringent","mean_rep","include_in_baitset")])
  
  #Visualize the covered mutations
  details$is.in.baitset <- details$include_in_baitset == 1
  tree=plot_tree(tree, cex.label = 0)
  add_annotation(tree=tree,
                 details=details,
                 list(),
                 annot_function=function(tree,details,matrices,node) {
                   add_binary_proportion(tree,details,matrices,node,bfield = "is.in.baitset",lwd = 1.5)
                 }
  )
})
dev.off()

#---------------------------WRITE THE FINAL MUTATION SETS TO INCLUDE---------------------------

#Bind the details matrices for the different individuals
details_comb=dplyr::bind_rows(all.muts)

included=details_comb%>%filter(include_in_baitset==1)%>%dplyr::select(Chrom,Pos,Ref,Alt)
write.table(included,file="All_included_muts",quote=F,sep="\t",row.names = F)

#Save sets to run at the different stringencies (combined across all individuals for this panel)
write_target_batches(details_comb[details_comb$bait_set=="most" & details_comb$include_in_baitset==1,],sampleID = "Comb_most_string",batch_size = 50000,coord_buffer=1,output_dir = output_dir)
write_target_batches(details_comb[details_comb$bait_set=="mod" & details_comb$include_in_baitset==1,],sampleID = "Comb_mod_string",batch_size = 50000,coord_buffer=1,output_dir = output_dir)
write_target_batches(details_comb[details_comb$bait_set=="least" & details_comb$include_in_baitset==1,],sampleID = "Comb_least_string",batch_size = 50000,coord_buffer=1,output_dir = output_dir)


#----FINAL CHECK FOR WHICH MUTATIONS ARE COVERED (after resubmitting final designs)------------------------------

cat_bed=read_in_bait_set_covered_regions(dir=paste0(output_dir),stringency="S3401023",file_suffix="_Covered.bed")
cat_bed$Chr<-gsub("chr","",cat_bed$Chr)
covered_regions_GR<-GenomicRanges::makeGRangesFromDataFrame(cat_bed,seqnames.field = "Chr",start.field = "Start",end.field = "End")

all.muts<-lapply(all.muts,function(details) {
  details$covered_final=0
  
  #Intersect each details coordinates with the final covered regions
  details_GR<-GenomicRanges::makeGRangesFromDataFrame(details,seqnames.field = "Chrom",start.field="Pos",end.field = "Pos")
  overlaps<-GenomicRanges::findOverlaps(covered_regions_GR,details_GR)
  details$covered_final[overlaps@to] <- 1 #Record that each mutation at that coordinate is covered by the final baitset
  return(details)
})

#---------------------------FINAL BAITSET SUMMARY AND VISUALIZATION---------------------------
library(ggplot2)
pdf(paste0(output_dir,"Baitset_trees_final.pdf"),width = 15,height = 8)
temp=Map(details=all.muts,tree=all.trees,pair=names(all.muts),f=function(details,tree,pair) {
  print("Below is a summary of the masking stringency of included mutations:")
  print(table(details$bait_set[details$covered_final==1]))
  print(paste("There are",sum(details$covered_final[details$node %in% 1:length(tree$tip.label)]),"private mutations in the baitset"))
  print(paste("There are",sum(details$covered_final[!details$node %in% 1:length(tree$tip.label)]),"shared mutations in the baitset"))
  print("Below is a summary of the coding mutations in driver genes with their masking stringencies:")
  print(details[details$coding_change_chip=="Coding change mutation in driver",c("Gene","covered_final","bait_set","include_in_baitset")])
  
  private_covered_prop=sapply(1:length(tree$tip.label), function(node) {
    total_muts=tree$edge.length[tree$edge[,2]==node]
    covered=sum(details$covered_final[details$node==node])
    return(covered/total_muts)
  })
  
  shared_covered_prop=sapply((length(tree$tip.label)+2):nrow(tree$edge), function(node) {
    total_muts=tree$edge.length[tree$edge[,2]==node]
    covered=sum(details$covered_final[details$node==node])
    return(covered/total_muts)
  })
  
  #Plot this info
  p1<-rbind(data.frame(proportion=private_covered_prop,type="private"),data.frame(proportion=shared_covered_prop,type="shared"))%>%
    ggplot(aes(x=proportion))+
    geom_histogram(fill="cadetblue3",col="black",size=0.5)+
    theme_classic()+
    facet_wrap(~type,ncol=2,scales="free")+
    labs(x="Proportion covered by bait set",
         y="Count",
         title=paste("Covered mutations by branch: Pair",pair))
  plot(p1)
  
  #Visualize the covered mutations
  details$is.in.baitset <- details$covered_final == 1
  tree=plot_tree(tree, cex.label = 0)
  add_annotation(tree=tree,
                 details=details,
                 list(),
                 annot_function=function(tree,details,matrices,node) {
                   add_binary_proportion(tree,details,matrices,node,bfield = "is.in.baitset",lwd = 1.5)
                 }
  )
})
dev.off()

##Write final bed files of mutations with co-ordinates covered by the baitset - across all covered individuals

all.muts%>%
  dplyr::bind_rows()%>%
  filter(covered_final==1 & Mut_type=="SNV")%>%
  dplyr::select(Chrom,Pos,Ref,Alt)%>%
  write.table(file="baitset3_SNVs.bed",quote=F,sep="\t",row.names=F,col.names = F)

all.muts%>%
  dplyr::bind_rows()%>%
  filter(covered_final==1 & Mut_type=="INDEL")%>%
  dplyr::select(Chrom,Pos,Ref,Alt)%>%
  write.table(file="baitset3_INDELs.bed",quote=F,sep="\t",row.names=F,col.names = F)

  
  