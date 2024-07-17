#1. Import the SNV trees that are corrected for sequencing coverage
#2. Reduce branch lengths by the proportion caused by non-ubiquitous mutational process i.e. APOBEC/ mutagenic exposures
#3. Reduce terminal branch lengths by the number of 'in vitro' or 'non linear' (e.g. mutations caused from increased replication through cell division during differentiation) additional mutations
#4. Then make ultrametric, scaling branches to the average mutation burden across all branches
#5. Estimate timing of coalescences during HSC life based on (1) 60 mutations at birth, (2) linear accumulation thereafter
#6. Create 'transplant time' shaded box based on time + poisson variation of mean mutation burden at that time

#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest")
bioconductor_packages=c()

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
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

#========================================#
# Set the ggplot2 theme for plotting ####
#========================================#

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 7),
                axis.title = element_text(size=8),
                axis.line = element_line(linewidth = 0.4),
                axis.ticks = element_line(linewidth = 0.3),
                legend.text = element_text(size=6),
                legend.title = element_text(size=8),
                strip.text = element_text(size=7),
                strip.background = element_rect(fill="lightgray",linewidth = 0.4),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#========================================#
# Set the root directory and read in the necessary files ####
#========================================#

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
source(paste0(root_dir,"/data/HSCT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")
tree_folder=paste0(root_dir,"/data/tree_files/")
annotated_muts_folder=paste0(root_dir,"/data/annot_files/")

#Read in Pair metadata data frame ----
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

#Define colour themes for the Pairs & DorR ----
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in the sample level metadata - this is compiled from:
#(1) information from Markus Manz's group about colony phenotype
#(2) CanApps project info &
#(3) DNA QC results
sample_metadata<-read.delim(paste0(root_dir,"/data/sample_level_metadata.tsv"))

#Generate information regarding loss-of-Y in male samples from X and Y coverage data ----
LOY_files=list.files(path=paste0(root_dir,"/data/SV_and_CNA_data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv("~/R_work/Zurich_HSCT/Data/Copy_number_changes.csv")

#Read in mutational signature extraction data
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)


#========================================#
# IMPORTING AND PREPARING THE DATA ####
#========================================#

##1. Import the SNV trees corrected for coverage & the details matrices ----
tree_paths=list.files(tree_folder,pattern=".tree",full.names = T)
all.trees<-lapply(exp_nos,function(exp_no){read.tree(grep(paste0("Pair",exp_no),tree_paths,value = T))})
names(all.trees)<-paste0("Pair",exp_nos)

annotated_muts_paths=list.files(annotated_muts_folder,pattern="annotated_mut_set",full.names = T)
all.muts<-lapply(exp_nos,function(exp_no){load(grep(paste0("Pair",exp_no),annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
names(all.muts)<-paste0("Pair",exp_nos)

#Create dataframe of the LOY events
Pair11_LOY_nodes=c(447, 468, 152, 407, 369, 56, 493, 539)
Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                    Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                    Protein="",Type="",SO_codes="",coding_change="Coding change",
                                    coding_change_chip="yes",
                                    ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))


##2. Re-annotate drivers using whichever custom list you want to use ----
myeloid_panel=read.csv(paste0(root_dir,"/data/Driver_lists/Illumina_trusight_myeloid_panel.csv"),header = F)$V1
bolton_panel=read.csv(paste0(root_dir,"/data/Driver_lists/Myeloid_drivers_bolton.csv"),header=F)$V1
all_ukbb_CH=read.csv(paste0(root_dir,"/data/Driver_lists/UKBB_drivers.csv"),header=F)$V1
driver_list<-c(myeloid_panel,bolton_panel,all_ukbb_CH)

all.muts<-lapply(all.muts,function(details) {
  details$coding_change <- ifelse(details$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                      "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                      "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                      "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                      "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                      "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss")|
                                    grepl("splice_site_variant",details$Type),
                                  "Coding change",
                                  "no")
  
  details$coding_change_chip<-ifelse(details$coding_change=="Coding change" & details$Gene%in%driver_list,"yes","no")
  details$ChromPos=paste(details$Chrom,details$Pos,sep="-")
  details$variant_ID=paste(details$Gene, details$Protein, sep = " ")
  
  return(details)
})

##Get annotated list of possible drivers with manual decisions for if these are "true" drivers ("Oncogenic" or "Possible")
#Read in previous annotations
annotated_driver_mut_set=read_csv(paste0(root_dir,"/data/Possible_drivers_annotated.csv"))[,c("mut_ref","node","Gene","variant_ID","Decision")]
#Get updated list of mutations in driver genes and join the previous annotations - write the updated file
Map(details=all.muts,pair=names(all.muts),function(details,pair) details%>%dplyr::filter(coding_change_chip=="yes")%>%dplyr::select(mut_ref,node,Gene,variant_ID)%>%mutate(Pair=pair))%>%
  dplyr::bind_rows()%>%
  left_join(annotated_driver_mut_set)%>%
  write_csv(file=paste0(root_dir,"/data/Possible_drivers.csv"))

#Manually update the "Possible_drivers.csv" spreadsheet with decisions & save as "Possible_drivers_annotated.csv" when complete. Now read in.
annotated_driver_mut_set=read_csv(paste0(root_dir,"/data/Possible_drivers_annotated.csv"))[,c("mut_ref","node","Gene","variant_ID","Decision")]

#3. Update the details matrices with the manual decisions on whether mutations are oncogenic/ possible oncogenic or VUS ----
all.muts<-lapply(all.muts,function(details) {
  details<-left_join(details,annotated_driver_mut_set%>%dplyr::select(mut_ref,Decision),by="mut_ref")
  return(details)
})


#4. Review duplicates and drop the duplicate samples ----
plot_duplicate_plate_map(tree=all.trees[[6]],
                         sample_metadata = sample_metadata%>%
                           mutate("position"=ifelse(nchar(Well)==3,Well,paste0(substr(Well,1,1),0,substr(Well,2,2))))%>%
                           mutate(plate_ID=paste(Zur_ID,Plate_no,sep="_")),
                         private_mut_threshold = 30,
                         labels="set")

#Spot the duplicates for dropping
drop.samples=Map(tree=all.trees,details=all.muts,function(tree,details) {
  private_mut_threshold=30
  print(paste("Threshold for determining duplicate samples is",private_mut_threshold,"mutations."))
  
  #Determine the sets of duplicate samples
  duplicate_samples=get_duplicate_sets(tree,mut_threshold = private_mut_threshold)
  lapply(duplicate_samples,function(x) print(paste(paste(x,collapse=" & "),"are recognised as duplicate samples")))
  
  #Choose one of each set of duplicate samples to keep - the one with the lower mutation burden (as suspect higher is usually in vitro acquired)
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
  return(drop_samples)
})

#Not all duplicates are recognised automatically, manually add additional samples
remove_samples_manual=list(Pair11=character(),
                           Pair13=character(),
                           Pair21=character(),
                           Pair24=character(),
                           Pair25=character(),
                           Pair28=c("PD45805b_lo0102","PD45805b_lo0146","PD45805b_lo0145","PD45805b_lo0118","PD45805b_lo0112"),
                           Pair31=character(),
                           Pair38=character(),
                           Pair40=character(),
                           Pair41=character())

remove_samples_list=Map(auto=drop.samples,manual=remove_samples_manual,function(auto,manual) {c(auto,manual)})

#Create a set of 'details' matrices without the duplicates and with nodes re-numbered
all.muts.nodups<-Map(tree=all.trees,details=all.muts,remove_samples=remove_samples_list,function(tree,details,remove_samples) {
  if(length(remove_samples)>0){
    output<-remove_samples_from_tree_and_update_details(remove_samples=remove_samples,tree=tree,details=details)
    return(output$details)
  } else {
    return(details)
  }
})

##5. Add sample-level data to the metadata table: SNV burden, peak VAF ----
sample_metadata=dplyr::bind_rows(lapply(names(all.trees),function(Pair_ID) {
  cat(paste0(Pair_ID,"\n"))
  tree<-all.trees[[Pair_ID]]
  annotated_muts_path=grep(Pair_ID,annotated_muts_paths,value=T)
  load(annotated_muts_path)
  
  
  sample_metadata_Ind<-sample_metadata%>%mutate(ID=paste0("Pair",Pair))%>%dplyr::filter(ID==Pair_ID)
  
  #Calculate the 'peak VAF' measure for the sample
  sample_metadata_Ind$peak_vaf=sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree$tip.label){
      return(NA)
    } else {
      max.density=function(x){
        dens<-density(x)
        return(dens$x[which.max(dens$y)])
      }
      mut_vafs=get_mut_vafs(SampleID,COMB_mats=filtered_muts$COMB_mats.tree.build,tree=tree)
      return(max.density(mut_vafs$vaf))
    }
  })
  
  #Calculate the mutation burden for the sample
  sample_metadata_Ind$SNV_burden_u=sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree$tip.label){
      return(NA)
    } else {
      sample_node=which(tree$tip.label==SampleID)
      sample_nodes=c(sample_node,getAncestors(tree = tree,node = sample_node,type="all"))
      return(sum(filtered_muts$COMB_mats.tree.build$mat$node%in%sample_nodes))
    }
  })
  
  #Calculate the uncorrected mutation burden for the sample
  sample_metadata_Ind$SNV_burden=sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree$tip.label){
      return(NA)
    } else {
      return(nodeheight(tree,node=which(tree$tip.label==SampleID)))
    }
  })
  
  return(sample_metadata_Ind)
}))

write.table(sample_metadata,file=paste0(root_dir,"/data/metadata_temp.tsv"),quote = F,sep = "\t",row.names = F)

#### HERE NEED TO RUN THE 'CORRECT FOR CLONALITY' SCRIPT TO ADD SOMATIC SNV SENSITIVITY INFORMATION
# See the script '01 Generating_the_phylogenies/programs_and_functions/Correct_for_clonality.R'
# However, this script also needs the 'mats_and_params' files, so you will need to download these from Mendeley Data in order to run

##6.Now calculate sensitivity based on both clonality & coverage ----
sample_metadata<-read.delim(paste0(root_dir,"/data/metadata_temp_updated.tsv"))
sens_df<-sample_metadata%>%
  dplyr::select(Sample,somatic_SNV_sensitivity)%>%
  dplyr::rename("SNV_sensitivity"=somatic_SNV_sensitivity)%>%
  rbind(data.frame(Sample="Ancestral",SNV_sensitivity=1))%>%
  filter(!duplicated(Sample))

get_corrected_tree=function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  #print(tree$tip.label[!tree$tip.label%in%sensitivity_df$Sample])
  if(!all(tree$tip.label%in%sensitivity_df$Sample)){stop(print("Not all samples in the tree are in the sensitivity dataframe"))}
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

##7. Get uncorrected set of trees ----
all.trees.uncorrected<-Map(tree=all.trees,details=all.muts,function(tree,details){
  tree$edge.length<-sapply(tree$edge[,2],function(node) return(sum(details$node==node & details$Mut_type=="SNV")))
  return(tree)
})
par(mfrow=c(5,2))
temp=lapply(all.trees.uncorrected,plot.phylo,show.tip.label=F,direction="downwards")

##8. Now get trees corrected for sensitivity & clonality ----
all.trees.cc.nodups<-Map(tree=all.trees.uncorrected,details=all.muts,remove_samples=remove_samples_list,function(tree,details,remove_samples) {
  tree_SNV<-get_subset_tree(tree,details,v.field = "Mut_type",value = "SNV")
  tree_SNV_c<-get_corrected_tree(tree_SNV,details,sensitivity_df = sens_df,include_SNVs = T,include_indels = F)
  tree_SNV_c<-drop.tip(tree_SNV_c,remove_samples)
})
par(mfrow=c(5,2))
temp=lapply(all.trees.cc.nodups,plot.phylo,show.tip.label=F,direction="downwards")

##9. Correct the branch lengths by removing mutations assigned to APOBEC/ platinum signature ----
all.trees.adj1=Map(exp_no=exp_nos,remove_samples=remove_samples_list,function(exp_no,remove_samples) {
  
  #Start from the raw tree with duplicates
  tree<-all.trees.uncorrected[[paste0("Pair",exp_no)]]
  
  #Remove APOBEC/ platinum signatures using the HDP output
  remove_sigs=paste0("N",c(3,4)) #Check that these are the right
  exposures_df_pair=exposures_df%>%dplyr::filter(Pair==paste0("Pair",exp_no))
  tree_adj<-tree
  tree_adj$edge.length=sapply(tree$edge[,2],function(node) {
    branch_length_original<-tree$edge.length[tree$edge[,2]==node]
    if(node %in% exposures_df_pair$node) {
      remove_prop<-sum(exposures_df_pair[exposures_df_pair$node==node,remove_sigs])
      return(branch_length_original*(1-remove_prop))
    } else {
      return(branch_length_original)
    }
  })
  tree_adj_c<-get_corrected_tree(tree_adj,details,sensitivity_df = sens_df,include_SNVs = T,include_indels = F,get_edge_from_tree=T)
  tree_adj_c.nodups<-drop.tip(tree_adj_c,remove_samples)
  return(tree_adj_c.nodups)
})
names(all.trees.adj1)<-names(all.trees)

if(check_plots) {
  par(mfrow=c(5,2))
  temp=lapply(all.trees.adj1,plot.phylo,show.tip.label=F,direction="downwards")
}

##10. Now correct terminal branches for in vitro/ differentiation-induced mutations by removing (in utero mutations - intercept) from linear regression  ----
in_vitro_mutation_burden=60
all.trees.adj2=lapply(all.trees.adj1,function(tree) {
  tree_adj<-tree
  tree_adj$edge.length<-sapply(1:nrow(tree$edge),function(i) {
    if(tree$edge[i,2]%in%1:length(tree$tip.label)) {
      if(tree$edge.length[i]<in_vitro_mutation_burden){
        print(paste("Branch length reduced to 1 for sample",tree$tip.label[tree$edge[i,2]]))
      }
      return(max(c(tree$edge.length[i] - in_vitro_mutation_burden,1)))
    } else {
      return(tree$edge.length[i])
    }
  })
  return(tree_adj)
})
names(all.trees.adj2)<-names(all.trees)

if(check_plots) {
  par(mfrow=c(5,2))
  temp=lapply(all.trees.adj2,plot.phylo,show.tip.label=F,direction="downwards")
}

##11. Add these adjusted mutation burdens to the metadata dataframe ----
sample_metadata=dplyr::bind_rows(lapply(names(all.trees),function(Pair_ID) {
  cat(paste0(Pair_ID,"\n"))
  tree<-all.trees[[Pair_ID]]
  details<-all.muts[[Pair_ID]]
  sample_metadata_Ind<-sample_metadata%>%dplyr::filter(ID==Pair_ID)
  
  tree_SNV<-get_subset_tree(tree,details,v.field = "Mut_type",value = "SNV")
  
  tree_adj1<-all.trees.adj1[[Pair_ID]]
  tree_adj2<-all.trees.adj2[[Pair_ID]]
  details_no_dups<-all.muts.nodups[[Pair_ID]]
  
  #plot(tree_SNV,show.tip.label=F,direction="downwards")
  var(get_mut_burden(tree_SNV))
  tree_SNV_c<-get_corrected_tree(tree_SNV,details,sensitivity_df = sens_df%>%dplyr::filter(Sample%in%tree$tip.label),include_SNVs = T,include_indels = F)
  var(get_mut_burden(tree_SNV_c))
  #plot(tree_SNV_c,show.tip.label=F,direction="downwards",cex=0.5)
  
  #Calculate the corrected mutation burden for the sample
  sample_metadata_Ind$SNV_burden=sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree_SNV_c$tip.label){
      return(NA)
    } else {
      return(nodeheight(tree_SNV_c,node=which(tree$tip.label==SampleID)))
    }
  })
  
  #Calculate the adjusted mutation burden - removing platinum & APOBEC signatures
  sample_metadata_Ind$SNV_burden_adj1=sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree_adj1$tip.label){
      return(NA)
    } else {
      return(nodeheight(tree_adj1,node=which(tree_adj1$tip.label==SampleID)))
    }
  })
  
  #Calculate the adjusted mutation burden - removing in vitro mutations
  sample_metadata_Ind$SNV_burden_adj2=sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree_adj2$tip.label){
      return(NA)
    } else {
      return(nodeheight(tree_adj2,node=which(tree_adj2$tip.label==SampleID)))
    }
  })
  
  return(sample_metadata_Ind)
}))

#12. Make the trees ultrametric ----
all.trees.ultra<-lapply(all.trees.adj2,function(tree) {
  mean_mutation_burden=mean(get_mut_burden(tree))
  tree.ultra<-make.ultrametric.tree(tree)
  tree.ultra$edge.length=tree.ultra$edge.length*mean_mutation_burden
  tree.ultra$edge.length[which(tree$edge[,2]==which(tree$tip.label=="Ancestral"))]<-0 #Set ancestral branch length back to 0
  return(tree.ultra)
})
names(all.trees.ultra)<-paste0("Pair",exp_nos)

#13.Edit the details matrices to have column indicating driver mutations on shared branches only ----
all.muts.nodups<-Map(details=all.muts.nodups,tree=all.trees.ultra,function(details,tree) {
  mutate(details,shared_coding_change_chip=ifelse(coding_change_chip!="no"&!node%in%1:length(tree$tip.label),"yes","no"))
})

#14. Get split Donor and Recipient tree and details objects ----
all.muts.nodups.D=Map(details=all.muts.nodups,tree=all.trees.ultra,function(details,tree) {
  DR_ids=get_DR_ids(tree=tree); remove_samples=grep(DR_ids[2],tree$tip.label,value = T)
  output<-remove_samples_from_tree_and_update_details(remove_samples=remove_samples,tree=tree,details=details)
  return(output$details)
})

all.muts.nodups.R=Map(details=all.muts.nodups,tree=all.trees.ultra,function(details,tree) {
  DR_ids=get_DR_ids(tree=tree); remove_samples=grep(DR_ids[1],tree$tip.label,value = T)
  output<-remove_samples_from_tree_and_update_details(remove_samples=remove_samples,tree=tree,details=details)
  return(output$details)
})

all.trees.ultra.D=lapply(all.trees.ultra,function(tree){
  DR_ids=get_DR_ids(tree=tree); remove_samples=grep(DR_ids[2],tree$tip.label,value = T)
  tree.D<-drop.tip(tree,remove_samples)
  tree.D$coords<-NULL
  return(tree.D)
})

all.trees.ultra.R=lapply(all.trees.ultra,function(tree){
  DR_ids=get_DR_ids(tree=tree); remove_samples=grep(DR_ids[1],tree$tip.label,value = T)
  tree.R<-drop.tip(tree,remove_samples)
  tree.R$coords<-NULL
  return(tree.R)
})

#Add mutational signature information to the sample metadata table before saving
#Get the signatures from individual samples - NB. this still only includes mutations on branches with >50 mutations
sigs<-dplyr::bind_rows(Map(tree=all.trees,Individual=names(all.trees),function(tree,Individual) {
  sigs=get_signatures_in_samples(tree = tree,signature_names=paste0("N",0:4),exposures_df = exposures_df%>%filter(Pair==Individual))
  return(sigs)
}))

#Join the signatures data to the sample metadata df
sample_metadata<-full_join(sigs,sample_metadata) #NB. this removes info from non-clonal samples/ other samples not included in the "trees with duplicates"
for(sig in paste0("N",0:4)) {
  sample_metadata[[paste0(sig,"_abs")]]<-sample_metadata[[sig]]*sample_metadata[["SNV_burden"]]
}
sample_metadata<-sample_metadata%>%
  filter(Sample!="Ancestral")
sample_metadata$Pair_new<-sapply(sample_metadata$ID,function(pair) return(Pair_metadata$Pair_new[Pair_metadata$Pair==pair]))

#Summarise the colony outcomes
contamination_colonies=c("PD45793b_lo0094","PD45793b_lo0095","PD45793b_lo0096",
                               "PD45803b_lo0081",
                               "PD45811b_lo0187","PD45811b_lo0188","PD45811b_lo0189","PD45811b_lo0190","PD45811b_lo0191","PD45811b_lo0192")
sample_metadata=sample_metadata%>%
  mutate(sample_status=ifelse(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)),"PASS",
                              ifelse(Sample%in%unlist(remove_samples_list),"Duplicate",
                                     ifelse(Coverage<4,"Low coverage",
                                            ifelse(is.na(Coverage),"Not Sequenced",
                                                   ifelse(Sample%in%contamination_colonies,"Different individual","Non-clonal"))))))

colony_outcomes<-c("PASS","Low coverage","Non-clonal","Different individual","Duplicate")
colony_outcome_cols=c("#33A02C","#A6CEE3","#1F78B4","#FB9A99","grey50")
names(colony_outcome_cols)<-colony_outcomes

colony_outcome_summary_plot<-sample_metadata%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  mutate(sample_status=factor(sample_status,levels=rev(colony_outcomes)))%>%
  filter(sample_status!="Not Sequenced")%>%
  ggplot(aes(x=DorR,fill=sample_status))+
  geom_bar(col="black",size=0.1)+
  theme_bw()+
  my_theme+
  theme(strip.text.y=element_text(angle=0,size=6))+
  facet_grid(rows=vars(Pair_new))+
  scale_fill_manual(values=colony_outcome_cols)+
  coord_flip()+
  labs(x="Donor or Recipient",y="Number of colonies",fill="Sample\noutcome")

ggsave(filename = paste0(plots_dir,"colony_summary_plot.pdf"),colony_outcome_summary_plot,width =4,height=3)

#Combine tree and details objects for saving
tree_objects<-c("all.trees",
                "all.trees.uncorrected",
                "all.trees.adj1",
                "all.trees.adj2",
                "all.trees.cc.nodups",
                "all.trees.ultra",
                "all.trees.ultra.D",
                "all.trees.ultra.R")
tree_lists<-lapply(tree_objects,get)
names(tree_lists)<-tree_objects

details_objects<-c("all.muts",
                   "all.muts.nodups",
                   "all.muts.nodups.D",
                   "all.muts.nodups.R")
details_lists<-lapply(details_objects,get)
names(details_lists)<-details_objects

#Save files
saveRDS(tree_lists,file=paste0(root_dir,"/data/tree_and_mutation_files/tree_lists.Rds"))
saveRDS(details_lists,file=paste0(root_dir,"/data/tree_and_mutation_files/details_lists.Rds"))
saveRDS(sample_metadata,file=paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
