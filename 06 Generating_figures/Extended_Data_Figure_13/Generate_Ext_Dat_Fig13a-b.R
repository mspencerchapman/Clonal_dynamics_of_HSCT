#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","data.table")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}

#========================================#
# Set the ggplot2 themes for plotting ####
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
vcf_header_path=paste0(root_dir,"/data/reference_files/vcfHeader.txt") #A dummy vcf header that can be used

## Read in Pair metadata data frame----
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

## Define colour themes for the Pairs & DorR----
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

## Read in other data objects ----
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))

#Read in the versions of the tree/ details object used for the Bayesian model - these are subtly different
all_pairs<-Pair_metadata$Pair
tree_folder=paste0(root_dir,"/data/trees_no_dups")
annotated_muts_folder=paste0(root_dir,"/data/annot_files_no_dups")
tree_paths=list.files(tree_folder,pattern="vaf_post_mix_post_dup.tree",full.names = T)
all.trees.cc.nodups<-lapply(all_pairs,function(PairID){read.tree(grep(PairID,tree_paths,value = T))})
names(all.trees.cc.nodups)<-all_pairs

annotated_muts_paths=list.files(annotated_muts_folder,pattern="vaf_post_mix_post_dup",full.names = T)
all.muts.nodups<-lapply(all_pairs,function(PairID){cat(PairID);load(grep(PairID,annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
names(all.muts.nodups)<-all_pairs

## Generate information regarding loss-of-Y in male samples from X and Y coverage data ----
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Create dataframe of the LOY events
Pair11_LOY_nodes=c(48, 373,409 ,450 ,156 ,155 ,493 ,541,626)
Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                    Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                    Protein="",Type="",SO_codes="",coding_change="Coding change",
                                    coding_change_chip="yes",
                                    ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))


## Read in the spreadsheet listing other copy number changes ----
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#========================================#
# Import the targeted sequencing metadata ####
#========================================#

## First read in the metadata for all the targeted sequencing samples
bulk_smry_all<-readr::read_csv(file = paste0(root_dir,"/data/targeted_sequencing_data/targeted_sequencing_metadata.csv"))

#Read in other data objects
annotated_drivers_df<-read_csv(paste0(root_dir,"/data/Possible_drivers_annotated.csv"))[,c("mut_ref","node","Gene","variant_ID","Decision")]
remove_names=function(x) {attr(x,"names")<-NULL;return(x)}

#========================================#
# Import the raw targeted sequencing count data ####
#========================================#

#Import the raw targeted sequencing results - count data across all samples/ mutations sites
#This is generated using alleleCounter across all the targeted sequencing bams
all_targeted_res<-readRDS(paste0(root_dir,"/data/Targeted_sequencing_data/all_targeted_res.Rds"))

#========================================#
# GET CELL FRACTIONS FROM THE MODEL ####
#========================================#

## Import the output of the phylogeny-aware Gibb's sampler that infers the posterior distribution of the cell fraction of each mutation ----
## For scripts to run the model using the count data and the phylogenies, please see separate directory 'Gibbs Sampler for Targ Seq'

posterior_cell_fracs_file=paste0(root_dir,"/data/Targeted_sequencing_data/posterior_cell_fracs.Rds")

if(file.exists(posterior_cell_fracs_file)) {
  cat("Reading in previously saved RDS file of the posterior cell fractions",sep="\n")
  posterior_cell_fracs<-readRDS(posterior_cell_fracs_file)  ## THIS FILE NEEDS TO BE DOWNLOADED FROM MENDELEY DATA
} else {
  cat("Importing posterior VAFs files",sep="\n")
  setwd(paste0(root_dir,"/data/Targeted_sequencing_data/raw_model_output"))  ## ALTERNATIVELY CAN DOWNLOAD THE RAW MODEL OUTPUT WHICH IS THEN PROCESSED IN THESE FUNCTIONS
  posterior_VAFs<-lapply(all_pairs,function(pair) {
    cat(pair,sep="\n")
    tissueIDs<-bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(tissueID)%>%unique()
    pair_posteriors<-lapply(tissueIDs,function(ID) {
      cat(ID,sep="\n")
      posterior_VAFs_file=paste0(ID,"/",ID,"_posterior_VAFs.txt")
      tissue_posteriors<-readr::read_delim(posterior_VAFs_file,delim = "\t",show_col_types = FALSE)
      return(tissue_posteriors)
    })
    names(pair_posteriors)<-tissueIDs
    return(pair_posteriors)
  })
  names(posterior_VAFs)<-all_pairs
  
  #Convert the VAF posterior to a clonal fraction posterior
  #Multiply VAFs of diploid chromosome mutations two
  cat("Converting VAFs to cell fractions",sep="\n")
  pair_sex=sapply(Pair_metadata$Pair,function(pair) {Pair_metadata$donor_sex[Pair_metadata$Pair==pair]})
  autosomes=as.character(1:22)
  posterior_cell_fracs<-Map(pair_posteriors=posterior_VAFs,sex=pair_sex,pair=all_pairs,function(pair_posteriors,sex,pair) {
    cat(pair,sep="\n")
    tissueIDs=names(pair_posteriors)
    pair_cell_fracs_posteriors<-Map(tissue_post=pair_posteriors,tissueID=tissueIDs,function(tissue_post,tissueID) {
      cat(tissueID,sep="\n")
      #Select & separate out the posterior values into individual columns
      VAF_distributions=tissue_post%>%
        dplyr::select(Posterior_VAFs)%>%
        separate(col="Posterior_VAFs",into=paste("n",1:100,sep="_"),sep=",")%>%
        mutate_all(as.numeric)
      
      if(sex=="Male"&F){ #The VAFs of XY mutations in males have now been divided by 2 in the original algorithm
        #For males, identify the XY mutations and only multiply the autosomal mutaiton VAFs by 2
        cat("Treating sample as male",sep="\n")
        Chrom=tissue_post%>%
          dplyr::mutate(Chrom=stringr::str_split(mutation_ID,pattern="-",simplify=T)[,1])%>%
          dplyr::pull(Chrom)
        cell_frac_distributions=VAF_distributions
        cell_frac_distributions[Chrom%in%autosomes,]<-2*VAF_distributions[Chrom%in%autosomes,]
      } else if(sex=="Female"|T){
        cat("Treating sample as female",sep="\n")
        cell_frac_distributions=2*VAF_distributions
      }
      cell_frac_post<-bind_cols((tissue_post%>%dplyr::select(Node_assignment,mutation_ID)),cell_frac_distributions)
      return(cell_frac_post)
    })
    return(pair_cell_fracs_posteriors)
  })
  names(posterior_cell_fracs)<-all_pairs
  
  saveRDS(posterior_cell_fracs,file = posterior_cell_fracs_file)
}

#========================================#
# GET CLONAL CONTRIBUTIONS OF CLONAL EXPANSIONS TO THE DIFFERENT COMPARTMENTS ####
#========================================#
# How to define a 'clone' when looking at the phylogenetic tree? To avoid 'double-counting' you have to define a single time point
# when you consider the earliest time point for clonal expansions to begin, and any expansion beyond that point you consider as 'subclones'.
# If you define this point too early, you will start including normal expansion within development.
# We define it here as 100 mutations of molecular time. The developmental 'burst' of branch points is over by then, but it will capture the root of most abnormal expansions.

#Convenience function to work out the branches that traverse the defined cut off point.
get_cutoff_branches=function(tree,cut_off) {
  heights=nodeHeights(tree)
  cutoff_branches=tree$edge[,2][heights[,1]<cut_off & heights[,2]>=cut_off]
  return(cutoff_branches)
}

#Get the posterior distributions for each 'clone' in each cell type defined by cutting tree at 100 mutations of molecular time
clone_cutoff=100
clone_sizes=Map(tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=all_pairs,function(tree,post,pair) {
  cat(pair,sep="\n")
  cutoff_branches=get_cutoff_branches(tree,cut_off = clone_cutoff)
  tissueIDs=names(post)
  tissue_clone_fractions<-lapply(tissueIDs,function(tissueID) {
    cat(tissueID,sep="\n")
    tissue_post=post[[tissueID]]
    clone_posteriors<-lapply(cutoff_branches,function(node) {
      
      #For each branch, need to work out which mutation most closely corresponds corresponds to the '100 mutations of molecular time' cutoff
      
      #1. Find the molecular time of the ends of the branch that crosses the time point
      branch_heights=nodeHeights(tree)[tree$edge[,2]==node,]
      min_height=branch_heights[1]
      max_height=branch_heights[2]
      
      #2. From this, work out how far down the branch the '100 mutations of molecular time' is
      prop=(clone_cutoff-min_height)/(max_height-min_height)
      
      #3. Find out how many mutations from this branch are covered by the targeted sequencing.
      # Assuming that the order of these mutations along the branch can be deduced by ordering them by decreasing VAF (i.e. their rank), work out which rank most closely
      # corresponds to the branch position corresponding to '100 mutations o molecular time'
      n_targseq_muts=sum(tissue_post$Node_assignment==node)
      which_branch_rank=max(round(prop*n_targseq_muts),1)
      
      #4. Pull out the posterior distributions of mutations on this branch only, and convert into a matrix
      node_cell_frac_distributions=tissue_post%>%
        filter(Node_assignment==node)%>%
        dplyr::select(-Node_assignment,-mutation_ID)%>%
        as.matrix()
      
      #5. Sort the posterior cell fractions by their median values (to approximate the order) & select the posterior which corresponds to the rank defined above
      if(nrow(node_cell_frac_distributions)>1) {
        
        #Keeps values of specific mutations linked together & just orders them by their median values
        rank=order(apply(node_cell_frac_distributions,1,median),decreasing = T)
        node_cell_frac_distributions_sorted<-node_cell_frac_distributions[rank,]
        
      } else {
        node_cell_frac_distributions_sorted<-node_cell_frac_distributions
      }
      
      clone_post=node_cell_frac_distributions_sorted[which_branch_rank,]
      return(clone_post)
    })
    names(clone_posteriors)<-paste("node",cutoff_branches,sep="_")
    return(clone_posteriors)
  })
  names(tissue_clone_fractions)<-tissueIDs
  return(tissue_clone_fractions)
})
names(clone_sizes)<-all_pairs

#Plot clone sizes across cells types
clones_df=Map(pair=clone_sizes,name=names(clone_sizes),function(pair,name) {
  tissues=names(pair)
  pair_fracs<-Map(tissue=pair,ID=tissues,function(tissue,ID) {
    data.frame(ID=ID,node=names(tissue),cell_frac=sapply(tissue,median))
  })%>%dplyr::bind_rows()%>%mutate(Pair=name)
  return(pair_fracs)
})

#Get a list of all clones with a fraction over 1%in donor or recipient
include_nodes_0.01<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type,time_point)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  filter(time_point==0)%>%
  #filter(cell_type%in%c("Monocytes","Granulocytes"))%>%
  mutate(uid=paste(node,Pair,sep="_"))%>%
  filter(cell_frac>0.01)%>%
  arrange(cell_frac)%>%
  pull(uid)%>%
  unique(fromLast=T)
length(include_nodes_0.01) #114 total clones over 1% clonal fraction in at least one compartment

#========================================#
# COMPARE CLONAL FRACTIONS between LYMPHOID & MYELOID compartments ####
#========================================#

#Get plots where significance of difference is included by point size
relative_output_B_cells_df<-Map(clones=clones_df,clone_size=clone_sizes,post=posterior_cell_fracs,pair=names(all_targeted_res),function(clones,clone_size,post,pair) {
  cat(pair,sep="\n")
  
  #Include nodes that are >1% in any cell fraction in donor or recipient
  clone_differences<-lapply(c("Donor","Recipient"),function(DorR) {
    clones_enhanced<-clones%>%
      left_join(bulk_smry_all%>%dplyr::select(tissueID,individual_type,cell_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
      mutate(uid=paste(node,Pair,sep="_"))%>%
      filter(cell_type%in%c("Monocytes","B_cells") & individual_type==DorR & uid%in%include_nodes_0.01)%>%
      dplyr::select(-ID)%>%
      pivot_wider(names_from="cell_type",values_from="cell_frac")
    
    if(nrow(clones_enhanced)==0) {stop(return(NULL))}
    
    #Test if B and myeloid fractions in donor are significantly different
    clones_enhanced$signif=sapply(clones_enhanced$node,function(clone_number) {
      B_cell_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="B_cells" & time_point==0)%>%pull(tissueID)%>%.[1]
      Mono_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="Monocytes" & time_point==0)%>%pull(tissueID)%>%.[1]
      
      ratio_post=clone_size[[B_cell_donor_id]][[clone_number]]/clone_size[[Mono_donor_id]][[clone_number]]
      quants=quantile(ratio_post,c(0.025,0.5,0.975))
      signif=ifelse(quants['2.5%']>1|quants['97.5%']<1,T,F)
      return(signif)
    })
    return(clones_enhanced)
  })%>%dplyr::bind_rows()
  
  return(clone_differences)
})%>%dplyr::bind_rows()

## Generate Fig. 13a ----
B_myeloid_output_comparison_plot<-relative_output_B_cells_df%>%
  mutate(size=ifelse(signif,1.5,0.75))%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(x=B_cells,y=Monocytes,col=Pair,size=size))+
  geom_point(alpha=0.5)+
  scale_x_log10(limits=c(0.0009,1))+
  scale_y_log10(limits=c(0.0009,1))+
  scale_size_identity(breaks=c(0.75,1.5),labels=c("No diff. myeloid - lymphoid","Signif. diff. myeloid - lymphoid"),guide = "legend")+
  geom_abline(linetype=2)+
  facet_grid(~individual_type)+
  scale_color_manual(values=Pair_cols,drop=T)+
  guides(col=guide_legend(order=2),size=guide_legend(order=1))+
  theme_bw()+
  my_theme+
  labs(x="B cell clonal fraction",y="Monocyte clonal fraction",size="",col="")+
  theme(legend.key.size = unit(3,"mm"))
ggsave(paste0(plots_dir,"ExtDatFig13a.pdf"),B_myeloid_output_comparison_plot,width=6,height=2.1)

relative_output_T_cells_df<-Map(clones=clones_df,clone_size=clone_sizes,post=posterior_cell_fracs,pair=names(all_targeted_res),function(clones,clone_size,post,pair) {
  cat(pair,sep="\n")
  
  #Include nodes that are >1% in any cell fraction in donor or recipient
  clone_differences<-lapply(c("Donor","Recipient"),function(DorR) {
    clones_enhanced<-clones%>%
      left_join(bulk_smry_all%>%dplyr::select(tissueID,individual_type,cell_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
      mutate(uid=paste(node,Pair,sep="_"))%>%
      filter(cell_type%in%c("Monocytes","T_cells") & individual_type==DorR & uid%in%include_nodes_0.01)%>%
      dplyr::select(-ID)%>%
      pivot_wider(names_from="cell_type",values_from="cell_frac")
    
    if(nrow(clones_enhanced)==0) {stop(return(NULL))}
    
    #Test if B and myeloid fractions in donor are significantly different
    clones_enhanced$signif=sapply(clones_enhanced$node,function(clone_number) {
      T_cell_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="T_cells" & time_point==0)%>%pull(tissueID)%>%.[1]
      Mono_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="Monocytes" & time_point==0)%>%pull(tissueID)%>%.[1]
      
      ratio_post=clone_size[[T_cell_donor_id]][[clone_number]]/clone_size[[Mono_donor_id]][[clone_number]]
      quants=quantile(ratio_post,c(0.025,0.5,0.975))
      signif=ifelse(quants['2.5%']>1|quants['97.5%']<1,T,F)
      return(signif)
    })
    return(clones_enhanced)
  })%>%dplyr::bind_rows()
  
  return(clone_differences)
})%>%dplyr::bind_rows()

## Generate Fig. 13b ----
T_myeloid_output_comparison_plot<-relative_output_T_cells_df%>%
  mutate(size=ifelse(signif,1.5,0.75))%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(x=T_cells,y=Monocytes,col=Pair,size=size))+
  geom_point(alpha=0.5)+
  scale_x_log10(limits=c(0.00075,1))+
  scale_y_log10(limits=c(0.00075,1))+
  scale_size_identity(breaks=c(0.75,1.5),labels=c("No diff. myeloid - lymphoid","Signif. diff. myeloid - lymphoid"),guide = "legend")+
  geom_abline(linetype=2)+
  facet_grid(~individual_type)+
  scale_color_manual(values=Pair_cols,drop=T)+
  guides(col=guide_legend(order=2),size=guide_legend(order=1))+
  theme_bw()+
  my_theme+
  labs(x="T cell clonal fraction",y="Monocyte clonal fraction",size="",col="")+
  theme(legend.key.size = unit(3,"mm"))
ggsave(paste0(plots_dir,"ExtDatFig13b.pdf"),T_myeloid_output_comparison_plot,width=6,height=2.1)

#How many of the 114 clones have significant differences in either B/ T cell comparisons in either donor or recipient
bind_rows(relative_output_B_cells_df,relative_output_T_cells_df)%>%filter(signif)%>%pull(uid)%>%unique()%>%length()

#ANSWER: 106 of the 114 cones have significant differences in at least one comparison.
