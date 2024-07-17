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
LOY_files=list.files(path=paste0(root_dir,"/data/SV_and_CNA_data/LOY_files"),pattern="meanCoverage",full.names = T)
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
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

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
# DEFINE CUSTOM PLOTTING FUNCTIONS ####
#========================================#
##Contains my own functions for plotting the targeted seqeuncing data on the tree - built on top of the 'plot_tree' package from Nick Williams
add_var_col=function(tree, ##<< enhanced phylo returned from plot_tree
                     details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                     node,
                     var_field,
                     b.add.line=TRUE,
                     colours = c("black","green","red"),
                     scale_muts_to_branch=TRUE,
                     ...){
  
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  muts_on_edge=length(info$idx.in.details)
  edge_length=tree$edge.length[tree$edge[,2]==node]
  
  if(muts_on_edge > 0 & edge_length>0) {
    bdat=details[[var_field]][info$idx]
    if(is.null(bdat) || class(bdat)!="numeric"){
      stop("Error in provided bfield (does it exist and is it numeric?)")
    }
    
    bdat = sort(bdat, decreasing = TRUE)
    if(scale_muts_to_branch) {
      mut_unit_of_edge=edge_length/muts_on_edge
    } else {
      mut_unit_of_edge=1
    }
    ##Could add in a third category NA
    #missing=sum(is.na(bdat))
    if(b.add.line){
      y0_next = info$yt
      for(i in 1:muts_on_edge) {
        arrows(y0=y0_next,y1=(y0_next - mut_unit_of_edge),x0=info$x,length = 0,col=col.scale[ceiling(100*bdat[i])],lend=1,...)
        y0_next = y0_next - mut_unit_of_edge
      }
    }
  }
}

####
Gibbs_targ_seq_plots=function(SampleID,
                              tree,
                              details_targ,
                              pair_cell_fracs,
                              scale_muts_to_branch=TRUE,
                              colour.scale=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                              log_min=-5,
                              vaf_lwd=5,
                              ...) {
  
  posterior_cell_fracs.SampleID <-pair_cell_fracs[[SampleID]]
  
  get_median_cellfracs=function(post.df) {
    require(dplyr)
    median_cell_fracs=post.df%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()%>%
      apply(1,median)
    return(data.frame(mutation_ID=post.df$mutation_ID,median_cell_frac=median_cell_fracs))
  }
  
  details_targ_full<-left_join(details_targ,
                               get_median_cellfracs(post.df=posterior_cell_fracs.SampleID),
                               by=c("mut_ref"="mutation_ID"))
  
  #Generate the rescaled log cell fraction for plotting with contrast
  details_targ_full$log_median_cell_frac<-log(details_targ_full$median_cell_frac)
  details_targ_full$log_median_cell_frac[details_targ_full$log_median_cell_frac<log_min]<-log_min
  details_targ_full$log_median_cell_frac=plotrix::rescale(details_targ_full$log_median_cell_frac,newrange = c(0,1))
  
  ##Generate the plot
  tree=plot_tree(tree,cex.label=F,...)
  add_annotation(tree=tree,
                 details=details_targ_full,
                 matrices=NULL,
                 annot_function=function(tree,details,matrices,node) {
                   add_var_col(tree,
                               details,
                               node,
                               var_field = "log_median_cell_frac",
                               lwd = vaf_lwd,
                               colours=colour.scale,
                               scale_muts_to_branch=scale_muts_to_branch)
                 }
  )
}

#========================================#
# Now do the actual plots ####
#========================================#

## Generate Extended Data Fig. 11a-b ----
pdf(file=paste0(plots_dir,"ExtDatFig11.pdf"),width=7,height=10)
par(mfrow=c(2,1))
pair_new_id="Pair_7"
DorR="Donor"
selected_cell_type="Granulocytes"
tissue_id<-bulk_smry_all%>%filter(Pair_new==pair_new_id &individual_type==DorR & cell_type==selected_cell_type & time_point==0)%>%dplyr::select(Pair,tissueID)%>%filter(!duplicated(.))

temp=Gibbs_targ_seq_plots(SampleID=tissue_id$tissueID,
                     tree=all.trees.cc.nodups[[tissue_id$Pair]],
                     details_targ=all_targeted_res[[tissue_id$Pair]]$details_targ,
                     pair_cell_fracs=posterior_cell_fracs[[tissue_id$Pair]],
                     scale_muts_to_branch=TRUE,
                     colour.scale=RColorBrewer::brewer.pal(9,"RdPu")[c(2:9,9,9)],
                     log_min=-6,
                     vaf_lwd=2,
                     title=paste(pair_new_id,"combined tree: Deep targeted sequencing of",tolower(selected_cell_type),"from",tolower(DorR)),
                     plot_axis=F,
                     lwd=0.2,
                     cex.title=0.6)

DorR="Recipient"
tissue_id<-bulk_smry_all%>%filter(Pair_new==pair_new_id &individual_type==DorR & cell_type==selected_cell_type & time_point==0)%>%dplyr::select(Pair,tissueID)%>%filter(!duplicated(.))

temp=Gibbs_targ_seq_plots(SampleID=tissue_id$tissueID,
                     tree=all.trees.cc.nodups[[tissue_id$Pair]],
                     details_targ=all_targeted_res[[tissue_id$Pair]]$details_targ,
                     pair_cell_fracs=posterior_cell_fracs[[tissue_id$Pair]],
                     scale_muts_to_branch=TRUE,
                     colour.scale=RColorBrewer::brewer.pal(9,"RdPu")[c(2:9,9,9)],
                     log_min=-6,
                     vaf_lwd=2,
                     title=paste(pair_new_id,"combined tree: Deep targeted sequencing of",tolower(selected_cell_type),"from",tolower(DorR)),
                     plot_axis=F,
                     lwd=0.2,
                     cex.title=0.6)
dev.off()
