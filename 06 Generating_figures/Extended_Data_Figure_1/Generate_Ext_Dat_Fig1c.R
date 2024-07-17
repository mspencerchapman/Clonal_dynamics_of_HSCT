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
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
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
tree_folder=paste0(root_dir,"/data/tree_and_mutation_files/trees_no_dups")
annotated_muts_folder=paste0(root_dir,"/data/tree_and_mutation_files/annot_files_no_dups")
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

## Read in the spreadsheet listing other copy number changes ----
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

#Create dataframe of the LOY events
Pair11_LOY_nodes=c(48, 373,409 ,450 ,156 ,155 ,493 ,541,626)
Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                    Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                    Protein="",Type="",SO_codes="",coding_change="Coding change",
                                    coding_change_chip="yes",
                                    ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))

## Read in mutational signature extraction data----
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)



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
# CAPTURED CELL FRACTIONS THROUGH TIME ####
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

#Function to get the 'sum of cell fractions' across the posterior distribution for a tissue at any given molecular time
get_sum_of_frac=function(tissue_post,tree,cut_off,verbose=F) {
  
  if(verbose) {cat(paste("Clone cutoff of",cut_off,"being used"),sep="\n")}
  
  #Define the cutoff branches
  cutoff_branches=get_cutoff_branches(tree,cut_off = cut_off)
  
  if(verbose) {cat(paste(length(cutoff_branches)," branches at the cutoff level"),sep="\n")}
  
  #If the cutoff branches have no mutations in the targeted sequencing data, need to replace with the daughter branches of that node
  number_of_covered_muts=function(node,tissue_post) {sum(tissue_post$Node_assignment==node)} #convenience function
  
  while(any(sapply(cutoff_branches,number_of_covered_muts,tissue_post=tissue_post)==0)){
    updated_nodes<-unlist(lapply(cutoff_branches,function(node) {
      if(number_of_covered_muts(node,tissue_post=tissue_post)==0) {
        if(verbose) {cat(paste("No mutations on node",node,"=> replacing with daughter nodes."),sep="\n")}
        daughter_nodes<-tree$edge[,2][tree$edge[,1]==node]
        return(daughter_nodes)
      } else {
        return(node)
      }
    }))
    cutoff_branches<-updated_nodes
  }
  
  #Define the clone posteriors (as a list)
  clone_posteriors<-lapply(cutoff_branches,function(node) {
    if(verbose) {cat(node,sep="\n")}
    branch_heights=nodeHeights(tree)[tree$edge[,2]==node,]
    min_height=branch_heights[1]
    max_height=branch_heights[2]
    
    if(min_height>cut_off) {
      prop<-0
    } else {
      prop=(cut_off-min_height)/(max_height-min_height)
    }
    
    n_targseq_muts=sum(tissue_post$Node_assignment==node)
    
    which_branch_rank=max(round(prop*n_targseq_muts),1)
    node_cell_frac_distributions=tissue_post%>%
      filter(Node_assignment==node)%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()
    
    #Sort each iteration of posterior by cell fraction & choose that which corresponds to the rank
    if(nrow(node_cell_frac_distributions)>1) {
      #Original method - sorts cell fracs of each Gibbs sample, and takes the 10th rank of each. ?significant overestimates
      #node_cell_frac_distributions_sorted<-apply(node_cell_frac_distributions,2,sort,decreasing=T)
      
      #Second method - keeps values of specific mutations linked together & just orders them by their median values
      rank=order(apply(node_cell_frac_distributions,1,median),decreasing = T)
      node_cell_frac_distributions_sorted<-node_cell_frac_distributions[rank,]
      
    } else {
      node_cell_frac_distributions_sorted<-node_cell_frac_distributions
    }
    
    clone_post=node_cell_frac_distributions_sorted[which_branch_rank,]
    return(clone_post)
  })
  names(clone_posteriors)<-paste("node",cutoff_branches,sep="_")
  
  #Use these to get the sum of VAF
  n_iter=length(clone_posteriors[[1]])
  sum_of_frac_post<-sapply(1:100,function(j) {
    sum(sapply(clone_posteriors,function(x) x[j]))
  })
  return(sum_of_frac_post)
}

#Because this is a fairly lengthy function to run, there is a presaved version of the output in the 'Data' folder
#However, this can be regenerated if desired

sum_of_frac_list_file=paste0(root_dir,"/data/Targeted_sequencing_data/sum_of_fracs.Rds")
if(file.exists(sum_of_frac_list_file)){
  sum_of_frac_list<-readRDS(sum_of_frac_list_file)
} else {
  #This is slow
  sum_of_frac_list<-Map(tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=names(all_targeted_res),function(tree,post,pair) {
    cat(pair,sep="\n")
    tree.ultra<-make.ultrametric.tree(tree)
    tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(drop.tip(tree,"Ancestral")))
    tree.ultra$edge.length[tree.ultra$edge[,2]==which(tree.ultra$tip.label=="Ancestral")]<-0
    
    tissueIDs=names(post)
    pair_sum_of_frac_by_cutoff<-lapply(tissueIDs,function(tissueID){
      cat(tissueID,sep="\n")
      tissue_post<-post[[tissueID]]
      tissue_sum_of_frac_by_cutoff<-lapply(seq(2,100,2),function(clone_cutoff) {
        cat(clone_cutoff,sep="\n")
        get_sum_of_frac(tissue_post = tissue_post,cut_off = clone_cutoff,tree = tree.ultra,verbose=F)
      })
      return(tissue_sum_of_frac_by_cutoff)
    })
    names(pair_sum_of_frac_by_cutoff)<-tissueIDs
    return(pair_sum_of_frac_by_cutoff)
  })
  saveRDS(sum_of_frac_list,file=sum_of_frac_list_file)
}

## Extract the median and 95% CI of the sum of frac statistic ----
all_sof<-Map(pair=names(all_targeted_res),pair_sum_of_frac_by_cutoff=sum_of_frac_list,f=function(pair,pair_sum_of_frac_by_cutoff) {
  pair_cell_fracs<-Map(tissue=pair_sum_of_frac_by_cutoff,tissueID=names(pair_sum_of_frac_by_cutoff),f=function(tissue,tissueID) {
    cat(tissueID,sep="\n")
    tissue_df<-Map(cutoff=seq(2,100,2),cutoff_sof=tissue,function(cutoff,cutoff_sof) {
      data.frame(cutoff=cutoff,lowerCI=quantile(cutoff_sof,0.025),median=median(cutoff_sof),upperCI=quantile(cutoff_sof,0.975))
    })%>%dplyr::bind_rows()%>%
      mutate(tissueID=tissueID,.before=1)
    return(tissue_df)
  })%>%dplyr::bind_rows()%>%
    mutate(PairID=pair,.before=1)
  return(pair_cell_fracs)
})


## Generate plot Extended Data Fig 1c ----
SOF_plot<-all_sof%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,individual_type,cell_type)%>%filter(!duplicated(.)),by="tissueID")%>%
  dplyr::mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")),
                individual_type=factor(individual_type))%>%
  dplyr::mutate(Pair_new=factor(Pair_new,levels=paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(ymin=lowerCI,ymax=upperCI,y=median,x=cutoff,col=individual_type,fill=individual_type))+
  geom_hline(yintercept = 1,linetype=2)+
  geom_ribbon(alpha=0.4,col=NA)+
  geom_line(size=0.3)+
  facet_grid(Pair_new~cell_type)+
  scale_y_continuous(breaks=seq(0,1,0.5),limits=c(0,1))+
  scale_fill_manual(values=remove_names(DorR_cols))+
  scale_color_manual(values=remove_names(DorR_cols))+
  theme_bw()+
  my_theme+
  labs(x="Molecular Time",y="Captured cell fraction",col="",fill="")

ggsave(filename=paste0(plots_dir,"ExtDatFig1c.pdf"),SOF_plot,width=4.8,height=5)
