
####################################################################################################
# 
# After completion of ABC computation at the end of "round 1", we have
# a sample of parameter vectors from the (approximate) posterior distribution;
# 
# At the start of "round 2", we take this pre-existing sample of parameter vectors from this (approximate) posterior distribution
# (OR a sub-sample drawn withOUT replacment from this pre-existing sample of parameter vectors),
# and for each parameter vector in this sample, we perform (a large number of) simulations
# from conditional distribution at each parameter vector;
# (by running script: resim_conditionals_transplant.R);
#
# This sample (the output from the script: resim_conditionals_transplant.R)
# is the input for the present script (compute_ppc_results_seq_transplant.R);
#
# The present script (compute_ppc_results_seq_transplant.R) uses this large sample of simulated data sets (from each conditional distribution)
# to compute moments of each conditional distribution;
#
# The present script (compute_ppc_results_seq_transplant.R) also collects and pools all these samples from the conditional distributions
# to obtain a stratified sample from the posterior predictive distribution;
#
####################################################################################################

## Run on farm5
args = commandArgs(TRUE)
args

source_directory = toString(args[1])
resim_directory = toString(args[2])
output_directory = toString(args[3])
posterior_sample_file_name = toString(args[4])
obs_stats_directory = toString(args[5])
obs_stats_file_name = toString(args[6])

min_resims_per_post_obs = as.integer(args[7])

model_name = toString(args[8])
stat_set_name = toString(args[9])
ABC_method = toString(args[10])
pair_ID = toString(args[11])

cat( "pair_ID = ", pair_ID, "\n", sep = "\t" )
# readLines("stdin",n=1)

####################################################################################################

# seq_resims_per_post_obs = seq( from = 50, to = 1000, by = 50 )

####################################################################################################
#!/software/R-4.1.0/bin/Rscript
library(optparse)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(phytools)
library(rsimpop)
library(abc)
library(ggridges)

####################################################################################################

# my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/"
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch126/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
setwd(tree_mut_dir); source("treemut.R") # ;setwd(my_working_directory)
visualize=F

####################################################################################################

###-----------------SPECIFIC FUNCTIONS FOR THIS SCRIPT-----------------
get_ltt = function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  ltt_tree = sapply(time_points, function(x) {
    sum(nodeheights[,1] < x & !nodeheights[,2] < x)
  })
  return(ltt_tree)
}

get_cutoff_branches=function(tree,cut_off) {
  heights=nodeHeights(tree)
  cutoff_branches=tree$edge[,2][heights[,1]<cut_off & heights[,2]>=cut_off]
  return(cutoff_branches)
}

get_coalescences = function(ltt) {
  coals=sapply(2:length(ltt), function(i) {return(ltt[i]-ltt[i-1])})
  return(coals)
}

coals_within_time_window_per_clone=function(tree,define_clone_height=5,time_points) {
  clone_nodes=get_cutoff_branches(tree,cut_off = define_clone_height)
  nodes_within_time_points=tree$edge[,1][nodeHeights(tree)[,1]>time_points[1] & nodeHeights(tree)[,1]<time_points[2]]
  
  coals_within_time_window_by_clone=sapply(clone_nodes,function(node) {
    clone_daughters<-get_all_node_children(node,tree)
    coals_within_time_points=intersect(clone_daughters,nodes_within_time_points)
    return(length(coals_within_time_points))
  })
  return(coals_within_time_window_by_clone)
}

get_DR_ids=function(tree){
  tree=ape::drop.tip(tree,"Ancestral")
  PD_IDs=unique(substr(tree$tip.label,1,8))
  #Get number elements only
  PD_numbers=readr::parse_number(PD_IDs)
  names(PD_IDs)<-sapply(PD_numbers,function(n) ifelse(n%%2==0,"donor_ID","recip_ID"))
  return(PD_IDs)
}

get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}

#Bulk totals for other expansions
find_latest_acquisition_node=function(tree,pos_samples){
  #Get list of ancestral nodes for all samples
  ancestral_nodes_list=lapply(pos_samples,function(Sample) {
    get_ancestral_nodes(node = which(tree$tip.label==Sample),edge=tree$edge)
  })
  #Find nodes that are ancestral to all the samples
  common_nodes=Reduce(intersect,ancestral_nodes_list)
  #Which of these is the most recent (i.e. has the maximum node height)
  nodeheights<-sapply(common_nodes,function(node) nodeheight(tree = tree,node = node))
  MRCA_node<-common_nodes[which.max(nodeheights)]
  return(MRCA_node)
}

#Function to extract clonal fractions from the total population based on a subsampled tree
#Need to make sure that the tip labels match
extract_bulk_cell_fractions=function(nodes,sub_pop,full_pop,states) {
  #Get a vector of the 'states' of the all tree tips (i.e. which compartment)
  full_pop_tip_states=full_pop$state[full_pop$edge[,2]<=length(full_pop$tip.label)]
  
  #Now go through each node from the sub_pop tree & find the full population fraction of the MRCA of that clade
  res<-lapply(nodes,function(node) {
    cat(node,sep = "\n")
    raw_tips<-stringr::str_split(getTips(sub_pop,node),pattern = "_",simplify = T)[,1] #Get the sample names from the clade
    bulk_node=find_latest_acquisition_node(tree=full_pop,pos_samples = raw_tips) #Get the MRCA of those samples from the full population
    clade_states=full_pop_tip_states[which(full_pop$tip.label%in%getTips(full_pop,node=bulk_node))] #Get the states of enclosed samples
    bulk_fracs=as.data.frame(table(clade_states)[as.character(states)]/table(full_pop_tip_states)[as.character(states)])%>%
      dplyr::mutate(clade_states=states)%>%
      tidyr::replace_na(replace = list(Freq=0))%>%
      pivot_wider(names_from = "clade_states",values_from="Freq",names_prefix="state_")
  })%>%dplyr::bind_rows()%>%
    mutate(node=nodes,.before=1)
  
  return(res)
}

#Do updated 'get_subsampled_tree' function that maintains original tip labels (needed to extract bulk clonal fractions)
get_subsampled_tree2=function (tree, N, tips = tree$edge[c(which(tree$state == 0 & 
                                                                   tree$edge[, 2] <= length(tree$tip.label)), sample(which(tree$state != 
                                                                                                                             0 & tree$edge[, 2] <= length(tree$tip.label)), N)), 2]) 
{
  N = length(tips)
  tip.labels<-tree$tip.label[sort(tips)]
  
  #Create the tree with the tips kept in the same order
  tree_same_order<-keep.tip(tree,tip=tip.labels)
  tmp = rsimpop:::C_subsample_pop(tree, sort(tips))
  tmp$tip.label = sprintf("s%d", 1:N)
  class(tmp) = c("simpop", "phylo")
  
  #Now reorder the tips
  node_trans<-all.equal(tree_same_order,tmp,use.tip.label=F,index.return = T)
  node_trans<-node_trans[node_trans[,2]<=length(tree_same_order$tip.label),]
  new_tips<-tip.labels[node_trans[,2]]
  tmp$tip.label<-new_tips
  
  #tmp$tip.label = tip.labels
  
  checkValidPhylo(tmp)
  tmp$is_combined = tree$is_combined
  tmp
}

####################################################################################################

###-----------------Pair metadata & trees-----------------

exp_nos<-c(11,13,21,24,25,28,31,38,40,41)
Pair_metadata=data.frame(Pair=paste0("Pair",exp_nos),
                         Pair_new=factor(c("Pair_9","Pair_7","Pair_5","Pair_1","Pair_4","Pair_10","Pair_6","Pair_8","Pair_3","Pair_2"),levels=paste0("Pair_",1:10)),
                         Age=c(74.8,65.5,64.5,34.2,58.4,79.9,65.2,65.8, 51.9,42.4),
                         Age_at_transplant=c(66,36,50,18,47,63,43,35,35,30),
                         MNC_dose=c(2.66,4.17,16.28,NA,10.94,13.9,NA,4.05,2.48,15.01),
                         CD34_dose=c(1.56,NA,7.1,7.9,8.97,2.4,NA,NA,NA,4.51),
                         stem_cell_source=c("BM","BM","PBSC","PBSC","PBSC","PBSC","BM","BM","BM","PBSC"),
                         conditioning=c("MAC","MAC","RIC","MAC","RIC","RIC","MAC","MAC","MAC","MAC"))

#Import the trees (the ones with duplicates removed)
# tree_folder="/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/filtering_runs2/tree_files"
tree_folder="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs2/tree_files"
tree_paths=list.files(tree_folder,pattern="vaf_post_mix_post_dup.tree",full.names = T)
all.trees<-lapply(exp_nos,function(exp_no){
  pair_tree_paths=grep(paste0("Pair",exp_no),tree_paths,value = T)
  if(length(pair_tree_paths)==1) {
    tree<-read.tree(pair_tree_paths)
  } else if(length(pair_tree_paths)>1){
    a_j_path<-grep("a_j",pair_tree_paths,value=T)
    tree<-read.tree(a_j_path)
  } else {
    stop(cat(paste("No tree for Pair",exp_no),sep="\n"))
  }
  return(tree)
})
names(all.trees)<-paste0("Pair",exp_nos)

ABC.trees<-readRDS("/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/trees_for_ABC.Rds")
plots_dir="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/plots/"

####################################################################################################

####################################################################################################
#
# Get:
# no_of_WGS_recip;
# no_of_WGS_donor;
###-----------------Set up parameters for this simulation-----------------

#Set the HSCT parameters
# PairID=paste0("Pair",j)
# PairID=paste0("Pair",pair_index)
PairID=pair_ID

cat( "pair_ID = ", pair_ID, "\n", sep = "\t" )

age_of_donor_at_HSCT=Pair_metadata$Age_at_transplant[Pair_metadata$Pair==PairID]
age_of_donor_at_sampling=Pair_metadata$Age[Pair_metadata$Pair==PairID]

tree.data<-all.trees[[PairID]]
DR_ids<-get_DR_ids(tree.data)

donor_id<-DR_ids['donor_ID'] # kjd
recip_id<-DR_ids['recip_ID'] # kjd

cat( "donor_id = ", donor_id, "\n", sep = "\t" ) # kjd
cat( "recip_id = ", recip_id, "\n", sep = "\t" ) # kjd


no_of_WGS_recip=sum(grepl(DR_ids['recip_ID'],tree.data$tip.label))
no_of_WGS_donor=sum(grepl(DR_ids['donor_ID'],tree.data$tip.label))

####################################################################################################

####################################################################################################

###--------------Get summary stats from the data-----------------

##EXTRACT SUMMARY STATISTICS
#(1) 3 LARGEST CLADES
# all.ss<-Map(tree=ABC.trees,pair=names(ABC.trees),function(tree,pair) {

# tree=ABC.trees[[ pair_index ]]
# pair=names(ABC.trees)[ pair_index ]

tree=ABC.trees[[ PairID ]]
pair=PairID

cat( "PairID = ", PairID, "\n", sep = "\t" ) # kjd
cat( "pair = ", pair, "\n", sep = "\t" ) # kjd
  
  donor_id<-get_DR_ids(tree)['donor_ID'] # kjd
  recip_id<-get_DR_ids(tree)['recip_ID'] # kjd
  
cat( "donor_id = ", donor_id, "\n", sep = "\t" ) # kjd
cat( "recip_id = ", recip_id, "\n", sep = "\t" ) # kjd
  
  age_of_donor_at_HSCT<-Pair_metadata$Age_at_transplant[Pair_metadata$Pair==pair]
  age_of_donor_at_sampling<-Pair_metadata$Age[Pair_metadata$Pair==pair]
  
  D_tree<-keep.tip(tree,tip=grep(donor_id,tree$tip.label))
  R_tree<-keep.tip(tree,tip=grep(recip_id,tree$tip.label))
  
  largest_clades_D<-get_expanded_clade_nodes(D_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]
  largest_clades_R<-get_expanded_clade_nodes(R_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]
  
  #(2) Number of singletons
  n_singletons_D<-get_expanded_clade_nodes(D_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    dplyr::filter(n_samples==1)%>%nrow(.)
  n_singletons_R<-get_expanded_clade_nodes(R_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    dplyr::filter(n_samples==1)%>%nrow(.)
  
  #(3) Peri-HSCT LTT/ coalescences
  
  tree.time<-tree
  tree.time$edge.length<-tree$edge.length*(age_of_donor_at_sampling/median(get_mut_burden(tree)))
  D_tree.time<-keep.tip(tree.time,tip=grep(donor_id,tree$tip.label))
  R_tree.time<-keep.tip(tree.time,tip=grep(recip_id,tree$tip.label))
  
  # peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
  if( stat_set_name == "original" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
  
  }else if( stat_set_name == "peri_interval_narrower" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-2.5,x+2.5)))
  
  }else if( stat_set_name == "peri_interval_wider" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-7.5,x+7.5)))
  
  }else if( stat_set_name == "pre_interval_divided" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5))) # same as "original"!
  
  }else if( stat_set_name == "peri_interval_divided" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x,x+5))) # notice extra point (at "x" = age_of_donor_at_HSCT) in vector c(x-5,x,x+5);
  
  }
  
  DR_ltt_peri<-get_ltt(tree=tree.time,time_points=peri_HSCT_time_points)
  D_ltt_peri<-get_ltt(tree=D_tree.time,time_points=peri_HSCT_time_points)
  R_ltt_peri<-get_ltt(tree=R_tree.time,time_points=peri_HSCT_time_points)
  DR_coals_peri<-get_coalescences(DR_ltt_peri)
  D_coals_peri<-get_coalescences(D_ltt_peri)
  R_coals_peri<-get_coalescences(R_ltt_peri)
  
  #(4) Pre-HSCT LTT/ coalescences
  
  # pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
  if( stat_set_name == "original" )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
  
  }else if( stat_set_name == "peri_interval_narrower" )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1]) # same as "original"; except that peri_HSCT_time_points[1] has been moved up!
  
  }else if( stat_set_name == "peri_interval_wider"  )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1]) # same as "original"; except that peri_HSCT_time_points[1] has been moved down!
  
  }else if( stat_set_name == "pre_interval_divided" )
  {
	mid_point = 5 + ((peri_HSCT_time_points[1] - 5)/2)
	pre_HSCT_time_points=c(5,mid_point,peri_HSCT_time_points[1]) # notice extra point ("mid_point");
  
  }else if( stat_set_name == "peri_interval_divided" )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1]) # same as "original"!
  
  }
  
  DR_ltt_pre<-get_ltt(tree=tree.time,time_points=pre_HSCT_time_points)
  D_ltt_pre<-get_ltt(tree=D_tree.time,time_points=pre_HSCT_time_points)
  R_ltt_pre<-get_ltt(tree=R_tree.time,time_points=pre_HSCT_time_points)
  DR_coals_pre<-get_coalescences(DR_ltt_pre)
  D_coals_pre<-get_coalescences(D_ltt_pre)
  R_coals_pre<-get_coalescences(R_ltt_pre)
  
  #(5) Maximum peri-transplant coalescences within single recipient clade (this is a new stat)
  R_max_coals_in_single_clade=max(coals_within_time_window_per_clone(tree=R_tree.time,define_clone_height = 5,time_points = peri_HSCT_time_points))
  
  visualize=F
  if(visualize){
    par(mfrow=c(2,1))
    zz=plot_tree(D_tree.time,cex.label=0)#,default_edge_color = "#11a0aa80")
    rect(xleft = 0,xright = length(D_tree.time$tip.label),ybottom = (age_of_donor_at_sampling-age_of_donor_at_HSCT-5),ytop = (age_of_donor_at_sampling-age_of_donor_at_HSCT+5),col="#D3D3D380",border = "black")
    yy=plot_tree(R_tree.time,cex.label=0)#,default_edge_color = "#c8256580")
    rect(xleft = 0,xright = length(R_tree.time$tip.label),ybottom = (age_of_donor_at_sampling-age_of_donor_at_HSCT-5),ytop = (age_of_donor_at_sampling-age_of_donor_at_HSCT+5),col="#D3D3D380",border = "black")
  }
  
  #(6) Bulk stats
  mean_abs_log2FC=NA
  median_abs_log2FC=NA
  max_abs_log2FC=NA
  mean_abs_change=NA
  median_abs_change=NA
  max_abs_change=NA
  
  ##--COMBINE SUMSTATS INTO SINGLE VECTOR--
#  ss_names=c(paste("largest_clade_D",1:3,sep="_"),
#             paste("largest_clade_R",1:3,sep="_"),
#             "n_singletons_D",
#             "n_singletons_R",
#             paste(c("DR","D","R"),"coals_pre",sep="_"),
#             paste(c("DR","D","R"),"coals_peri",sep="_"),
#             "R_max_coals_in_single_clade",
#             "mean_abs_log2FC",
#             "median_abs_log2FC",
#             "max_abs_log2FC",
#             "mean_abs_change",
#             "median_abs_change",
#             "max_abs_change")
  
  
  if( stat_set_name %in% c( "original", "peri_interval_narrower", "peri_interval_wider" ) )
  {
	ss_names=c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             "DR_coals_pre",
             "D_coals_pre",
             "R_coals_pre",
             "DR_coals_peri",
             "D_coals_peri",
             "R_coals_peri",
             "R_max_coals_in_single_clade",
             "mean_abs_log2FC",
             "median_abs_log2FC",
             "max_abs_log2FC",
             "mean_abs_change",
             "median_abs_change",
             "max_abs_change")
	
  }else if( stat_set_name == "pre_interval_divided" )
  {
	ss_names=c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             "DR_coals_pre_1", # new
             "DR_coals_pre_2", # new
             "D_coals_pre_1", # new
             "D_coals_pre_2", # new
             "R_coals_pre_1", # new
             "R_coals_pre_2", # new
             "DR_coals_peri",
             "D_coals_peri",
             "R_coals_peri",
             "R_max_coals_in_single_clade",
             "mean_abs_log2FC",
             "median_abs_log2FC",
             "max_abs_log2FC",
             "mean_abs_change",
             "median_abs_change",
             "max_abs_change")
	
  }else if( stat_set_name == "peri_interval_divided" )
  {
	ss_names=c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             "DR_coals_pre",
             "D_coals_pre",
             "R_coals_pre",
             "DR_coals_peri_1", # new
             "DR_coals_peri_2", # new
             "D_coals_peri_1", # new
             "D_coals_peri_2", # new
             "R_coals_peri_1", # new
             "R_coals_peri_2", # new
             "R_max_coals_in_single_clade",
             "mean_abs_log2FC",
             "median_abs_log2FC",
             "max_abs_log2FC",
             "mean_abs_change",
             "median_abs_change",
             "max_abs_change")
	
  }
  
  if( stat_set_name %in% c( "original", "peri_interval_narrower", "peri_interval_wider" ) )
  {
	ss_comb=c(largest_clades_D,
            largest_clades_R,
            n_singletons_D,
            n_singletons_R,
            DR_coals_pre,
            D_coals_pre,
            R_coals_pre,
            DR_coals_peri,
            D_coals_peri,
            R_coals_peri,
            R_max_coals_in_single_clade,
            mean_abs_log2FC,
            median_abs_log2FC,
            max_abs_log2FC,
            mean_abs_change,
            median_abs_change,
            max_abs_change)
	  names(ss_comb)<-ss_names
	  
  }else if( stat_set_name == "pre_interval_divided" )
  {
	
	DR_coals_pre_1 = DR_coals_pre[ 1 ] # new
    DR_coals_pre_2 = DR_coals_pre[ 2 ] # new
    D_coals_pre_1 = D_coals_pre[ 1 ] # new
    D_coals_pre_2 = D_coals_pre[ 2 ] # new
    R_coals_pre_1 = R_coals_pre[ 1 ] # new
    R_coals_pre_2 = R_coals_pre[ 2 ] # new
    
	ss_comb=c(largest_clades_D,
            largest_clades_R,
            n_singletons_D,
            n_singletons_R,
            DR_coals_pre_1, # new
            DR_coals_pre_2, # new
            D_coals_pre_1, # new
            D_coals_pre_2, # new
            R_coals_pre_1, # new
            R_coals_pre_2, # new
            DR_coals_peri,
            D_coals_peri,
            R_coals_peri,
            R_max_coals_in_single_clade,
            mean_abs_log2FC,
            median_abs_log2FC,
            max_abs_log2FC,
            mean_abs_change,
            median_abs_change,
            max_abs_change)
	  names(ss_comb)<-ss_names
	  
  }else if( stat_set_name == "peri_interval_divided" )
  {
	
	DR_coals_peri_1 = DR_coals_peri[ 1 ] # new
    DR_coals_peri_2 = DR_coals_peri[ 2 ] # new
    D_coals_peri_1 = D_coals_peri[ 1 ] # new
    D_coals_peri_2 = D_coals_peri[ 2 ] # new
    R_coals_peri_1 = R_coals_peri[ 1 ] # new
    R_coals_peri_2 = R_coals_peri[ 2 ] # new
	
	ss_comb=c(largest_clades_D,
            largest_clades_R,
            n_singletons_D,
            n_singletons_R,
            DR_coals_pre,
            D_coals_pre,
            R_coals_pre,
            DR_coals_peri_1, # new
            DR_coals_peri_2, # new
            D_coals_peri_1, # new
            D_coals_peri_2, # new
            R_coals_peri_1, # new
            R_coals_peri_2, # new
            R_max_coals_in_single_clade,
            mean_abs_log2FC,
            median_abs_log2FC,
            max_abs_log2FC,
            mean_abs_change,
            median_abs_change,
            max_abs_change)
	  names(ss_comb)<-ss_names
	  
  }
  
#  return(ss_comb)
# })

####################################################################################################

obs_stats_table = as.data.frame( matrix( NA, nrow = 1, ncol = length( ss_names ) ) )
names( obs_stats_table ) <- ss_names

# obs_stats_table[ 1, cols_stats ] <- unlist(stats_data)

obs_stats_table[ 1, ss_names ] = ss_comb

####################################################################################################

if( model_name == "m1" )
{
  
  # param_names_vect=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "fitnessGammaFn", "dpcpd" )
  param_names_vect=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
  
}else if( model_name == "m2" )
{
  
  # param_names_vect=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
  param_names_vect=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd", "gamma_shape_engraftment", "gamma_rate_engraftment" )
  
}else if( model_name == "m3" )
{
  
  # param_names_vect=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
  param_names_vect=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd", "exaggeration_ratio", "exaggeration_prop" )
	
}

####################################################################################################
####################################################################################################
#
# Read in simulations;
#

resim_file_list = list.files( path = resim_directory, pattern = "resim_" )

n_resim_files = length( resim_file_list )

cat( "0.0. main: length( resim_file_list ) = ", length( resim_file_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)

setwd( resim_directory )
# resim_table_list = Map( function(x){ read.table( file = x, sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE ) }, sim_trees_file_list )

resim_table_list = list()
for( i in 1:length( resim_file_list ) )
{
	resim_file = resim_file_list[[ i ]]
	
	if( file_test( "-f", resim_file ) )
	{
		resim_result = try( resim_table <- read.table( file = resim_file, sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE ) )
		if( class( resim_result ) != "try-error" )
		{
			resim_table_list = append( resim_table_list, list( resim_table ) )
			
		}
		
	}
	
}
n_resim_tables = length( resim_table_list )
cat( "0.0. main: length( resim_table_list ) = ", length( resim_table_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)

if( length( resim_table_list ) > 0 ) 
{
	if( length( resim_table_list ) == 1 ) 
	{
		# pooled_resim_table = resim_table_list[[ 1 ]]
		stratified_resim_stats_table = resim_table_list[[ 1 ]]
		
	}else
	{
		# pooled_resim_table <- do.call( "rbind", resim_table_list )
		stratified_resim_stats_table <- do.call( "rbind", resim_table_list )
		
	}
	
} # if( length( resim_files ) > 0 )

cat( "1.0. compute_ppc_results_seq_transplant.R main: dim( stratified_resim_stats_table ) = ", dim( stratified_resim_stats_table ), "\n", sep = "\t" )
cat( "1.0. compute_ppc_results_seq_transplant.R main: names( stratified_resim_stats_table ) = ", names( stratified_resim_stats_table ), "\n", sep = "\t" )
# readLines("stdin",n=1)

# Free-up memory:
rm( resim_table_list )
gc()

###############################################################################################
#
# Count number of resims collected;
# 
# Write data outfile "resim_count_file.txt";
#

resim_count_outfile = paste( output_directory, "/", "resim_count_file", ".txt", sep = "" )
sink( file = resim_count_outfile )
# cat( "number of resim files = ", length( resim_file_list ), "\n", sep = "\t" )
cat( "number of resim files = ", n_resim_files, "\n", sep = "\t" )
cat("\n")
# cat( "number of resim tables = ", length( resim_table_list ), "\n", sep = "\t" )
cat( "number of resim tree tables = ", n_resim_tables, "\n", sep = "\t" )
cat("\n")
# cat( "number of resims completed = ", nrow( pooled_resim_table ), "\n", sep = "\t" )
cat( "number of resims completed = ", nrow( stratified_resim_stats_table ), "\n", sep = "\t" )
sink()

###############################################################################################
####################################################################################################
#
# Specify "sumstats_to_include";
#

if( stat_set_name %in% c( "original", "peri_interval_narrower", "peri_interval_wider" ) )
{
	sumstats_to_include = c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             "D_coals_pre","R_coals_pre","D_coals_peri","R_coals_peri",
             "R_max_coals_in_single_clade")
    
    R_only_stat_col_names_vect = c(paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_R",
             "R_coals_pre",
             "R_coals_peri",
             "R_max_coals_in_single_clade")
	
	
}else if( stat_set_name == "pre_interval_divided" )
{
	sumstats_to_include = c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             "DR_coals_pre_1", # new
             "DR_coals_pre_2", # new
             "D_coals_pre_1", # new
             "D_coals_pre_2", # new
             "R_coals_pre_1", # new
             "R_coals_pre_2", # new
             "D_coals_peri",
             "R_coals_peri",
             "R_max_coals_in_single_clade")
	
	sumstats_to_include = c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
            # "DR_coals_pre_1", # new
            # "DR_coals_pre_2", # new
             "D_coals_pre_1", # new
             "D_coals_pre_2", # new
             "R_coals_pre_1", # new
             "R_coals_pre_2", # new
             "D_coals_peri",
             "R_coals_peri",
             "R_max_coals_in_single_clade")
	
    R_only_stat_col_names_vect = c(paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_R",
             "R_coals_pre_1", # new
             "R_coals_pre_2", # new
             "R_coals_peri",
             "R_max_coals_in_single_clade")
	
}else if( stat_set_name == "peri_interval_divided" )
{
	sumstats_to_include = c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             "D_coals_pre",
             "R_coals_pre",
             "DR_coals_peri_1", # new
             "DR_coals_peri_2", # new
             "D_coals_peri_1", # new
             "D_coals_peri_2", # new
             "R_coals_peri_1", # new
             "R_coals_peri_2", # new
             "R_max_coals_in_single_clade")
	
	sumstats_to_include = c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             "D_coals_pre",
             "R_coals_pre",
            # "DR_coals_peri_1", # new
            # "DR_coals_peri_2", # new
             "D_coals_peri_1", # new
             "D_coals_peri_2", # new
             "R_coals_peri_1", # new
             "R_coals_peri_2", # new
             "R_max_coals_in_single_clade")
	
    R_only_stat_col_names_vect = c(paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_R",
             "R_coals_pre",
             "R_coals_peri_1", # new
             "R_coals_peri_2", # new
             "R_max_coals_in_single_clade")
	
}

###############################################################################################

selected_stat_names_vect = sumstats_to_include

####################################################################################################
# This function is used within the function:
# make_stats_residuals_table(...);

is_good_stat_col <- function( stat_col_name, stats_table )
{
	is_good_bool = FALSE
	
	stat_col_vect = stats_table[[ stat_col_name ]]
	
	is_finite_index_vect = which( is.finite( stat_col_vect ) )
	
	if( length( is_finite_index_vect ) > 1 )
	{
		finite_vals_vect = unique( stat_col_vect[ is_finite_index_vect ] )
		
		if( length( finite_vals_vect ) > 1 )
		{
			is_good_bool = TRUE
		}
	}
	
	return( is_good_bool )
}

####################################################################################################
# This function computes the empirical cumulative distribution function (c.d.f.)
# for the values in column "stat_col_name".
#
# This function is used within the function:
# make_stats_residuals_table(...);

calc_ecdf <- function( stat_col_name, stats_table )
{
	stat_col_vect = stats_table[[ stat_col_name ]]
	
	is_finite_index_vect = which( is.finite( stat_col_vect ) )
	
	stat_vect = stat_col_vect[ is_finite_index_vect ]
	
	CDFun <- ecdf( stat_vect )
	
	return( CDFun )
}

####################################################################################################
# This function is used within the function:
# make_stats_residuals_table(...);

est_sim_P_col_vect <- function( stat_col_name, CDFun, stats_table )
{
	n_size = nrow( stats_table )
	P_min = 1 / ( n_size + 1 )
	P_max = n_size / ( n_size + 1 )
	
	stat_col_vect = stats_table[[ stat_col_name ]]
	
	P_col_vect = CDFun( stat_col_vect )
	
	P_col_vect[ P_col_vect <= 0 ] <- P_min
	P_col_vect[ P_col_vect >= 1 ] <- P_max
	
	return( P_col_vect )	
}

####################################################################################################
# This function is used within the function:
# make_stats_residuals_table(...);

est_cumul_P <- function( stat, CDFun, n_size )
{
	P_min = 1 / ( n_size + 1 )
	P_max = n_size / ( n_size + 1 )
	
	cumul_P = CDFun( stat )
	
	if( cumul_P <= 0 ){ cumul_P = P_min }
	if( cumul_P >= 1 ){ cumul_P = P_max }
	
	return( cumul_P )	
}

####################################################################################################
# This function is used within the function:
# make_stats_residuals_table(...);

est_residual_col_vect <- function( cumul_P_col_name, stats_table )
{
	cumul_P_col_vect = stats_table[[ cumul_P_col_name ]]
	
	resids_col_vect = qnorm( cumul_P_col_vect, mean = 0, sd = 1, lower.tail = TRUE )
		
	return( resids_col_vect )	
}

####################################################################################################
####################################################################################################
# This function computes an empirical c.d.f. for each statistic (column of data frame "stats_table");
# and then fills the columns:
# sim_P_col_names_vect;
# # obs_P_col_names_vect;
# sim_residual_col_names_vect;
# # obs_residual_col_names_vect;
# "chisq_sim";
# # "chisq_obs";
# "diff_chisq";
#

fill_chisq_col <- function( stats_table, obs_stats_table, selected_stat_names_vect, R_only_stat_col_names_vect )	
{
	
	diff_chisq_col_vect = rep( 0, nrow( stats_table ) )
	
	if( nrow( stats_table ) > 1 )
	{
		# good_stat_col_bool_vect = sapply( stat_names_vect, function(x){ is_good_stat_col( x, stats_table ) } )
		good_stat_col_bool_vect = sapply( selected_stat_names_vect, function(x){ is_good_stat_col( x, stats_table ) } )
		
		good_stat_col_index_vect = which( good_stat_col_bool_vect )
		
		if( length( good_stat_col_index_vect ) > 0 )
		{
			#
			# For each statistic,
			# compute the empirical c.d.f.;
			#
			
			good_stat_col_names_vect = selected_stat_names_vect[ good_stat_col_index_vect ]
			
			good_stat_col_names_list = as.list( good_stat_col_names_vect )
			
			good_stat_CDF_list = lapply( good_stat_col_names_list, function(x){ calc_ecdf( x, stats_table ) } )
			
			
			#
			# Fill columns:
			# sim_P_col_names_vect;
			# # obs_P_col_names_vect;
			#
			
			sim_P_col_vect_list = Map( function(x,fun){ est_sim_P_col_vect( x, fun, stats_table ) }, good_stat_col_names_list, good_stat_CDF_list )
			
			#
			# Fill columns:
			# sim_residual_col_names_vect;
			# # obs_residual_col_names_vect;
			#
			
			sim_residual_col_vect_list = Map( function(x){ qnorm( x, mean = 0, sd = 1, lower.tail = TRUE ) }, sim_P_col_vect_list )
			
			# Free-up memory:
			rm( sim_P_col_vect_list )
			# gc()
			
			
			#
			# Estimate "obs_P_vect";
			#
			
			n_size = nrow( stats_table )
			
			obs_good_stat_list = as.list( obs_stats_table[ 1, good_stat_col_names_vect ] )
			
			obs_P_list = Map( function(x,fun){ est_cumul_P( x, fun, n_size ) }, obs_good_stat_list, good_stat_CDF_list )
			
			obs_P_vect = unlist( obs_P_list )
			
			#
			# Estimate "obs_residual_vect";
			#
			
			obs_residual_vect = qnorm( obs_P_vect, mean = 0, sd = 1, lower.tail = TRUE )
			
			#
			# Fill columns:
			# "chisq_sim";
			# # "chisq_obs";
			# "diff_chisq";
			#
			
			chisq_obs = sum( sapply( obs_residual_vect, function(x){ x^2 } ) )
			
			sim_resid_sqrd_col_vect_list = Map( function(x){ x^2 }, sim_residual_col_vect_list )
			
			# Free-up memory:
			rm( sim_residual_col_vect_list )
			# gc()
			
			chisq_sim_col_vect = rowSums( as.data.frame( sim_resid_sqrd_col_vect_list ) )
			
			diff_chisq_col_vect = chisq_sim_col_vect - chisq_obs
			
	###############################################################################################
			
			R_only_stat_index_vect = which( good_stat_col_names_list %in% R_only_stat_col_names_vect )
			
			#
			# Fill columns:
			# "chisq_sim";
			# # "chisq_obs";
			# "diff_chisq";
			#
			
			R_only_stat_chisq_obs = sum( sapply( obs_residual_vect[ R_only_stat_index_vect ], function(x){ x^2 } ) )
			
			R_only_stat_chisq_sim_col_vect = rowSums( as.data.frame( sim_resid_sqrd_col_vect_list[ R_only_stat_index_vect ] ) )
			
			
			# Free-up memory:
			rm( sim_resid_sqrd_col_vect_list )
			gc()
			
			R_only_stat_diff_chisq_col_vect = R_only_stat_chisq_sim_col_vect - R_only_stat_chisq_obs
			
	###############################################################################################
			
			
		} # if( length( good_stat_col_index_vect ) > 0 )
		
	} # if( n_trees > 1 )
	
	diff_chisq_table = data.frame( diff_chisq = diff_chisq_col_vect, R_only_stat_diff_chisq = R_only_stat_diff_chisq_col_vect )
	
	# return( stats_table )
	# return( diff_chisq_col_vect )
	return( diff_chisq_table )
	
} # fill_chisq_col <- function( sim_trees_file, post_obs_index, obs_tree, stat_names_vect, selected_stat_names_vect, sim_P_col_names_vect, sim_residual_col_names_vect )	

####################################################################################################
####################################################################################################

write_table_row_to_outfile <- function( sim, sim_table, sim_outfile )
{
	
	if( sim == 1 )
	{
		write.table( sim_table[ sim, ], file = sim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
		# write.table( sim_table[ sim, ], file = sim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
		
	}else
	{
		# write.table( sim_table[ sim, ], file = sim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
		write.table( sim_table[ sim, ], file = sim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
		
	}
	
}

####################################################################################################

####################################################################################################
# 
# Create data frame "ppc_result_table";
# 

seq_resims_per_post_obs = seq( from = 50, to = 1000, by = 50 )

ppc_result_col_names_vect = c( "n_post_obs_represented", "min_resims_per_post_obs", "n_sims", "n_sims_included", "n_reject", "p_reject", "p_formatted" )
ppc_result_table = as.data.frame( matrix( NA, nrow=length( seq_resims_per_post_obs ), ncol=length( ppc_result_col_names_vect) ) )
names( ppc_result_table ) <- ppc_result_col_names_vect

ppc_result_table[[ "min_resims_per_post_obs" ]] = seq_resims_per_post_obs

####################################################################################################
# 
# Create data frame "ppc_result_R_only_stat_table";
# 

ppc_result_R_only_stat_table = ppc_result_table

####################################################################################################
# 
# Write header of data frame "ppc_result_table";
#

ppc_result_outfile = paste( output_directory, "/", "ppc_result_seq_table", "_MODEL_", model_name, "_STAT_SET_", stat_set_name, "_R_and_D_stats", "_DONOR_", pair_ID, ".txt", sep = "" )
# ppc_result_outfile = paste( output_directory, "/", "ppc_result_seq_table", "_MODEL_", model_name, "_STAT_SET_", stat_set_name, "_ABC_METHOD_", ABC_method, "_DONOR_", pair_ID, ".txt", sep = "" )

write.table( ppc_result_table[ c(), ], file = ppc_result_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
# write.table( ppc_result_table[ c(), ], file = ppc_result_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )

####################################################################################################
# 
# Write header of data frame "ppc_result_R_only_stat_table";
#

ppc_result_R_only_stat_outfile = paste( output_directory, "/", "ppc_result_seq_table", "_MODEL_", model_name, "_STAT_SET_", stat_set_name, "_R_only_stats", "_DONOR_", pair_ID, ".txt", sep = "" )

write.table( ppc_result_R_only_stat_table[ c(), ], file = ppc_result_R_only_stat_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
# write.table( ppc_result_table[ c(), ], file = ppc_result_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )

####################################################################################################
# 
# Identify those posterior observations for which we can obtain a "good estimate"
# of the conditional density;
# 

post_obs_levels_vect = unique( stratified_resim_stats_table[[ "post_obs_index" ]] )

post_obs_levels_list = as.list( post_obs_levels_vect )

post_obs_index_vect_list = Map( function(x){ which( stratified_resim_stats_table[[ "post_obs_index" ]] == x ) }, post_obs_levels_list )

post_obs_resim_count_list = lapply( post_obs_index_vect_list , length )
post_obs_resim_count_vect = unlist( post_obs_resim_count_list )


####################################################################################################
# 
# For debugging!!!
#

debug_col_names_vect = c( "post_obs_index", "min_resims_per_post_obs", "n_sims", selected_stat_names_vect )
debug_table = as.data.frame( matrix( NA, nrow=1, ncol=length( debug_col_names_vect) ) )
names( debug_table ) <- debug_col_names_vect

# 
# Write header of data frame "ppc_result_table";
#

debug_outfile = paste( output_directory, "/", "debug_table", "_MODEL_", model_name, "_STAT_SET_", stat_set_name, "_DONOR_", pair_ID, ".txt", sep = "" )
# debug_outfile = paste( output_directory, "/", "debug_table", "_MODEL_", model_name, "_STAT_SET_", stat_set_name, "_ABC_METHOD_", ABC_method, "_DONOR_", pair_ID, ".txt", sep = "" )

write.table( debug_table[ c(), ], file = debug_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
# write.table( debug_table[ c(), ], file = debug_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )


####################################################################################################
# 
# Fill data frames "ppc_result_table";
#

for( row_index in 1:nrow( ppc_result_table ) )
{

	min_resims_per_post_obs = ppc_result_table[[ "min_resims_per_post_obs" ]][ row_index ]

	good_est_bool_vect = c( post_obs_resim_count_vect >= min_resims_per_post_obs )
	good_est_index_vect = which( good_est_bool_vect )

	good_est_post_obs_index_vect_list = post_obs_index_vect_list[ good_est_index_vect ]
	
####################################################################################################
	# 
	# Fill data frames "stratified_resim_stats_table" columns:
	# sim_P_col_names_vect;
	# sim_residual_col_names_vect;
	# "chisq_sim";
	# "diff_chisq";
	# 
	
	diff_chisq_table_list = Map( function(x){ fill_chisq_col( stratified_resim_stats_table[ x, ], obs_stats_table, selected_stat_names_vect, R_only_stat_col_names_vect ) }, good_est_post_obs_index_vect_list )
	
	# diff_chisq_col_vect = unlist( diff_chisq_vect_list )
	diff_chisq_table <- do.call( "rbind", diff_chisq_table_list )
	
	# Free-up memory:
	# rm( stratified_resim_stats_table )
	# gc()
	
####################################################################################################
	# 
	# Create data frame "ppc_result_table";
	# 
	
	n_post_obs_good_est = length( post_obs_resim_count_vect )

	#
	# Estimate (Bayesian) p-value:
	#

	diff_chisq_col_vect = diff_chisq_table[[ "diff_chisq" ]]
	
	n_post_obs_good_est = length( good_est_index_vect )
	
	n_sims = sum( post_obs_resim_count_vect )
	
	n_sims_included = length( diff_chisq_col_vect )

	n_reject = length( which( diff_chisq_col_vect > 0 ) )

	# p_reject = n_reject / n_sims
	p_reject = n_reject / n_sims_included

	p_formatted = round( p_reject, digits = 4 )

	ppc_result_table[[ "n_post_obs_represented" ]][ row_index ] = n_post_obs_good_est
	# ppc_result_table[[ "min_resims_per_post_obs" ]][ row_index ] = min_resims_per_post_obs
	ppc_result_table[[ "n_sims" ]][ row_index ] = n_sims
	ppc_result_table[[ "n_sims_included" ]][ row_index ] = n_sims_included
	ppc_result_table[[ "n_reject" ]][ row_index ] = n_reject
	ppc_result_table[[ "p_reject" ]][ row_index ] = p_reject
	ppc_result_table[[ "p_formatted" ]][ row_index ] = p_formatted
	
	
	# 
	# Write row (row_index") of data frame "ppc_result_table";
	#
	
	# write.table( ppc_result_table[ row_index, ], file = ppc_result_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	write.table( ppc_result_table[ row_index, ], file = ppc_result_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
	
	###############################################################################################
	
	#
	# Estimate (Bayesian) p-value:
	#

	# diff_chisq_col_vect = diff_chisq_table[[ "diff_chisq" ]]
	diff_chisq_col_vect = diff_chisq_table[[ "R_only_stat_diff_chisq" ]]
	
	n_post_obs_good_est = length( good_est_index_vect )
	
	n_sims = sum( post_obs_resim_count_vect )
	
	n_sims_included = length( diff_chisq_col_vect )

	n_reject = length( which( diff_chisq_col_vect > 0 ) )

	# p_reject = n_reject / n_sims
	p_reject = n_reject / n_sims_included

	p_formatted = round( p_reject, digits = 4 )

	ppc_result_R_only_stat_table[[ "n_post_obs_represented" ]][ row_index ] = n_post_obs_good_est
	# ppc_result_R_only_stat_table[[ "min_resims_per_post_obs" ]][ row_index ] = min_resims_per_post_obs
	ppc_result_R_only_stat_table[[ "n_sims" ]][ row_index ] = n_sims
	ppc_result_R_only_stat_table[[ "n_sims_included" ]][ row_index ] = n_sims_included
	ppc_result_R_only_stat_table[[ "n_reject" ]][ row_index ] = n_reject
	ppc_result_R_only_stat_table[[ "p_reject" ]][ row_index ] = p_reject
	ppc_result_R_only_stat_table[[ "p_formatted" ]][ row_index ] = p_formatted
	
	
	# 
	# Write row (row_index") of data frame "ppc_result_R_only_stat_table";
	#
	
	# write.table( ppc_result_R_only_stat_table[ row_index, ], file = ppc_result_R_only_stat_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	write.table( ppc_result_R_only_stat_table[ row_index, ], file = ppc_result_R_only_stat_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
	
	###############################################################################################
	
	
} # for( row_index in 1:nrow( ppc_result_table ) )

####################################################################################################

cat( "done!", "\n" )




















