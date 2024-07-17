## Run on farm5
args = commandArgs(TRUE)
args

source_directory = toString(args[1])
sim_directory = toString(args[2])
output_directory = toString(args[3])
posterior_sample_file_name = toString(args[4])
obs_stats_directory = toString(args[5])
obs_stats_file_name = toString(args[6])
pooled_sim_trees_directory = toString(args[7])
n_sims_accept = as.integer(args[8])
n_sims_max = as.integer(args[9])
model_name = toString(args[10])
stat_set_name = toString(args[11])
ABC_method = toString(args[12])
pair_ID = toString(args[13])

cat( "pair_ID = ", pair_ID, "\n", sep = "\t" )
# readLines("stdin",n=1)

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

my_working_directory = output_directory
setwd(my_working_directory)

system("mkdir pdfs")

####################################################################################################
####################################################################################################
# This function maps the values of a parameter "theta" (restricted to range [a,b])
# onto the real line, using the "logit" function:
#
# z = logit( p ) = ln( p/(1-p) )
#

param_to_z <- function( theta, a, b, epsilon )
{
	
	p = ( theta - a ) / ( b - a )
	
	if( p <= 0 )
	{
		p = epsilon
		
	}else if( 1 <= p )
	{
		p = 1 - epsilon
		
	}
	
	z = log(p) - log(1-p)
	
	return( z )
	
}

####################################################################################################
# This function maps the values of a variable "z" (anywhere on the real line)
# onto the parameter range [a,b], using the "logistic" function:
#
# p = exp(z)/(1+exp(z)) = 1/(1+exp(-z))
#

z_to_param <- function( z, a, b, epsilon )
{
	
	p = 1/( 1 + exp(-z) )
	
	if( p <= 0 )
	{
		p = epsilon
		
	}else if( 1 <= p )
	{
		p = 1 - epsilon
		
	}
	
	theta = a + ( ( b - a ) * p )
	
	return( theta )
	
}

####################################################################################################
epsilon = 0.0001
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

# Save data frame "obs_stats_table" to directory "obs_stats_directory";

obs_stats_file = paste0( obs_stats_directory, "/", obs_stats_file_name )
write.table( obs_stats_table, file = obs_stats_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

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

sim_file_list = list.files( path = sim_directory, pattern = "sim_" )

cat( "0.0. main: length( sim_file_list ) = ", length( sim_file_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)

setwd( sim_directory )
sim_table_list = list()
for( i in 1:length( sim_file_list ) )
{
	sim_file = sim_file_list[[ i ]]
	
	if( file_test( "-f", sim_file ) )
	{
		sim_result = try( sim_table <- read.table( file = sim_file, sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE ) )
		if( class( sim_result ) != "try-error" )
		{
			sim_table_list = append( sim_table_list, list( sim_table ) )
			
		}
		
	}
	
}
n_sim_tables = length( sim_table_list )
cat( "0.0. main: length( sim_table_list ) = ", length( sim_table_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)

if( length( sim_table_list ) > 0 ) 
{
	if( length( sim_table_list ) == 1 ) 
	{
		# pooled_sim_table = sim_table_list[[ 1 ]]
		sim_table = sim_table_list[[ 1 ]]
		
	}else
	{
		# pooled_sim_table <- do.call( "rbind", sim_table_list )
		sim_table <- do.call( "rbind", sim_table_list )
		
	}
	
} # if( length( sim_files ) > 0 )

cat( "1.0. compute_conditionals_drivers.R main: dim( sim_table ) = ", dim( sim_table ), "\n", sep = "\t" )
cat( "1.0. compute_conditionals_drivers.R main: names( sim_table ) = ", names( sim_table ), "\n", sep = "\t" )
# readLines("stdin",n=1)

# Free-up memory:
rm( sim_table_list )
gc()


###############################################################################################
#
# Count number of sims collected;
# 
# Write data outfile "sim_count_file.txt";
#

sim_count_outfile = paste( output_directory, "/", "sim_count_file", ".txt", sep = "" )
sink( file = sim_count_outfile )
cat( "number of sim tree files = ", length( sim_file_list ), "\n", sep = "\t" )
cat("\n")
# cat( "number of sim tree tables = ", length( sim_table_list ), "\n", sep = "\t" )
cat( "number of sim tree tables = ", n_sim_tables, "\n", sep = "\t" )
cat("\n")
cat( "number of sim trees completed = ", nrow( sim_table ), "\n", sep = "\t" )
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
	
}

###############################################################################################
#
# Remove simulations which have missing data;
#

abc_cols = c( param_names_vect, sumstats_to_include )

complete_rows_index_vect = complete.cases( sim_table[ , abc_cols ] )
sim_table = sim_table[ complete_rows_index_vect, ]

cat( "1.1. run_ABC_drift.R main: dim( sim_table ) = ", dim( sim_table ), "\n", sep = "\t" )
cat( "1.1. run_ABC_drift.R main: names( sim_table ) = ", names( sim_table ), "\n", sep = "\t" )
# readLines("stdin",n=1)

if( nrow( sim_table ) > n_sims_max )
{
  sim_table = sim_table[ 1:n_sims_max, ]
}
n_sims = nrow( sim_table ) # n_sims = nrow( parameters.sim )
tol = n_sims_accept / n_sims

####################################################################################################
#
# Specify input to abc(...) function;
#
# Transform parameters:
#

parameters.sim=sim_table[ , param_names_vect ]


HSC_pop_size_min = min( parameters.sim[[ "HSC_pop_size" ]] )
HSC_pop_size_max = max( parameters.sim[[ "HSC_pop_size" ]] )

parameters.sim[[ "HSC_pop_size" ]] = sapply( parameters.sim[[ "HSC_pop_size" ]], function(x){ param_to_z( x, a = HSC_pop_size_min, b = HSC_pop_size_max, epsilon ) } )

HSCT_bottleneck_min = min( parameters.sim[[ "HSCT_bottleneck" ]] )
HSCT_bottleneck_max = max( parameters.sim[[ "HSCT_bottleneck" ]] )

parameters.sim[[ "HSCT_bottleneck" ]] = sapply( parameters.sim[[ "HSCT_bottleneck" ]], function(x){ param_to_z( x, a = HSCT_bottleneck_min, b = HSCT_bottleneck_max, epsilon ) } )

number_drivers_per_year_min = min( parameters.sim[[ "number_drivers_per_year" ]] )
number_drivers_per_year_max = max( parameters.sim[[ "number_drivers_per_year" ]] )
# number_drivers_per_year_min = 1
# number_drivers_per_year_max = 1000

parameters.sim[[ "number_drivers_per_year" ]] = sapply( parameters.sim[[ "number_drivers_per_year" ]], function(x){ param_to_z( x, a = number_drivers_per_year_min, b = number_drivers_per_year_max, epsilon ) } )

gamma_shape_min = min( parameters.sim[[ "gamma_shape" ]] )
gamma_shape_max = max( parameters.sim[[ "gamma_shape" ]] )
# gamma_shape_min = 0.1
# gamma_shape_max = 4

parameters.sim[[ "gamma_shape" ]] = sapply( parameters.sim[[ "gamma_shape" ]], function(x){ param_to_z( x, a = gamma_shape_min, b = gamma_shape_max, epsilon ) } )

gamma_rate_min = min( parameters.sim[[ "gamma_rate" ]] )
gamma_rate_max = max( parameters.sim[[ "gamma_rate" ]] )
# gamma_rate_min = 5
# gamma_rate_max = 120

parameters.sim[[ "gamma_rate" ]] = sapply( parameters.sim[[ "gamma_rate" ]], function(x){ param_to_z( x, a = gamma_rate_min, b = gamma_rate_max, epsilon ) } )

dpcpd_min = min( parameters.sim[[ "dpcpd" ]] )
dpcpd_max = max( parameters.sim[[ "dpcpd" ]] )
# dpcpd_max = 1e-5

parameters.sim[[ "dpcpd" ]] = sapply( parameters.sim[[ "dpcpd" ]], function(x){ param_to_z( x, a = dpcpd_min, b = dpcpd_max, epsilon ) } )


if( model_name == "m2" )
{
	
	gamma_shape_engraftment_min = 0.1
	gamma_shape_engraftment_max = 1.5
	
	parameters.sim[[ "gamma_shape_engraftment" ]] = sapply( parameters.sim[[ "gamma_shape_engraftment" ]], function(x){ param_to_z( x, a = gamma_shape_engraftment_min, b = gamma_shape_engraftment_max, epsilon ) } )
	
	gamma_rate_engraftment_min = 0.1
	gamma_rate_engraftment_max = 1.5
	
	parameters.sim[[ "gamma_rate_engraftment" ]] = sapply( parameters.sim[[ "gamma_rate_engraftment" ]], function(x){ param_to_z( x, a = gamma_rate_engraftment_min, b = gamma_rate_engraftment_max, epsilon ) } )
	
}else if( model_name == "m3" )
{
	
	exaggeration_ratio_min = 1.5
	exaggeration_ratio_max = 6.0
	
	parameters.sim[[ "exaggeration_ratio" ]] = sapply( parameters.sim[[ "exaggeration_ratio" ]], function(x){ param_to_z( x, a = exaggeration_ratio_min, b = exaggeration_ratio_max, epsilon ) } )
	
	exaggeration_prop_min = 0.1
	exaggeration_prop_max = 0.3
	
	parameters.sim[[ "exaggeration_prop" ]] = sapply( parameters.sim[[ "exaggeration_prop" ]], function(x){ param_to_z( x, a = exaggeration_prop_min, b = exaggeration_prop_max, epsilon ) } )
	
}

####################################################################################################
#
# Specify input to abc(...) function;
#
# Select summary statistics:
#

# sumstat=summary_stats.sim[,sumstats_to_include]
sumstat=sim_table[ , sumstats_to_include ]

cat( "2.0. run_ABC_drift.R main: dim( sumstat ) = ", dim( sumstat ), "\n", sep = "\t" )
cat( "2.0. run_ABC_drift.R main: names( sumstat ) = ", names( sumstat ), "\n", sep = "\t" )
# readLines("stdin",n=1)

# target=stats_data[sumstats_to_include]
target = unlist( obs_stats_table[ 1, sumstats_to_include ] )

cat( "3.0. run_ABC_drift.R main: length( target ) = ", length( target ), "\n", sep = "\t" )
cat( "3.0. run_ABC_drift.R main: names( target ) = ", names( target ), "\n", sep = "\t" )
# readLines("stdin",n=1)

abc_out = abc(target = target, param = parameters.sim, sumstat = sumstat, tol = tol, method = ABC_method)

####################################################################################################

dist_table = as.data.frame(abc_out$dist)
weights_table = as.data.frame(abc_out$weights)
# stats_accept_table = as.data.frame(abc_out$ss)
unadj_z_table = as.data.frame(abc_out$unadj.values)

if( ABC_method %in% c( "loclinear", "ridge", "neuralnet" ) )
{
	adj_z_table = as.data.frame(abc_out$adj.values)
	residuals_table = as.data.frame(abc_out$residuals)
	
}

if( ABC_method %in% c( "loclinear", "ridge", "neuralnet" ) )
{
	unadj_posterior_sample_table = as.data.frame(abc_out$unadj.values)
	posterior_sample_table = as.data.frame(abc_out$adj.values)
	
}else # ABC_method == "rejection"
{
	posterior_sample_table = as.data.frame(abc_out$unadj.values)
	
}

####################################################################################################
#
# Inverse transform parameters:
#

posterior_sample_table[[ "HSC_pop_size" ]] = sapply( posterior_sample_table[[ "HSC_pop_size" ]], function(x){ z_to_param( x, a = HSC_pop_size_min, b = HSC_pop_size_max, epsilon ) } )

posterior_sample_table[[ "HSCT_bottleneck" ]] = sapply( posterior_sample_table[[ "HSCT_bottleneck" ]], function(x){ z_to_param( x, a = HSCT_bottleneck_min, b = HSCT_bottleneck_max, epsilon ) } )

posterior_sample_table[[ "number_drivers_per_year" ]] = sapply( posterior_sample_table[[ "number_drivers_per_year" ]], function(x){ z_to_param( x, a = number_drivers_per_year_min, b = number_drivers_per_year_max, epsilon ) } )

posterior_sample_table[[ "gamma_shape" ]] = sapply( posterior_sample_table[[ "gamma_shape" ]], function(x){ z_to_param( x, a = gamma_shape_min, b = gamma_shape_max, epsilon ) } )

posterior_sample_table[[ "gamma_rate" ]] = sapply( posterior_sample_table[[ "gamma_rate" ]], function(x){ z_to_param( x, a = gamma_rate_min, b = gamma_rate_max, epsilon ) } )

posterior_sample_table[[ "dpcpd" ]] = sapply( posterior_sample_table[[ "dpcpd" ]], function(x){ z_to_param( x, a = dpcpd_min, b = dpcpd_max, epsilon ) } )

if( model_name == "m2" )
{
	posterior_sample_table[[ "gamma_shape_engraftment" ]] = sapply( posterior_sample_table[[ "gamma_shape_engraftment" ]], function(x){ z_to_param( x, a = gamma_shape_engraftment_min, b = gamma_shape_engraftment_max, epsilon ) } )
	posterior_sample_table[[ "gamma_rate_engraftment" ]] = sapply( posterior_sample_table[[ "gamma_rate_engraftment" ]], function(x){ z_to_param( x, a = gamma_rate_engraftment_min, b = gamma_rate_engraftment_max, epsilon ) } )
	
}else if( model_name == "m3" )
{
	posterior_sample_table[[ "exaggeration_ratio" ]] = sapply( posterior_sample_table[[ "exaggeration_ratio" ]], function(x){ z_to_param( x, a = exaggeration_ratio_min, b = exaggeration_ratio_max, epsilon ) } )
	posterior_sample_table[[ "exaggeration_prop" ]] = sapply( posterior_sample_table[[ "exaggeration_prop" ]], function(x){ z_to_param( x, a = exaggeration_prop_min, b = exaggeration_prop_max, epsilon ) } )
	
}

# posterior_sample_file = paste0( output_directory, "/", "posterior_sample.txt" )
posterior_sample_file = paste0( output_directory, "/", posterior_sample_file_name )
write.table( posterior_sample_table, file = posterior_sample_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )


if( ABC_method %in% c( "loclinear", "ridge", "neuralnet" ) )
{
	unadj_posterior_sample_table[[ "HSC_pop_size" ]] = sapply( unadj_posterior_sample_table[[ "HSC_pop_size" ]], function(x){ z_to_param( x, a = HSC_pop_size_min, b = HSC_pop_size_max, epsilon ) } )
	
	unadj_posterior_sample_table[[ "HSCT_bottleneck" ]] = sapply( unadj_posterior_sample_table[[ "HSCT_bottleneck" ]], function(x){ z_to_param( x, a = HSCT_bottleneck_min, b = HSCT_bottleneck_max, epsilon ) } )
	
	unadj_posterior_sample_table[[ "number_drivers_per_year" ]] = sapply( unadj_posterior_sample_table[[ "number_drivers_per_year" ]], function(x){ z_to_param( x, a = number_drivers_per_year_min, b = number_drivers_per_year_max, epsilon ) } )
	
	unadj_posterior_sample_table[[ "gamma_shape" ]] = sapply( unadj_posterior_sample_table[[ "gamma_shape" ]], function(x){ z_to_param( x, a = gamma_shape_min, b = gamma_shape_max, epsilon ) } )
	
	unadj_posterior_sample_table[[ "gamma_rate" ]] = sapply( unadj_posterior_sample_table[[ "gamma_rate" ]], function(x){ z_to_param( x, a = gamma_rate_min, b = gamma_rate_max, epsilon ) } )
	
	unadj_posterior_sample_table[[ "dpcpd" ]] = sapply( unadj_posterior_sample_table[[ "dpcpd" ]], function(x){ z_to_param( x, a = dpcpd_min, b = dpcpd_max, epsilon ) } )
	
	if( model_name == "m2" )
	{
		unadj_posterior_sample_table[[ "gamma_shape_engraftment" ]] = sapply( unadj_posterior_sample_table[[ "gamma_shape_engraftment" ]], function(x){ z_to_param( x, a = gamma_shape_engraftment_min, b = gamma_shape_engraftment_max, epsilon ) } )
		unadj_posterior_sample_table[[ "gamma_rate_engraftment" ]] = sapply( unadj_posterior_sample_table[[ "gamma_rate_engraftment" ]], function(x){ z_to_param( x, a = gamma_rate_engraftment_min, b = gamma_rate_engraftment_max, epsilon ) } )
		
	}else if( model_name == "m3" )
	{
		unadj_posterior_sample_table[[ "exaggeration_ratio" ]] = sapply( unadj_posterior_sample_table[[ "exaggeration_ratio" ]], function(x){ z_to_param( x, a = exaggeration_ratio_min, b = exaggeration_ratio_max, epsilon ) } )
		unadj_posterior_sample_table[[ "exaggeration_prop" ]] = sapply( unadj_posterior_sample_table[[ "exaggeration_prop" ]], function(x){ z_to_param( x, a = exaggeration_prop_min, b = exaggeration_prop_max, epsilon ) } )
		
	}
	
	unadj_posterior_sample_file = paste0( output_directory, "/", "unadj_posterior_sample", ".txt" )
	write.table( unadj_posterior_sample_table, file = unadj_posterior_sample_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	
}

####################################################################################################

dist_file = paste0( output_directory, "/", "dist", ".txt" )
write.table( dist_table, file = dist_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

weights_file = paste0( output_directory, "/", "weights", ".txt" )
write.table( weights_table, file = weights_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

unadj_z_file = paste0( output_directory, "/", "unadj_z", ".txt" )
write.table( unadj_z_table, file = unadj_z_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

if( ABC_method %in% c( "loclinear", "ridge", "neuralnet" ) )
{
	adj_z_file = paste0( output_directory, "/", "adj_z", ".txt" )
	write.table( adj_z_table, file = adj_z_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	
	residuals_file = paste0( output_directory, "/", "residuals", ".txt" )
	write.table( residuals_table, file = residuals_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	
}

###############################################################################################
#
# Create data frame "prior_sample_table";
#

n_prior_sample = min( 10*n_sims_accept, nrow( sim_table ) )

row_index_vect = sample( x = nrow( sim_table ), size = n_prior_sample, replace = FALSE )

prior_sample_table = sim_table[ row_index_vect, param_names_vect ]

# posterior_sample_file = paste0( output_directory, "/", "posterior_sample.txt" )
prior_sample_file = paste0( output_directory, "/", "prior_sample.txt" )
write.table( prior_sample_table, file = prior_sample_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

###############################################################################################
#
# Create data frame "posterior_and_prior_table";
#

posterior_table = posterior_sample_table

prior_table = prior_sample_table

posterior_table$type="post"
prior_table$type="prior"
posterior_and_prior_table = rbind( posterior_table, prior_table )

posterior_and_prior_file = paste0( output_directory, "/", "posterior_and_prior", ".txt" )
write.table( posterior_and_prior_table, file = posterior_and_prior_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

###############################################################################################


cat( "done!", "\n" )


































