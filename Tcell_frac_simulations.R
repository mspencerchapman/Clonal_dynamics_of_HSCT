library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(phytools)
library(rsimpop)
library(abc)

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

root_dir<-ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/Clonal_dynamics_of_HSCT","/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT")
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch119/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
source(paste0(root_dir,"/data/HSCT_functions.R"))
setwd(tree_mut_dir); source("treemut.R");setwd(root_dir)
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")

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


###-----------------Pair metadata & trees-----------------

#Read in Pair metadata data frame
exp_nos<-c(11,13,21,24,25,28,31,38,40,41)
Pair_metadata=data.frame(Pair=paste0("Pair",exp_nos),
                         Pair_new=factor(c("Pair_9","Pair_7","Pair_5","Pair_1","Pair_4","Pair_10","Pair_6","Pair_8","Pair_3","Pair_2"),levels=paste0("Pair_",1:10)),
                         Age=c(74.8,65.5,64.5,34.2,58.4,79.9,65.2,65.8, 51.9,42.4),
                         Age_at_transplant=c(66,36,50,18,47,63,43,35,35,30),
                         MNC_dose=c(2.66,4.17,16.28,NA,10.94,13.9,NA,4.05,2.48,15.01),
                         CD34_dose=c(1.56,NA,7.1,7.9,8.97,2.4,NA,NA,NA,4.51),
                         stem_cell_source=c("BM","BM","PBSC","PBSC","PBSC","PBSC","BM","BM","BM","PBSC"),
                         conditioning=c("MAC","MAC","RIC","MAC","RIC","RIC","MAC","MAC","MAC","MAC"))
new_pair_names=paste("Pair",1:nrow(Pair_metadata),sep = "_") #Useful named vector to swap from old to new names
names(new_pair_names)<-Pair_metadata%>%arrange(Age)%>%pull(Pair)

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#-------------IMPORT PARAMETERS AND SUMMARY STATISTICS FOR THE SIMULATIONS & PERFORM ABC ------------------------------

setwd(paste0(root_dir,"/data/ABC_simulation_results/Tcell_fraction_simulations/"))

#setwd("/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Tcell_duration")

#Now extract the parameters/ summary statistics
all_stats<-lapply(Pair_metadata%>%pull(Pair),function(pair){
  cat(pair,sep="\n")
  
  sim.files=list.files(path=pair,pattern=".Rds",full.names = T)
  params_file=paste0("all_params_",pair,".Rds")
  sumstats_file=paste0("all_sumstats_",pair,".Rds")
  
  if(file.exists(params_file) & file.exists(sumstats_file)){
    cat("Reading in saved files",sep="\n\n")
    params<-readRDS(params_file)
    sumstats<-readRDS(sumstats_file)
  } else {
    cat("Importing individual results files and concatenating",sep="\n\n")
    sim.output<-lapply(sim.files,function(file) {if(which(sim.files==file)%%100 == 0) {print(which(sim.files==file))};readRDS(file)})
    
    params<-lapply(sim.output,function(res) unlist(res$params))%>%dplyr::bind_rows()%>%mutate(idx=1:length(sim.output))
    sumstats<-lapply(1:length(sim.output),function(i) {
      if(i%%100 == 0) {print(i)}
      df<-sim.output[[i]]$ratios%>%mutate(idx=i)
      return(df)
      })%>%dplyr::bind_rows()
    
    #temp<-lapply(1:length(sim.output),function(i) if(is.null(sim.output[[i]]$ratios)){print(i)})
    
    saveRDS(params,file=params_file)
    saveRDS(sumstats,file=sumstats_file)
  }
  
  df<-left_join(params,sumstats,by="idx")%>%mutate(Pair=pair)
  return(df)
})%>%bind_rows()

Tcell_fraction_simulation_plot<-all_stats%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  filter(Pair%in%paste("Pair",3:10,sep="_") & Myeloid>0.05)%>%
  ggplot(aes(x=T_cell_clone_duration,y=`T/M_ratio`))+
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method="lm",col="black",fullrange=T)+
  theme_bw()+
  my_theme+
  facet_wrap(~Pair,nrow=2)+
  geom_hline(yintercept=1,linetype=2)+
  annotate(geom="rect", xmin = 8, xmax = 15, ymin = 0, ymax = 1.2,alpha = .2)+
  scale_x_continuous(limits=c(0,20))+
  scale_y_continuous(limits=c(0,1.2),breaks=seq(0,2,0.2))+
  labs(x="Simulated lifespan of T cell clones\n(Average T cell age at sampling is half of this)",
       y="Expanded clone T-cell fraction/\nExpanded clone myeloid fraction")

ggsave(filename = paste0(plots_dir,"Tcell_fraction_simulation_plot.pdf"),Tcell_fraction_simulation_plot,device = "pdf",width=7,height = 3)

