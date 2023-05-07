##Examine the NULL and selection models for Type 1 and Type 2 selection parameters
library(stringr)
library(ape)
library(seqinr)
library(data.table)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)


root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch119/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
source(paste0(root_dir,"/data/HSCT_functions.R"))
setwd(tree_mut_dir); source("treemut.R");setwd(root_dir)
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

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


####---------------------SELECTION TYPE ANALYSIS--------------------
exp_nos<-c(11,13,21,24,25,28,31,38,40,41)
Pair_metadata=data.frame(Pair=paste0("Pair",exp_nos),
                         Pair_new=factor(c("Pair_9","Pair_7","Pair_5","Pair_1","Pair_4","Pair_10","Pair_6","Pair_8","Pair_3","Pair_2"),levels=paste0("Pair_",1:10)),
                         Age=c(74.8,65.5,64.5,34.2,58.4,79.9,65.2,65.8, 51.9,42.4),
                         Age_at_transplant=c(66,36,50,18,47,63,43,35,35,30),
                         MNC_dose=c(2.66,4.17,16.28,NA,10.94,13.9,NA,4.05,2.48,15.01),
                         CD34_dose=c(1.56,NA,7.1,7.9,8.97,2.4,NA,NA,NA,4.51),
                         stem_cell_source=c("BM","BM","PBSC","PBSC","PBSC","PBSC","BM","BM","BM","PBSC"),
                         conditioning=c("MAC","MAC","RIC","MAC","RIC","RIC","MAC","MAC","MAC","MAC"))

NULL_model_folder<-"/lustre/scratch119/realdata/mdt1/team154/ms56/Zur_HSCT/ABC_models/ABC_Dec2022"
Type1_model_folder<-"/lustre/scratch119/realdata/mdt1/team154/ms56/Zur_HSCT/ABC_models/ABC_Engraftment_selection"
Type2_model_folder<-"/lustre/scratch119/realdata/mdt1/team154/ms56/Zur_HSCT/ABC_models/ABC_Postengraftment_selection"

###--------------Get summary stats from the data-----------------
ABC.trees<-readRDS("/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/trees_for_ABC.Rds")

##EXTRACT SUMMARY STATISTICS
#(1) 3 LARGEST CLADES
all.ss<-Map(tree=ABC.trees,pair=names(ABC.trees),function(tree,pair) {
  donor_id<-get_DR_ids(tree)['donor_ID']
  recip_id<-get_DR_ids(tree)['recip_ID']
  
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
  
  peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
  DR_ltt_peri<-get_ltt(tree=tree.time,time_points=peri_HSCT_time_points)
  D_ltt_peri<-get_ltt(tree=D_tree.time,time_points=peri_HSCT_time_points)
  R_ltt_peri<-get_ltt(tree=R_tree.time,time_points=peri_HSCT_time_points)
  DR_coals_peri<-get_coalescences(DR_ltt_peri)
  D_coals_peri<-get_coalescences(D_ltt_peri)
  R_coals_peri<-get_coalescences(R_ltt_peri)
  
  #(4) Pre-HSCT LTT/ coalescences
  pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
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
  ss_names=c(paste("largest_clade_D",1:3,sep="_"),
             paste("largest_clade_R",1:3,sep="_"),
             "n_singletons_D",
             "n_singletons_R",
             paste(c("DR","D","R"),"coals_pre",sep="_"),
             paste(c("DR","D","R"),"coals_peri",sep="_"),
             "R_max_coals_in_single_clade",
             "mean_abs_log2FC",
             "median_abs_log2FC",
             "max_abs_log2FC",
             "mean_abs_change",
             "median_abs_change",
             "max_abs_change")
  
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
  return(ss_comb)
})

###---------------------Run the ABC for each model to get the uids of trees in the posterior distribution----------------

library(abc)

selected_uids_by_model<-Map(model_dir=list(NULL_model_folder,Type1_model_folder,Type2_model_folder),f=function(model_dir) {
  setwd(model_dir)
  cat(model_dir,sep="\n")
  
  model_selected_uids<-lapply(names(ABC.trees),function(pair){
    cat(pair,sep="\n")
    
    sim.files=list.files(path=pair,pattern=".Rds",full.names = T)
    params_file=paste0("all_params_",pair,".Rds")
    sumstats_file=paste0("all_sumstats_",pair,".Rds")
    
    if(file.exists(params_file) & file.exists(sumstats_file)){
      params<-readRDS(params_file)
      sumstats<-readRDS(sumstats_file)
    } else {
      sim.output<-lapply(sim.files,function(file) {if(which(sim.files==file)%%100 == 0) {print(which(sim.files==file))};readRDS(file)})
      
      sim.file.uids<-stringr::str_sub(sim.files,-36,-5)
      
      params<-lapply(sim.output,function(res) unlist(res$params))%>%dplyr::bind_rows()%>%bind_cols(data.frame(uid=sim.file.uids))
      sumstats<-lapply(sim.output,function(res) unlist(res$sumstats))%>%dplyr::bind_rows()
      
      saveRDS(params,file=params_file)
      saveRDS(sumstats,file=sumstats_file)
    }
    
    #Can run with different sets of summary stats
    all_tree_sumstats=1:15 #All the phylo statistics (not including targeted sequencing results)
    tree_sumstats_no_Rpre_coals=c(1,2,3,4,6,8,9,10,11,12,13,14)
    
    #For assessing engrafting HSC number ignore the R pre-HSCT coalescences (which distorts the posterior as ABC is unable to fit this)
    which_sumstats<-tree_sumstats_no_Rpre_coals
    
    #Run the ABC - log transform the parameters for the regression step. Use tolerance of 5% (= ~ top 1000 for 20,000 simulations)
    abc.res<-abc(target = all.ss[[pair]][which_sumstats],param = params[,6:7],sumstat=sumstats[,which_sumstats],tol = 0.05,transf = c("log","log"),method="neuralnet")
    
    selected_uids<-params$uid[abc.res$region]
    
    data.ss.df<-data.frame(sumstat=names(all.ss[[pair]][which_sumstats]),value=all.ss[[pair]][which_sumstats])
    post.params<-as.data.frame(abc.res$adj.values)%>%gather(key="parameter",value="value")%>%mutate(PairID=pair,.before=1)
    
    return(selected_uids)
  })
  
  return(model_selected_uids)
})









#########


##Import all the trees from the different simulation models  
all.sim.trees<-Map(Pair_ID=Pair_metadata$Pair,
                   NULL_uids=selected_uids_by_model[[1]],
                   Type1_uids=selected_uids_by_model[[2]],
                   Type2_uids=selected_uids_by_model[[3]],
                   f=function(Pair_ID,NULL_uids,Type1_uids,Type2_uids) {
  cat(Pair_ID,sep="\n")
  
  NULL_model_files<-paste0(NULL_model_folder,"/",Pair_ID,"/tree_files/tree_",Pair_ID,"_",NULL_uids,".tree")
  cat(paste("Importing",length(NULL_model_files),"NULL model trees (no HSCT-specific selection)"),sep="\n")
  NULL_model_trees<-lapply(NULL_model_files,read.tree)
  
  Type1_model_files<-paste0(Type1_model_folder,"/",Pair_ID,"/tree_files/tree_",Pair_ID,"_",Type1_uids,".tree")
  cat(paste("Importing",length(Type1_model_files),"Type 1 model trees (Engraftment selection)"),sep="\n")
  Type1_model_trees<-lapply(Type1_model_files,read.tree)
  
  Type2_model_files<-paste0(Type2_model_folder,"/",Pair_ID,"/tree_files/tree_",Pair_ID,"_",Type2_uids,".tree")
  cat(paste("Importing",length(Type2_model_files),"Type 2 model trees (Post-engraftment selection)"),sep="\n")
  Type2_model_trees<-lapply(Type2_model_files,read.tree)
  
  return(list(NULL_model_trees=NULL_model_trees,
              Type1_model_trees=Type1_model_trees,
              Type2_model_trees=Type2_model_trees))
})
names(all.sim.trees)<-Pair_metadata$Pair

#Classify 'type' of selection by looking at pattern of coalescences in D and R trees

all_pair_stats<-Map(pair_trees=all.sim.trees,pair=Pair_metadata$Pair,f=function(pair_trees,pair) {
  cat(pair,sep="\n")
  names(pair_trees)<-c("NULL_model_trees",
                         "Type1_model_trees",
                         "Type2_model_trees")
  pair_stats<-Map(model_trees=pair_trees,model=names(pair_trees),function(model_trees,model) {
    cat(model,sep="\n")
    model_tree_stats<-Map(tree=model_trees,idx=1:length(model_trees),function(tree,idx) {
      cat(idx,sep="\n")
      age_of_donor_at_sampling<-Pair_metadata%>%filter(Pair==pair)%>%pull(Age)
      age_of_donor_at_HSCT<-Pair_metadata%>%filter(Pair==pair)%>%pull(Age_at_transplant)
      
      tree.time<-make.ultrametric.tree(tree)
      tree.time$edge.length[which(tree$edge[,2]==grep("donor|recip",tree.time$tip.label,value=F,invert=T))]<-0
      tree.time<-drop.tip(tree.time,grep("donor|recip",tree.time$tip.label,value=T,invert=T))
      tree.time$edge.length<-tree.time$edge.length*age_of_donor_at_sampling
      D_tree.time<-keep.tip(tree.time,tip=grep("donor",tree.time$tip.label))
      R_tree.time<-keep.tip(tree.time,tip=grep("recip",tree.time$tip.label))
      
      #Find the expanded clades at 100 muts of molecular time
      expanded_clades=get_expanded_clade_nodes(tree=tree.time,height_cut_off = 5,min_samples = 3)
      if(nrow(expanded_clades)==0) {stop(return(NULL))}
      
      #Get the donor & recipient samples within each expanded clade
      clade_samples=lapply(expanded_clades$nodes,function(node) getTips(tree.time,node=node))
      D_clade_samples=lapply(clade_samples,function(samples) {return(grep("donor",samples,value=T))})
      R_clade_samples=lapply(clade_samples,function(samples) {return(grep("recip",samples,value=T))})
      
      #Update the expanded clades dataframe with extra stats around the numbers of donor/ recipient samples
      expanded_clades$Pair=pair
      expanded_clades$D_total=length(D_tree.time$tip.label)
      expanded_clades$R_total=length(R_tree.time$tip.label)
      expanded_clades$D_nodes=sapply(D_clade_samples,function(samples) if(length(samples)==0) {return(NA)} else {find_latest_acquisition_node(tree=D_tree.time,pos_samples = samples)})
      expanded_clades$R_nodes=sapply(R_clade_samples,function(samples) if(length(samples)==0) {return(NA)} else {find_latest_acquisition_node(tree=R_tree.time,pos_samples = samples)})
      expanded_clades$n_D_samples=sapply(D_clade_samples,length)
      expanded_clades$n_R_samples=sapply(R_clade_samples,length)
      
      #Define the "peri-HSCT" time points - this is the window +/-5 years around the known time of transplant to reflect error in the timing due to poisson error
      peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
      
      #Select out the nodes (in donor & recipient) that fall into this time zone
      D_peri_HSCT_nodes<-D_tree.time$edge[nodeHeights(D_tree.time)[,2]>peri_HSCT_time_points[1]&nodeHeights(D_tree.time)[,2]<peri_HSCT_time_points[2],2]
      R_peri_HSCT_nodes<-R_tree.time$edge[nodeHeights(R_tree.time)[,2]>peri_HSCT_time_points[1]&nodeHeights(R_tree.time)[,2]<peri_HSCT_time_points[2],2]
      
      #Define pre-HSCT time points
      pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
      
      #Select out the nodes (in donor & recipient) that fall into this time zone
      D_pre_HSCT_nodes<-D_tree.time$edge[nodeHeights(D_tree.time)[,2]>pre_HSCT_time_points[1]&nodeHeights(D_tree.time)[,2]<pre_HSCT_time_points[2],2]
      R_pre_HSCT_nodes<-R_tree.time$edge[nodeHeights(R_tree.time)[,2]>pre_HSCT_time_points[1]&nodeHeights(R_tree.time)[,2]<pre_HSCT_time_points[2],2]
      
      #Now get list of the post development nodes (post 5 years) in donor & recipient for bootstrapping
      D_all_PDNs <-D_tree.time$edge[nodeHeights(D_tree.time)[,2]>5 & !D_tree.time$edge[,2]%in%1:length(D_tree.time$tip.label),2]
      R_all_PDNs <-R_tree.time$edge[nodeHeights(R_tree.time)[,2]>5 & !R_tree.time$edge[,2]%in%1:length(R_tree.time$tip.label),2]
      
      #Now get a list of the pre-HSCT lineages
      R_node_heights=nodeHeights(R_tree.time)
      D_node_heights=nodeHeights(D_tree.time)
      
      R_pre_lineages=R_tree.time$edge[,2][which(R_node_heights[,1]<peri_HSCT_time_points[1] & R_node_heights[,2]>peri_HSCT_time_points[1])]
      D_pre_lineages=D_tree.time$edge[,2][which(D_node_heights[,1]<peri_HSCT_time_points[1] & D_node_heights[,2]>peri_HSCT_time_points[1])]
      
      #Now go through the clades and define these statistics
      expanded_clades_full<-bind_cols(expanded_clades,
                                      lapply(1:nrow(expanded_clades),function(i) {
                                        D_clade_nodes<-c(expanded_clades$D_nodes[i],get_all_node_children(expanded_clades$D_nodes[i],tree=D_tree.time))
                                        R_clade_nodes<-c(expanded_clades$R_nodes[i],get_all_node_children(expanded_clades$R_nodes[i],tree=R_tree.time))

                                        #Define the node sets that are (1) in the tested clade & (2) at the relevant time
                                        D_clade_peri_nodes<-intersect(D_peri_HSCT_nodes,D_clade_nodes)
                                        R_clade_peri_nodes<-intersect(R_peri_HSCT_nodes,R_clade_nodes)
                                        
                                        D_clade_pre_nodes<-intersect(D_pre_HSCT_nodes,D_clade_nodes)
                                        R_clade_pre_nodes<-intersect(R_pre_HSCT_nodes,R_clade_nodes)
                                        
                                        #Define the actual stats based on these
                                        n_peri.D=sum(D_all_PDNs%in%D_clade_peri_nodes)
                                        n_peri.R=sum(R_all_PDNs%in%R_clade_peri_nodes)
                                        
                                        n_pre.D=sum(D_all_PDNs%in%D_clade_pre_nodes)
                                        n_pre.R=sum(R_all_PDNs%in%R_clade_pre_nodes)
                                        
                                        #Pre-HSCT lineages in the clade in question
                                        n_pre_lineages_in_clade.D=sum(D_pre_lineages%in%D_clade_nodes)
                                        n_pre_lineages_in_clade.R=sum(R_pre_lineages%in%R_clade_nodes)
                                        
                                        #Proportion of peri-HSCT coalescences in the clade
                                        R_prop_pre_lineages_in_clade<-n_pre_lineages_in_clade.R/length(R_pre_lineages)
                                        R_prop_of_peri_nodes_in_clade<-length(R_clade_peri_nodes)/length(R_peri_HSCT_nodes)
                                        
                                        #Probability of all peri-HSCT coalescences being from within this clade by chance
                                        p_clade_skewing<-binom.test(x = n_peri.R,n = length(R_peri_HSCT_nodes),p = R_prop_pre_lineages_in_clade,alternative="g")$p.value
                                        
                                        selection_stats_df<-data.frame(n_peri.D=n_peri.D,
                                                                       n_peri.R=n_peri.R,
                                                                       n_pre.D=n_pre.D,
                                                                       n_pre.R=n_pre.R,
                                                                       Type1_selection_stat=((n_pre.R+1)/(expanded_clades$R_total[i]+1))/((n_pre.D+1)/(1+expanded_clades$D_total[i])),
                                                                       Type2_selection_stat=((n_peri.R+1)/(expanded_clades$R_total[i]+1))/((n_peri.D+1)/(1+expanded_clades$D_total[i])),
                                                                       p_skew=p_clade_skewing)
                                        return(selection_stats_df)
                                        
                                      })%>%dplyr::bind_rows())
    })
    names(model_tree_stats)<-paste0("Posterior_tree_",1:length(model_tree_stats))
  })
  return(pair_stats)
})

model_all<-Map(df=model_tree_stats,sim=names(model_tree_stats),f=function(df,sim) {
  df<-df%>%mutate(sim=sim)
  })%>%
  dplyr::bind_rows()

include_ids<-model_all%>%
  mutate(ratio=((1+n_R_samples)/(1+R_total))/((1+n_D_samples)/(1+D_total)))%>%
  filter((n_D_samples+n_R_samples)>9 & ratio>1.25)%>%
  mutate(id=paste(sim,nodes,sep="_"))%>%pull(id)

plot.selection.type<-model_all%>%
  mutate(id=paste(sim,nodes,sep="_"))%>%
  filter(id%in%include_ids)%>%
  ggplot(aes(x=log(Type1_selection_stat),y=log(Type2_selection_stat)))+
  geom_abline(slope = 1)+
  scale_x_continuous(limits=c(-1,5.1))+
  scale_y_continuous(limits=c(-1,4.6))+
  geom_jitter(height=0.2,width=0.2)+
  theme_classic()+
  my_theme+
  labs(x="Log of the Type 1 selection statistic",y="Log of the Type 2 selection statistic")
  # geom_polygon(data=confidence_interval_df%>%filter(id%in%include_ids),aes(x=x,y=y,fill=Pair_new,group=id),alpha=0.2,inherit.aes = F)+

  ggrepel::geom_label_repel(size=1.5,show.legend = F)+
  theme_classic()+
  my_theme+
  #geom_vline(xintercept=0)+
  scale_color_manual(values=Pair_cols)+
  scale_fill_manual(values=Pair_cols)+
  labs(col="Pair",fill="Pair",x="Log of the Type 1 selection statistic",y="Log of the Type 2 selection statistic")
