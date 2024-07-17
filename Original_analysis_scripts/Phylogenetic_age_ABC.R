library(optparse)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(phytools)
library(rsimpop)


root_dir<-ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/Clonal_dynamics_of_HSCT","/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT")
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch119/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
source(paste0(root_dir,"/data/HSCT_functions.R"))
setwd(tree_mut_dir); source("treemut.R");setwd(root_dir)
plots_dir=paste0(root_dir,"/plots/")

#plots_dir="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Age/plots/"

###-----------------SPECIFIC FUNCTIONS FOR THIS SCRIPT-----------------
get_ltt = function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  ltt_tree = sapply(time_points, function(x) {
    sum(nodeheights[,1] < x & !nodeheights[,2] < x)
  })
  return(ltt_tree)
}

get_coalescences = function(ltt) {
  coals=sapply(2:length(ltt), function(i) {return(ltt[i]-ltt[i-1])})
  return(coals)
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

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

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
tree_folder=paste0(root_dir,"/data/trees_no_dups")
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


colony_numbers_df<-Map(tree=all.trees,pair=names(all.trees),function(tree,pair) {
  DR_ids<-get_DR_ids(tree)
  
  no_of_WGS_recip=sum(grepl(DR_ids['recip_ID'],tree$tip.label))
  no_of_WGS_donor=sum(grepl(DR_ids['donor_ID'],tree$tip.label))
  
  return(data.frame(PairID=pair,Donor=no_of_WGS_donor,Recip=no_of_WGS_recip))
})%>%dplyr::bind_rows()%>%
  gather(-PairID,key="D_or_R",value="nWGS")%>%
  mutate(ID=paste(PairID,D_or_R,sep = "_"))

setwd(root_dir)


#combined_output_file="/lustre/scratch126/realdata/mdt1/team154/ms56/Zur_HSCT/ABC_models/ABC_Age/combined_output.Rds"
res_by_pair_file=paste0(root_dir,"/data/ABC_simulation_results/Phylogenetic_age_ABC_simulation_results/res_by_pair.Rds")


if(file.exists(res_by_pair_file)) {
  res_by_pair=readRDS(res_by_pair_file)
} else {
  
  ##Importing full original ABC output
  if(file.exists(combined_output_file)) {
    all_sim_files=readRDS(combined_output_file)
  } else {
    #Import all the age simulation output files - this takes a while
    sim_files=list.files(pattern=".Rds")
    length(sim_files)
    all_sim_files=lapply(sim_files,function(file) {
      {print(which(sim_files==file))}
      res<-readRDS(file)
      return(res)
    })
    saveRDS(all_sim_files,file=combined_output_file)
  }
  
  ##Get sumstats for each individual
  res_by_pair<-lapply(1:nrow(colony_numbers_df),function(i) {
    cat(i,sep="\n")
    pair_ss<-lapply(all_sim_files,function(res) {
      pair_res<-res[[i]]
      pair_res$sumstats
    })%>%bind_rows()
    
    pair_params<-lapply(all_sim_files,function(res) {
      res[[i]]$params
    })%>%bind_rows()
    
    return(list(sumstats=pair_ss,params=pair_params))
  })
  names(res_by_pair)<-colony_numbers_df$ID
  
  ##Get trees for each individual
  phylos_by_pair<-lapply(1:nrow(colony_numbers_df),function(i) {
    cat(i,sep="\n")
    pair_trees<-lapply(all_sim_files,function(res) {
      res[[i]]$tree
    })
    class(pair_trees)<-"multiPhylo"
    return(pair_trees)
  })
}



# setwd("pair_output")
# temp=Map(phylos=phylos_by_pair,ID=names(res_by_pair),function(phylos,ID) {
#   write.tree(phylos,file=paste0("trees_",ID,".trees"))
# })
# 
# temp=Map(pair_res=res_by_pair,ID=names(res_by_pair),function(pair_res,ID){
#   write.table(pair_res$sumstats,file = paste0("sumstats_",ID,".tsv"),quote = F,row.names = F,sep = "\t")
#   write.table(pair_res$params,file = paste0("params_",ID,".tsv"),quote = F,row.names = F,sep = "\t")
# })

data.sumstats<-lapply(1:nrow(colony_numbers_df),function(i) {
  PairID=colony_numbers_df$PairID[i]
  D_or_R=colony_numbers_df$D_or_R[i]
  nWGS=colony_numbers_df$nWGS[i]
  
  cat(paste(PairID,D_or_R,sep="_"),sep="\n")
  
  tree<-all.trees[[PairID]]
  tree.ultra<-make.ultrametric.tree(tree=tree)
  PD_IDs<-get_DR_ids(tree = tree)
  
  if(D_or_R=="Donor"){
    sub_tree<-drop.tip(tree,tip=grep(PD_IDs['recip_ID'],tree$tip.label,value=T))
    sub_tree.ultra<-drop.tip(tree.ultra,tip=grep(PD_IDs['recip_ID'],tree$tip.label,value=T))
    
  } else if(D_or_R=="Recip") {
    sub_tree<-drop.tip(tree,tip=grep(PD_IDs['donor_ID'],tree$tip.label,value=T))
    sub_tree.ultra<-drop.tip(tree.ultra,tip=grep(PD_IDs['donor_ID'],tree$tip.label,value=T))
  }
  
  
  ##Get sumstats for the data
  ##EXTRACT SUMMARY STATISTICS
  #(1) 3 LARGEST CLADES
  largest_clades<-get_expanded_clade_nodes(sub_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]
  
  #(2) Number of singletons
  n_singletons<-get_expanded_clade_nodes(sub_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    dplyr::filter(n_samples==1)%>%nrow(.)
  
  #(3) Peri-HSCT LTT/ coalescences
  ltt<-get_ltt(tree=sub_tree.ultra,time_points=c(0.2,0.4,0.6))
  coals<-get_coalescences(ltt)
  
  #(4) Bulk totals for expanded populations to look at clonal shifts
  expanded_df=get_expanded_clade_nodes(sub_tree,min_clonal_fraction=0.001,height_cut_off=50,min_samples=3)
  
  if(nrow(expanded_df)>0){
    expanded_df<-expanded_df%>%
      mutate(count=sapply(nodes,function(node) {length(getTips(sub_tree,node))}))%>%
      mutate(prop=count/nWGS)
    
    prop_from_expansions=sum(expanded_df$prop)
    
  } else {
    prop_from_expansions=0
  }
  
  ##--COMBINE SUMSTATS INTO SINGLE VECTOR--
  ss_names=c(paste("largest_clades",1:3,sep="_"),
             "n_singletons",
             paste("coals",1:2,sep="_"),
             "prop_from_expansions")
  
  ss_comb=c(largest_clades,
            n_singletons,
            coals,
            prop_from_expansions)
  
  names(ss_comb)<-ss_names
  return(ss_comb)
})

# data.sumstats<-bind_rows(data.sumstats)%>%as.matrix()
# rownames(data.sumstats)<-paste(colony_numbers_df$PairID,colony_numbers_df$D_or_R,sep="_")

#Perform the ABC over all donor/ recipient trees
library(abc)
abc.res.all<-list()
method="neuralnet"
for(i in 1:length(res_by_pair)) {
  res<-res_by_pair[[i]]
  target<-data.sumstats[[i]]
  ID<-names(res_by_pair)[i]
  
  cat(ID,sep="\n")
  params<-res$params
  sumstats<-res$sumstats
  
  if(grepl("Pair41",ID)|grepl("Pair40_Donor",ID)){
    #Need to take a broader tolerance for Pair41/ Pair40 donor as otherwise there is inadequate variance of the summary statistics to perform the regression.
    abc.res.all[[i]]<-abc(target=target,param=params[,1],sumstat=sumstats,method = method,tol=0.1,transf = "none")
  } else {
    abc.res.all[[i]]<-abc(target=target,param=params[,1],sumstat=sumstats,method = method,tol=0.05,transf = "none")
  }
  
  #If want to do rejection method (but then also need to change to 'unadj.values' in next code snippet)
  # abc.res.all[[i]]<-abc(target=target,param=params,sumstat=sumstats,method = "rejection",tol=0.05,transf = "none")
}

#Extract the age posteriors
age.post<-Map(res=abc.res.all,ID=names(res_by_pair),function(res,ID) {
  post<-res[[ifelse(method=="rejection","unadj.values","adj.values")]]
  return(data.frame(ID=ID,tree_age=post))
})
names(age.post)<-names(res_by_pair)

#Create output in format for linear mixed effects model
age.post.for.lme<-age.post%>%
  bind_rows()%>%
  separate(ID,into=c("PairID","D_or_R"),sep="_")%>%
  left_join(Pair_metadata%>%dplyr::select(Pair,Age),by=c("PairID"="Pair"))

##Perform the LME
library(lme4)

lmer.age_only <- lmer(tree_age ~ Age + (1|PairID), data = age.post.for.lme)
summary(lmer.age_only)

lmer.age_and_DR_status <- lmer(tree_age ~ Age + D_or_R +(1|PairID), data = age.post.for.lme)
lme.summary<-summary(lmer.age_and_DR_status)
confint(lmer.age_and_DR_status)

lme.plot<-ggplot(age.post.for.lme,aes(x=Age,y=tree_age,col=factor(D_or_R)))+
  #geom_abline(slope=1,intercept=0,linetype=2,col="black")+
  geom_jitter(size=0.05,alpha=0.01,height=0,width = 1)+
  #geom_smooth(method="lm")+
  geom_abline(slope = lme.summary$coefficients['Age','Estimate'],intercept = lme.summary$coefficients['(Intercept)','Estimate'],col=DorR_cols[1])+
  geom_abline(slope = lme.summary$coefficients['Age','Estimate'],intercept = lme.summary$coefficients['D_or_RRecip','Estimate']+lme.summary$coefficients['(Intercept)','Estimate'],col=DorR_cols[2])+
  scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(20,100),breaks=seq(0,100,10))+
  scale_y_continuous(limits=c(0,110),breaks=seq(0,110,10))+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1)))+
  theme_classic()+
  my_theme+
  labs(col="",x="Actual age",y="Tree age posterior")
ggsave(filename=paste0(plots_dir,"lme.plot.pdf"),plot = lme.plot,width=3,height = 2.5)

##Other plots of differences
age.post.for.lme%>%
  ggplot(aes(x=tree_age))+
  geom_histogram()+
  facet_grid(D_or_R~PairID)+
  theme_classic()

tree_age_vs_DorR<-age.post%>%
  bind_rows()%>%
  group_by(ID)%>%
  summarise(median=median(tree_age),lowerCI=quantile(tree_age,0.025),upperCI=quantile(tree_age,0.975))%>%
  separate(ID,into=c("PairID","D_or_R"),sep="_",remove=F)%>%
  mutate(D_or_R=substr(D_or_R,1,1))%>%
  mutate(Pair_new=factor(new_pair_names[PairID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=D_or_R,y=median,ymin=lowerCI,ymax=upperCI,col=D_or_R))+
  geom_point()+
  geom_errorbar(width=0.2)+
  facet_grid(~Pair_new,scales="free",space="free")+
  scale_color_manual(values=DorR_cols)+
  theme_classic()+
  my_theme+
  theme(legend.position="none")+
  labs(x="Donor (D) or Recipient (R)",y="Predicted age from phylogeny")

tree_age_vs_age_plot<-age.post%>%
  bind_rows()%>%
  group_by(ID)%>%
  summarise(median=median(tree_age),lowerCI=quantile(tree_age,0.025),upperCI=quantile(tree_age,0.975))%>%
  separate(ID,into=c("PairID","D_or_R"),sep="_",remove=F)%>%
  mutate(Pair_new=factor(new_pair_names[PairID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  left_join(Pair_metadata,by="Pair_new")%>%
  ggplot(aes(x=Age,y=median,ymin=lowerCI,ymax=upperCI,col=Pair_new))+
  geom_abline(slope=1,intercept=0,linetype=2)+
  geom_point(size=1,alpha=0.5)+
  geom_errorbar(alpha=0.5,width=2)+
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20))+
  scale_y_continuous(limits=c(0,110),breaks=seq(0,100,20))+
  scale_color_manual(values=Pair_cols)+
  facet_grid(~D_or_R)+
  theme_classic()+
  my_theme+
  theme(legend.key.height = unit(3,"mm"))+
  labs(x="Actual age",y="Predicted age from phylogeny",col="")
ggsave(filename=paste0(plots_dir,"tree_age_vs_age_plot.pdf"),plot = tree_age_vs_age_plot,width=5,height = 2.5)

age.post%>%
  bind_rows()%>%
  group_by(ID)%>%
  summarise(age=median(tree_age))%>%
  separate(ID,into=c("PairID","D_or_R"),sep="_")%>%
  pivot_wider(names_from="D_or_R",values_from="age")%>%
  mutate(change=Recip-Donor)


#Summarise the difference by subtracting the posteriors
change_in_recip_tree_age<-lapply(Pair_metadata$Pair,function(pair){
  recip.post<-age.post[[paste(pair,"Recip",sep="_")]]
  donor.post<-age.post[[paste(pair,"Donor",sep="_")]]
  
  return(data.frame(PairID=pair,R_post=(recip.post%>%dplyr::pull(age)),D_post=(donor.post%>%dplyr::pull(age)))%>%
    mutate(Diff=R_post-D_post))
})%>%dplyr::bind_rows()%>%
  group_by(PairID)%>%
  summarise(median=median(Diff),lowerCI=quantile(Diff,0.25),upperCI=quantile(Diff,0.75))%>%
  mutate(PairID=factor(PairID,levels=Pair_metadata$Pair[order(Pair_metadata$Age)]))%>%
  ggplot(aes(x=PairID,y=median,ymin=lowerCI,ymax=upperCI,col=PairID))+
  geom_point(size=0.6)+
  geom_errorbar(width=0.2,size=0.4)+
  geom_hline(yintercept=0,linetype=2)+
  scale_color_brewer(palette="Paired")+
  theme_classic()+
  my_theme+
  theme(legend.key.height = unit(3,"mm"),axis.text.x=element_text(angle=90))+
  labs(x="Pair ID",y="Increase in recipient 'tree age'",col="Pair ID")
ggsave(filename=paste0(plots_dir,"change_in_recip_tree_age.pdf"),plot = change_in_recip_tree_age,width=2.5,height = 2)



#Calculate the mean of medians (i.e. the mean most likely shift in "tree age" from ABC)
lapply(Pair_metadata$Pair,function(pair){
  recip.post<-age.post[[paste(pair,"Recip",sep="_")]]
  donor.post<-age.post[[paste(pair,"Donor",sep="_")]]
  
  return(data.frame(PairID=pair,R_post=(recip.post%>%dplyr::pull(age)),D_post=(donor.post%>%dplyr::pull(age)))%>%
           mutate(Diff=R_post-D_post))
})%>%dplyr::bind_rows()%>%
  group_by(PairID)%>%
  summarise(median=median(Diff),lowerCI=quantile(Diff,0.25),upperCI=quantile(Diff,0.75))%>%
  summarise(mean(median))

