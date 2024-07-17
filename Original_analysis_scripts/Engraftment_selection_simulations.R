#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("abc","ggridges","ggplot2","tidyr","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lme4","lmerTest")
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
if(!require("rsimpop", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("NickWilliamsSanger/rsimpop")
  library("rsimpop",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Set the ggplot2 theme for plotting
#----------------------------------

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#----------------------------------
# Set the root directory and read in the necessary files
#----------------------------------

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
source(paste0(root_dir,"/data/HSCT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in the trees list for the ABC
ABC.trees<-readRDS(paste0(root_dir,"/data/trees_for_ABC.Rds"))

#----------------------------------
# SPECIFIC FUNCTIONS FOR THIS SCRIPT
#----------------------------------

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

#----------------------------------
# GET SUMMARY STATISTICS FROM THE DATA
#----------------------------------

##EXTRACT SUMMARY STATISTICS
#(1) 3 LARGEST CLADES
all.ss<-Map(tree=ABC.trees,pair=names(ABC.trees),function(tree,pair) {
  cat(pair,sep="\n")
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

#----------------------------------
# IMPORT PARAMETERS AND SUMMARY STATISTICS FOR THE SIMULATIONS & PERFORM ABC
#----------------------------------

setwd(paste0(root_dir,"/data/ABC_simulation_results/Engrafting_cell_number_ABC_Engraftment_selection_simulation_results"))

#Can run with different sets of summary stats
all_tree_sumstats=1:15 #All the phylo statistics (not including targeted sequencing results)
tree_sumstats_no_Rpre_coals=c(1,2,3,4,6,8,9,10,11,12,13,14)

#Now extract the parameters/ summary statistics and run the ABC
param_posteriors<-lapply(names(ABC.trees),function(pair){
  cat(pair,sep="\n")
  
  sim.files=list.files(path=pair,pattern=".Rds",full.names = T)
  params_file=paste0("all_params_",pair,".Rds")
  sumstats_file=paste0("all_sumstats_",pair,".Rds")
  
  if(file.exists(params_file) & file.exists(sumstats_file)){
    params<-readRDS(params_file)
    sumstats<-readRDS(sumstats_file)
  } else {
    sim.output<-lapply(sim.files,function(file) {if(which(sim.files==file)%%100 == 0) {print(which(sim.files==file))};readRDS(file)})
    
    params<-lapply(sim.output,function(res) unlist(res$params))%>%dplyr::bind_rows()
    sumstats<-lapply(sim.output,function(res) unlist(res$sumstats))%>%dplyr::bind_rows()
    
    saveRDS(params,file=params_file)
    saveRDS(sumstats,file=sumstats_file)
  }
  
  
  #For assessing engrafting HSC number ignore the R pre-HSCT coalescences (which distorts the posterior as ABC is unable to fit this)
  which_sumstats<-c(1:8,10:11,13:15)
  #which_sumstats<-tree_sumstats_no_Rpre_coals
  
  #Run the ABC - log transform the parameters for the regression step. Use tolerance of 5% (= ~ top 1000 for 20,000 simulations)
  abc.res<-abc(target = all.ss[[pair]][which_sumstats],param = params[,6:7],sumstat=sumstats[,which_sumstats],tol = 0.05,transf = c("log","log"),method="neuralnet")
  
  closest_fit<-which.min(abc.res$dist)
  as.matrix(sumstats[closest_fit,])
  
  data.ss.df<-data.frame(sumstat=names(all.ss[[pair]][which_sumstats]),value=all.ss[[pair]][which_sumstats])
  
  return(as.data.frame(abc.res$adj.values)%>%gather(key="parameter",value="value")%>%mutate(PairID=pair,.before=1))
})%>%bind_rows()

#----------------------------------
# VISUALIZE RESULTS - INCLUDING CORRELATION WITH OTHER HSCT PARAMETERS
#----------------------------------

HSC_pop_posteriors_plot<-param_posteriors%>%
  mutate(PairID=factor(new_pair_names[PairID],levels=Pair_metadata$Pair_new[order(Pair_metadata$Age)]))%>%
  filter(parameter=="HSC_pop_size")%>%
  ggplot(aes(x=value,y=PairID,fill=PairID))+
  ggridges::geom_density_ridges(scale=5,size=0.3)+
  scale_fill_manual(values=Pair_cols)+
  scale_x_continuous(breaks=seq(5e4,4e5,5e4),limits=c(4e4,3.5e5),labels=scales::label_comma())+
  #scale_x_log10(breaks=c(300,1000,3000,10000,30000,100000),labels=scales::label_comma())+
  theme_classic()+
  my_theme+
  labs(x="HSC population size",y="Density")+
  theme(legend.key.size = unit(3,"mm"),axis.text.x = element_text(angle=90))

HSCT_bottleneck_posteriors_plot<-param_posteriors%>%
  mutate(PairID=factor(new_pair_names[PairID],levels=Pair_metadata$Pair_new[order(Pair_metadata$Age)]))%>%
  filter(parameter=="HSCT_bottleneck")%>%
  ggplot(aes(x=value,y=PairID,fill=PairID))+
  ggridges::geom_density_ridges(scale=5,size=0.3)+
  scale_fill_manual(values=Pair_cols)+
  scale_x_log10(breaks=c(300,1000,3000,10000,30000,100000),labels=scales::label_comma())+
  theme_classic()+
  my_theme+
  labs(x="# of engrafting HSCs",y="Density",fill="Pair")+
  theme(legend.key.size = unit(3,"mm"),axis.text.x = element_text(angle=90))

ggsave(filename = paste0(plots_dir,"HSC_pop_posteriors_plot_Esim.pdf"),HSC_pop_posteriors_plot,width=2.5,height=2)
ggsave(filename = paste0(plots_dir,"HSCT_bottleneck_posteriors_plot_Esim.pdf"),HSCT_bottleneck_posteriors_plot,width=2.5,height=2)

HSCT_correlates<-param_posteriors%>%
  mutate(PairID=factor(new_pair_names[PairID],levels=Pair_metadata$Pair_new[order(Pair_metadata$Age)]))%>%
  filter(parameter=="HSCT_bottleneck")%>%
  group_by(PairID)%>%
  summarise(median=median(value),lowerCI=quantile(value,0.025),upperCI=quantile(value,0.975))%>%
  left_join(Pair_metadata,by=c("PairID"="Pair_new"))%>%
  mutate(stem_cell_source=ifelse(stem_cell_source=="BM","Bone marrow","GCSF-mobilized PB"),
         time_since_HSCT=Age-Age_at_transplant)

ECN_by_age_at_sampling<-HSCT_correlates%>%
  ggplot(aes(x=Age,y=median,ymin=lowerCI,ymax=upperCI,col=PairID,shape=stem_cell_source))+
  geom_point()+
  geom_errorbar(width=1)+
  scale_color_manual(values=Pair_cols)+
  scale_y_log10(labels=scales::label_comma())+
  labs(x="Age at time of sampling",y="# of engrafting HSCs",col="Pair",shape="HSC source")+
  theme_classic()+
  my_theme+
  theme(legend.key.height = unit(3,"mm"))

ECN_by_age_at_HSCT<-HSCT_correlates%>%
  ggplot(aes(x=Age_at_transplant,y=median,ymin=lowerCI,ymax=upperCI,col=PairID,shape=stem_cell_source))+
  geom_point()+
  geom_errorbar(width=1)+
  scale_color_manual(values=Pair_cols)+
  scale_y_log10(labels=scales::label_comma())+
  labs(x="Age at HSCT",y="# of engrafting HSCs",col="Pair",shape="HSC source")+
  theme_classic()+
  my_theme+
  theme(legend.key.height = unit(3,"mm"))

ECN_by_time_since_HSCT<-HSCT_correlates%>%
  ggplot(aes(x=time_since_HSCT,y=median,ymin=lowerCI,ymax=upperCI,col=PairID,shape=stem_cell_source))+
  geom_point()+
  geom_errorbar(width=1)+
  scale_color_manual(values=Pair_cols)+
  scale_y_log10(labels=scales::label_comma())+
  labs(x="Years since HSCT",y="# of engrafting HSCs",col="Pair",shape="HSC source")+
  theme_classic()+
  my_theme+
  theme(legend.key.height = unit(3,"mm"))


HSCT_regression_df<-param_posteriors%>%
  mutate(PairID=factor(new_pair_names[PairID],levels=Pair_metadata$Pair_new[order(Pair_metadata$Age)]))%>%
  filter(parameter=="HSCT_bottleneck")%>%
  left_join(Pair_metadata,by=c("PairID"="Pair_new"))%>%
  mutate(time_since_HSCT=Age-Age_at_transplant)

lme.age<-lmerTest::lmer(log(value)~Age+(1|PairID),data=HSCT_regression_df)
lme.age_at_HSCT<-lmerTest::lmer(log(value)~Age_at_transplant+(1|PairID),data=HSCT_regression_df)
lme.CD34<-lmerTest::lmer(log(value)~CD34_dose+(1|PairID),data=HSCT_regression_df)
lme.HSC_source<-lmerTest::lmer(log(value)~stem_cell_source+(1|PairID),data=HSCT_regression_df)
lme.conditioning<-lmerTest::lmer(log(value)~conditioning+(1|PairID),data=HSCT_regression_df)
lme.time_since_HSCT<-lmerTest::lmer(log(value)~time_since_HSCT+(1|PairID),data=HSCT_regression_df)

summary(lme.age);confint(lme.age)
summary(lme.age_at_HSCT);confint(lme.age_at_HSCT)
summary(lme.CD34);confint(lme.CD34)
summary(lme.HSC_source);confint(lme.HSC_source)
summary(lme.conditioning);confint(lme.conditioning)
summary(lme.time_since_HSCT);confint(lme.time_since_HSCT)


ECN_by_CD34<-HSCT_correlates%>%
  ggplot(aes(x=CD34_dose,y=median,ymin=lowerCI,ymax=upperCI,col=PairID,shape=stem_cell_source))+
  geom_point()+
  geom_errorbar(width=0.3)+
  scale_color_manual(values=Pair_cols)+
  scale_y_log10(labels=scales::label_comma())+
  labs(x="CD34 dose (million/ kg)",y="# of engrafting HSCs",col="Pair",shape="HSC source")+
  theme_classic()+
  my_theme+
  theme(legend.key.height = unit(3,"mm"))

ECN_by_MNC<-HSCT_correlates%>%
  ggplot(aes(x=MNC_dose,y=median,ymin=lowerCI,ymax=upperCI,col=PairID,shape=stem_cell_source))+
  geom_point()+
  geom_errorbar(width=0.3)+
  scale_color_manual(values=Pair_cols)+
  scale_y_log10(labels=scales::label_comma())+
  labs(col="Pair",x="MNC dose (million/ kg)",y="# of engrafting HSCs",shape="HSC source")+
  theme_classic()+
  my_theme+
  theme(legend.key.height = unit(3,"mm"))

ggsave(filename = paste0(plots_dir,"ECN_by_age_at_sampling_ESim.pdf"),ECN_by_age_at_sampling,width=2.5,height=2)
ggsave(filename = paste0(plots_dir,"ECN_by_age_at_HSCT_ESim.pdf"),ECN_by_age_at_HSCT,width=2.5,height=2)
ggsave(filename = paste0(plots_dir,"ECN_by_CD34_ESim.pdf"),ECN_by_CD34,width=2.5,height=2)
ggsave(filename = paste0(plots_dir,"ECN_by_MNC_ESim.pdf"),ECN_by_MNC,width=2.5,height=2)
ggsave(filename = paste0(plots_dir,"ECN_by_time_since_HSCT_ESim.pdf"),ECN_by_time_since_HSCT,width=3,height=2)

#----------------------------------
# NOW DO POSTERIOR CHECKS OF SUMMARY STATISTICS
#----------------------------------

#Function to condense the pre/ peri-HSCT coalesnces summary stats into D/R ratios
condense_sumstats=function(sumstats) {
  condense_columns=c("R_coals_pre","D_coals_pre","R_coals_peri","D_coals_peri")
  if(is.null(dim(sumstats))) {
    if(!all(condense_columns%in%names(sumstats))) {
      stop(print("The columns to condense are not in the target sumstats"))
    } else {
      R_D_coals_pre_ratio=log((1+sumstats["R_coals_pre"])/(1+sumstats["D_coals_pre"]))
      R_D_coals_peri_ratio=log((1+sumstats["R_coals_peri"])/(1+sumstats["D_coals_peri"]))
      names(R_D_coals_pre_ratio)<-"R_D_coals_pre_ratio"
      names(R_D_coals_peri_ratio)<-"R_D_coals_peri_ratio"
      
      new_sumstats<-c(sumstats[-which(names(sumstats)%in%condense_columns)],R_D_coals_pre_ratio,R_D_coals_peri_ratio)
      return(new_sumstats)
    }
  } else {
    if(!all(condense_columns%in%colnames(sumstats))) {
      stop(print("The columns to condense are not in the simulation sumstats"))
    } else {
      R_D_coals_pre_ratio=data.frame(log((1+sumstats[,"R_coals_pre"])/(1+sumstats[,"D_coals_pre"])))
      R_D_coals_peri_ratio=data.frame(log((1+sumstats[,"R_coals_peri"])/(1+sumstats[,"D_coals_peri"])))
      colnames(R_D_coals_pre_ratio)<-"R_D_coals_pre_ratio"
      colnames(R_D_coals_peri_ratio)<-"R_D_coals_peri_ratio"
      
      new_sumstats<-cbind(sumstats%>%dplyr::select(-all_of(condense_columns)),R_D_coals_pre_ratio,R_D_coals_peri_ratio)
      return(new_sumstats)
    }
  }
}

#Get the summary stats for posterior 
posterior_checks<-lapply(names(ABC.trees),function(pair){
  cat(pair,sep="\n\n")
  
  sim.files=list.files(path=pair,pattern=".Rds",full.names = T)
  params_file=paste0("all_params_",pair,".Rds")
  sumstats_file=paste0("all_sumstats_",pair,".Rds")
  trees_file=paste0("all_trees_",pair,".tree")
  
  if(file.exists(params_file) & file.exists(sumstats_file)){
    cat("Reading in existing parameter and summary statistics files",sep="\n")
    params<-readRDS(params_file)
    sumstats<-readRDS(sumstats_file)
  } else {
    cat("Extracting parameter and summary statistics from simulation results",sep="\n")
    sim.output<-lapply(sim.files,function(file) {if(which(sim.files==file)%%100 == 0) {print(which(sim.files==file))};readRDS(file)})
    
    params<-lapply(sim.output,function(res) unlist(res$params))%>%dplyr::bind_rows()
    sumstats<-lapply(sim.output,function(res) unlist(res$sumstats))%>%dplyr::bind_rows()
    multi.phylo=lapply(sim.output,function(res) res$DR_tree)
    class(multi.phylo)<-"multiPhylo"
    
    saveRDS(params,file=params_file)
    saveRDS(sumstats,file=sumstats_file)
    write.tree(multi.phylo,file=trees_file)
  }
  
  
  #Choose the summary statistic set to use
  which_sumstats<-all_tree_sumstats #These are the phylo summary stats (not the targeted sequencing)
  
  #Run the ABC - log transform the parameters for the regression step. Use tolerance of 5% (= ~ top 1000 for 20,000 simulations)
  cat("Running the abc",sep="\n\n")
  abc.res<-abc(target = condense_sumstats(all.ss[[pair]][which_sumstats]),param = params[,6:7],sumstat=condense_sumstats(sumstats[,which_sumstats]),tol = 0.05,transf = c("log","log"),method="neuralnet")
  
  #How close are the closest fit??
  closest_fit<-which.min(abc.res$dist)
  as.matrix(sumstats[closest_fit,])
  
  data.ss.df<-data.frame(sumstat=names(all.ss[[pair]][which_sumstats]),value=all.ss[[pair]][which_sumstats])
  
  return(as.data.frame(abc.res$ss)%>%
           #gather(key="sumstat",value="value")%>%
           mutate(PairID=pair,.before=1))
})%>%bind_rows()%>%
  mutate(PairID=factor(new_pair_names[PairID],levels=paste("Pair",1:10,sep="_")))

#Get the data summary stats in a format for plotting with ggplot2
which_sumstats<-all_tree_sumstats
sumstats_df<-Map(ss=all.ss,pair=names(all.ss),function(ss,pair) {
  ss<-condense_sumstats(ss[which_sumstats])
  #ss<-ss[which_sumstats]
  df<-data.frame(sumstat=names(ss),value=ss)%>%mutate(PairID=pair,.before=1)
})%>%bind_rows()%>%tibble::remove_rownames()%>%
  mutate(PairID=factor(new_pair_names[PairID],levels=paste("Pair",1:10,sep="_")))
  

#Posterior checks plot (large plot look over all pairs and summary statistics)
ppc_plot<-as.data.frame(posterior_checks)%>%
  gather(-PairID,key="sumstat",value="value")%>%
  ggplot(aes(x=value))+
  geom_histogram()+
  facet_grid(PairID~sumstat,scales="free")+
  geom_vline(data = sumstats_df, aes(xintercept=value),col="red")+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(size=5))
ggsave(filename = paste0(plots_dir,"PPC_plot.pdf"),ppc_plot,width=15,height=15)

#visualize specifically the summary stats that are evidence of HSCT-specific selection
sumstats_to_visualize<-c("n_singletons_R","R_D_coals_pre_ratio","R_D_coals_peri_ratio","R_max_coals_in_single_clade")
IDs_to_visualize=c("Pair_9","Pair_7","Pair_5")
ppc_selection_sumstats_plot<-as.data.frame(posterior_checks)%>%
  gather(-PairID,key="sumstat",value="value")%>%
  dplyr::filter(sumstat%in%sumstats_to_visualize)%>%
  ggplot(aes(x=value))+
  geom_histogram()+
  facet_grid(PairID~sumstat,scales="free")+
  geom_vline(data = sumstats_df%>%dplyr::filter(sumstat%in%sumstats_to_visualize), aes(xintercept=value),col="red")+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(size=5),strip.text.y = element_text(size=5))
ggsave(filename = paste0(plots_dir,"ppc_selection_sumstats_plot.pdf"),ppc_selection_sumstats_plot,width=15,height=15)

#Visualize the 'n singletons in Recipient' and 'R D coals pre-HSCT ratio' sumstats as 2d density plot & position of data
#(or could do posterior as dot plot with high alpha??)
sumstats_to_visualize<-c("n_singletons_R","R_D_coals_pre_ratio")
IDs_to_visualize=c("Pair_7","Pair_5","Pair_10","Pair_6","Pair_3")
ppc_selection_sumstats_plot<-as.data.frame(posterior_checks)%>%
  dplyr::select(PairID,all_of(sumstats_to_visualize))%>%
  dplyr::filter(PairID%in%IDs_to_visualize)%>%
  ggplot(aes(x=get(sumstats_to_visualize[1]),y=get(sumstats_to_visualize[2])))+
  #geom_density_2d()+
  #stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_jitter(alpha=0.05,size=1,height = 0,width=0.1)+
  facet_wrap(~PairID,scales="fixed",nrow=1)+
  scale_y_continuous(limits=c(-1.5,3))+
  geom_point(data = sumstats_df%>%dplyr::filter(sumstat%in%sumstats_to_visualize & PairID%in%IDs_to_visualize)%>%pivot_wider(names_from = "sumstat",values_from="value"),size=1.5, shape=18,col="red")+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(size=6))+
  labs(x="Number of singletons in recipient phylogeny",y="Ratio of numbers of\npre-HSCT coalescences\n(Recipient:Donor)")
ggsave(filename = paste0(plots_dir,"ppc_selection_sumstats_plot_Esim.pdf"),ppc_selection_sumstats_plot,width=6,height=2)

sumstats_to_visualize<-c("DR_coals_peri","R_max_coals_in_single_clade")
IDs_to_visualize=c("Pair_9","Pair_5")
ppc_postHSCT_division_select_sumstats_plot<-as.data.frame(posterior_checks)%>%
  dplyr::select(PairID,all_of(sumstats_to_visualize))%>%
  dplyr::filter(PairID%in%IDs_to_visualize)%>%
  ggplot(aes(x=get(sumstats_to_visualize[1]),y=get(sumstats_to_visualize[2])))+
  #geom_density_2d()+
  #stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_jitter(alpha=0.05,size=1,height=0.5,width=0.5)+
  facet_wrap(~PairID,scales="free_x",nrow=1)+
  geom_point(data = sumstats_df%>%dplyr::filter(sumstat%in%sumstats_to_visualize & PairID%in% IDs_to_visualize)%>%pivot_wider(names_from = "sumstat",values_from="value"),size=1.5, shape=18,col="red")+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(size=6))+
  labs(x="Number of singletons in recipient phylogeny",y="Maximum number of peri-HSCT\ncoalescences in single clade")
ggsave(filename = paste0(plots_dir,"ppc_postHSCT_division_select_sumstats_plot_Esim.pdf"),ppc_postHSCT_division_select_sumstats_plot,width=6,height=2)

