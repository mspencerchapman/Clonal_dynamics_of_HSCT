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
                legend.text = element_text(size=7),
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

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in the trees list for the ABC
#These are the versions of the trees used for the ABC - see the 'Compile_data_objects.R' script for how these are produced
ABC.trees<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/trees_for_ABC.Rds"))

#========================================#
# SPECIFIC FUNCTIONS FOR THIS SCRIPT ####
#========================================#

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

#========================================#
# GET SUMMARY STATISTICS FROM THE DATA ####
#========================================#

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

data.ss<-Map(ss=all.ss,Pair=names(all.ss),function(ss,Pair) {
  data.frame(ss)%>%t()%>%as.data.frame()%>%
    mutate(Pair=Pair,.before=1)%>%
    mutate(ratio=log(R_coals_pre)-log(D_coals_pre))
})%>%dplyr::bind_rows()%>%
  mutate(type="Data")
#========================================#
# Review the summary statistics from models M1 and M2 ####
#========================================#

## Generate Extended Data Fig. 14a-b ----
m1_sumstats=lapply(Pair_metadata$Pair,function(this_pair) {
  #These files are generated as shown in the scripts in the folder 05 ABC_simulation_scripts/
  dat=read.delim(paste0(root_dir,"/data/ABC_simulation_results/Engrafting_cell_number_ABC_m1/output_",this_pair,"/stats_sample.txt"),header=T)
  return(dat%>%mutate(Pair=this_pair,.before=1))
})%>%dplyr::bind_rows()

m1_sumstat_plot<-m1_sumstats%>%
  mutate(ratio=log(R_coals_pre)-log(D_coals_pre),type="Posterior distribution")%>%
  bind_rows(data.ss)%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste0("Pair_",1:10)))%>%
  filter(Pair_new%in%paste0("Pair_",c(3,5,6,7,10)))%>%
  mutate(alpha=ifelse(type=="Data",1,0.1))%>%
  ggplot(aes(x=n_singletons_R,y=ratio,col=type,size=factor(type,levels=c("Posterior distribution","Data")),alpha=alpha))+
  geom_point(aes(shape=factor(type,levels=c("Posterior distribution","Data"))))+
  #scale_size_discrete(limits=factor(c(1,0.1)))+
  guides(col=guide_legend(override.aes = list(alpha=1,size=2)),shape="none",alpha="none",size="none")+
  facet_grid(~Pair_new)+
  scale_color_discrete(type=c("red","grey"))+
  scale_alpha_continuous(range=c(0.3,1))+
  scale_size_discrete(range=c(0.5,1.5))+
  scale_y_continuous(limits=c(-2,4))+
  scale_x_continuous(breaks=seq(0,100,50))+
  theme_bw()+
  my_theme+
  theme(legend.title = element_blank())+
  labs(x="Number of singletons in recipient phylogeny",
       y="Log ratio of numbers of\npre-HCT coalescences\n(Recipient:Donor)")

m2_sumstats=lapply(Pair_metadata$Pair,function(this_pair) {
  #These files are generated as shown in the scripts in the folder 05 ABC_simulation_scripts/
  dat=read.delim(paste0(root_dir,"/data/ABC_simulation_results/Engrafting_cell_number_ABC_m2/output_",this_pair,"/stats_sample.txt"),header=T)
  return(dat%>%mutate(Pair=this_pair,.before=1))
})%>%dplyr::bind_rows()

m2_sumstat_plot<-m2_sumstats%>%
  mutate(ratio=log(R_coals_pre)-log(D_coals_pre),type="Posterior distribution")%>%
  bind_rows(data.ss)%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste0("Pair_",1:10)))%>%
  filter(Pair_new%in%paste0("Pair_",c(3,5,6,7,10)))%>%
  mutate(alpha=ifelse(type=="Data",1,0.1))%>%
  ggplot(aes(x=n_singletons_R,y=ratio,col=type,size=factor(type,levels=c("Posterior distribution","Data")),alpha=alpha))+
  geom_point(aes(shape=factor(type,levels=c("Posterior distribution","Data"))))+
  #scale_size_discrete(limits=factor(c(1,0.1)))+
  guides(col=guide_legend(override.aes = list(alpha=1,size=2)),shape="none",alpha="none",size="none")+
  facet_grid(~Pair_new)+
  scale_color_discrete(type=c("red","grey"))+
  scale_alpha_continuous(range=c(0.3,1))+
  scale_size_discrete(range=c(0.5,1.5))+
  scale_y_continuous(limits=c(-2,4))+
  scale_x_continuous(breaks=seq(0,100,50))+
  theme_bw()+
  my_theme+
  theme(legend.title = element_blank())+
  labs(x="Number of singletons in recipient phylogeny",
       y="Log ratio of numbers of\npre-HCT coalescences\n(Recipient:Donor)")

Fig14ab<-gridExtra::arrangeGrob(grobs = list(m1_sumstat_plot,m2_sumstat_plot),nrow=2)
ggsave(filename = paste0(plots_dir,"Fig14a-b.pdf"),Fig14ab,width=7,height=4)

## Generate Extended Data Fig. 14c ----
m1_sumstat_plot2<-m1_sumstats%>%
  mutate(ratio=log(R_coals_pre)-log(D_coals_pre),type="Posterior distribution")%>%
  bind_rows(data.ss)%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste0("Pair_",1:10)))%>%
  filter(Pair_new%in%paste0("Pair_",c(5,9)))%>%
  mutate(alphav=ifelse(type=="Data",1,0.8))%>%
  ggplot(aes(x=n_singletons_R,y=R_max_coals_in_single_clade,col=type,size=factor(type,levels=c("Posterior distribution","Data")),alpha=alphav))+
  geom_point(aes(shape=factor(type,levels=c("Posterior distribution","Data"))))+
  guides(col=guide_legend(override.aes = list(alpha=1,size=3)),shape="none",alpha="none",size="none")+
  facet_grid(~Pair_new)+
  scale_alpha_continuous(range=c(0.5,1))+
  scale_size_discrete(range=c(0.5,2))+
  scale_color_discrete(type=c("red","grey"))+
  scale_x_continuous(breaks=seq(0,100,50))+
  theme_bw()+
  my_theme+
  theme(legend.title = element_blank())+
  labs(x="Number of singletons in recipient phylogeny",
       y="Maximum number of peri-HCT \ncoalescences in single recipient clade")

ggsave(filename = paste0(plots_dir,"Fig14c.pdf"),m1_sumstat_plot2,width=5,height=2)

#========================================#
# Review the p-values from the Posterior Predictive Checks ####
#========================================#

## Generate Extended Data Fig. 14d ----
all_models=c("m1","m2","m3")
PPC_pvalues=lapply(all_models,function(model) {
  #These files are generated as shown in the scripts in the folder 05 ABC_simulation_scripts/
  dat=read.delim(paste0(root_dir,"/data/ABC_simulation_results/Engrafting_cell_number_ABC_",model,"/PPCs/ppc_result_table_MODEL_",model,"_STAT_SET_original_R_and_D_stats_ABC_METHOD_rejection_RESIMS_PER_POST_OBS_1000.txt"),header=F)
  colnames(dat)<-c("idx","Pair","nWGS","nsim","x","y","z","p_value_full","p-value")
  return(dat%>%dplyr::select(-idx)%>%mutate(model=model))
})%>%dplyr::bind_rows()

all_plots<-lapply(all_models,function(this_model) {
  p<-PPC_pvalues%>%
    filter(model==this_model)%>%
    arrange(`p-value`)%>%
    mutate(idx=1:nrow(.))%>%
    mutate(Quantile=(idx-0.5)/nrow(.),Pair_new=new_pair_names[Pair])%>%
    ggplot(aes(x=Quantile,y=`p-value`,col=Pair_new))+
    geom_point()+
    scale_color_manual(values=Pair_cols)+
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1),expand = c(0,0))+
    scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1),expand = c(0,0))+
    geom_abline(col="skyblue")+
    geom_hline(yintercept = 0.05,col="darkgray",alpha=0.5)+
    theme_bw()+
    my_theme+
    theme(legend.title = element_blank(),panel.grid.minor=element_blank())+
    guides(col=guide_legend(overide.aes=list(size=3)))+
    labs(y="Posterior p value")
  
  if(this_model=="m3") {
    return(p+theme(axis.title.y = element_blank()))
  } else if(this_model=="m1"){
    return(p+theme(legend.position = "none"))
  } else {
    return(p+theme(axis.title.y = element_blank(),legend.position = "none"))
  }
})

plot_comb<-gridExtra::arrangeGrob(grobs=all_plots,ncol=3,widths = c(1.1,1,1.4))
ggsave(filename=paste0(plots_dir,"ExtDatFig14d.pdf"),plot_comb,width=7,height=2.4)
