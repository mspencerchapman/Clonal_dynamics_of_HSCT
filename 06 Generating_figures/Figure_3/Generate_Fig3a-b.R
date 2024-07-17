#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","tidyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest")
bioconductor_packages=c("clusterProfiler")

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
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

#========================================#
# Set the ggplot2 theme for plotting ####
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
trees_list<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/details_lists.Rds"))

## Extract objects from these lists in a 'for' loop ----
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

## Generate information regarding loss-of-Y in male samples from X and Y coverage data ----
LOY_files=list.files(path=paste0(root_dir,"/data/SV_and_CNA_data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Given that there are clearly many separate LOY events, manually annotate the nodes where these have been acquired
Pair11_LOY_nodes=c(447, 468, 152, 407, 369, 56, 493, 539)
Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                    Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                    Protein="",Type="",SO_codes="",coding_change="Coding change",
                                    coding_change_chip="yes",
                                    ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))

## Read in the spreadsheet listing other copy number changes ----
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

#Read in mutational signature extraction data
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)


#========================================#
# Do an "oligo-clonal" contribution diagram ####
#========================================#

#Define function to collect the expanded clades and their clonal fractions from the tree
get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        !tree$edge[,2]%in%1:length(tree$tip.label) &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction]
  df=data.frame(nodes=nodes,clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}

## Look at clones that are higher/ lower in donor vs recipient----
#Can we use these to get a handle on patterns of differential selection during transplant??
expanded_df=lapply(all.trees.ultra,get_expanded_clade_nodes,min_clonal_fraction=0.01,height_cut_off=50)
expanded_df_full=dplyr::bind_rows(Map(function(Pair,df) {
  print(Pair)
  if(nrow(df)==0) {stop(return(df))}
  tree<-all.trees.ultra[[Pair]]
  donor_ID=get_DR_ids(tree)['donor_ID']
  recip_ID=get_DR_ids(tree)['recip_ID']
  
  total_donor=sum(grepl(donor_ID,tree$tip.label))
  total_recip=sum(grepl(recip_ID,tree$tip.label))
  
  df$D_count=sapply(df$nodes,function(node) {sum(grepl(donor_ID,getTips(tree,node)))})
  df$R_count=sapply(df$nodes,function(node) {sum(grepl(recip_ID,getTips(tree,node)))})
  
  return(df%>%mutate(Pair=Pair,D_total=total_donor,R_total=total_recip)%>%mutate(D_frac=D_count/D_total,R_frac=R_count/R_total))
},Pair=names(expanded_df),df=expanded_df))

## Add driver status of each expansion----
expanded_df_full$driver_status<-sapply(1:nrow(expanded_df_full),function(i) {
  this_node<-expanded_df_full$nodes[i]
  pair<-expanded_df_full$Pair[i]
  details<-all.muts.nodups[[pair]]
  
  this_node_and_ancestors<-c(this_node,getAncestors(tree=all.trees.ultra[[pair]],node=this_node,type="all"))
  
  #Add the "loss of Y" details data frame for Pair 11, as these can be considered driver events
  if(pair=="Pair11"){details<-bind_rows(details,Pair11_loss_of_Y_details)}
  
  #Get the "driver status" of all mutations in that clone
  driver_status=details%>%
    filter(node%in%this_node_and_ancestors)%>%
    pull(coding_change_chip)
  return(ifelse(any(driver_status=="yes"),"Known driver","No known driver"))
})

expanded_df_full%>%
  filter(D_frac>0.02|R_frac>0.02)%>%
  group_by(driver_status)%>%
  summarise(n_clones=n())

expanded_clones_DvsR_plot<-expanded_df_full%>%
  dplyr::select(Pair,driver_status,D_frac,R_frac)%>%
  tidyr::gather(-Pair,-driver_status,key="DorR",value="clonal_fraction")%>%
  mutate(DorR=gsub("_frac","",DorR))%>%
  filter(clonal_fraction>0.02)%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")),
         driver_status=factor(driver_status,levels=c("No known driver","Known driver")))%>%
  arrange(Pair_new,-clonal_fraction)%>% 
  ggplot(aes(x=DorR,y=clonal_fraction,fill=driver_status,group=forcats::fct_inseq(factor(clonal_fraction))))+
  geom_bar(stat="identity",position="stack",col="black",linewidth=0.3)+
  geom_text(data=Pair_metadata,aes(x=1.5,y=0.95,label=paste0("Age=",Age)),size=2,inherit.aes = F)+
  facet_grid(cols=vars(Pair_new),drop = F)+
  scale_fill_brewer(palette = "Set2")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),panel.grid.minor = element_blank())+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction of\nexpanded clones",fill="")
ggsave(filename = paste0(plots_dir,"Fig3a.pdf"),expanded_clones_DvsR_plot,width=6.5,height=2)

#Calculate the mean change in oligoclonality in the older trees
expanded_df_full%>%
  dplyr::select(Pair,driver_status,D_frac,R_frac)%>%
  gather(-Pair,-driver_status,key="DorR",value="clonal_fraction")%>%
  mutate(DorR=gsub("_frac","",DorR))%>%
  filter(clonal_fraction>0.02)%>%
  group_by(Pair,DorR)%>%
  summarise(oligoclonality=sum(clonal_fraction),.groups="drop")%>%
  pivot_wider(id_cols="Pair",names_from="DorR",values_from="oligoclonality")%>%
  replace_na(list(D=0))%>%
  mutate(diff=R-D)%>%
  summarise(mean_diff=mean(diff))

#Add the p.value from proportions test to see which are robustly different between donor & recipient
expanded_df_full$p.value=sapply(1:nrow(expanded_df_full),function(i) {
  prop.test(x=c(expanded_df_full$D_count[i],expanded_df_full$R_count[i]),n=c(expanded_df_full$D_total[i],expanded_df_full$R_total[i]))$p.value
})

#========================================#
# CALCULATE A RANGE OF PHYLOGENETIC DIVERSITY MEASURES TO COMPARE THE DONOR & RECIPIENT TREES ####
#========================================#

#Similar to sharedness, calculate the "Mean nearest taxon distance" for donor & recipient trees
#Mean nearest taxonomic distance
calculate_MNTD=function(tree) {
  tree<-drop.tip(tree,"Ancestral")
  tree$edge.length<-tree$edge.length/nodeheight(tree,1) #normalize the tree to overall height of 1 (so that each is comparable)
  cophen=cophenetic.phylo(tree)
  NTDs=apply(cophen,1,function(x) min(x[x>0]))
  return(mean(NTDs))
}

#Mean pairwise distance
calculate_MPD=function(tree) {
  require(ape)
  tree<-drop.tip(tree,"Ancestral")
  tree$edge.length<-tree$edge.length/nodeheight(tree,1) #normalize the tree to overall height of 1 (so that each is comparable)
  return(mean(cophenetic.phylo(tree)))
}

#Combine all functions into a single list that can be applied over the trees
diversity_functions=list(MPD=calculate_MPD,MNTD=calculate_MNTD)

#Input trees should be ultrametric, but not normalized - such that diversity index can correctly label clones post embryonic period
Diversity_stats_df=dplyr::bind_rows(Map(tree=all.trees.ultra,pair=names(all.trees.ultra),function(tree,pair) {
  tree<-drop.tip(tree,"Ancestral") #Remove the ancestral branch
  
  #Extract the donor/ recipient Ids and the donor/ recipient trees
  donor_ID=get_DR_ids(tree)['donor_ID']; D_tree<-keep.tip(tree,grep(donor_ID,tree$tip.label,value=T))
  recip_ID=get_DR_ids(tree)['recip_ID']; R_tree<-keep.tip(tree,grep(recip_ID,tree$tip.label,value=T))
  
  #Diversity indices are sensitive to tree size (or 'richness'), therefore subsample the larger tree so that it is the same size as the smaller tree
  min_samples<-min(sum(grepl(donor_ID,tree$tip.label)),sum(grepl(recip_ID,tree$tip.label)))
  if(length(D_tree$tip.label)>min_samples){D_tree<-keep.tip(D_tree,sample(D_tree$tip.label,size=min_samples))}
  if(length(R_tree$tip.label)>min_samples){R_tree<-keep.tip(R_tree,sample(R_tree$tip.label,size=min_samples))}
  
  #Apply each stat function over the 
  D_stats<-sapply(diversity_functions,function(FUNC) {return(FUNC(D_tree))})
  R_stats<-sapply(diversity_functions,function(FUNC) {return(FUNC(R_tree))})
  
  return(data.frame(Pair=rep(pair,2),D_or_R=c("D","R"))%>%cbind(bind_rows(D_stats,R_stats)))
}))

plot.MNTD.MPD.FC<-Diversity_stats_df%>%
  gather(-Pair,-D_or_R,key="Stat",value="value")%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste0("Pair_",1:10)))%>%
  mutate(Age=sapply(Pair_new,function(pair) Pair_metadata$Age[Pair_metadata$Pair_new==pair]))%>%
  mutate(D_or_R=ifelse(D_or_R=="D","Donor","Recipient"))%>%
  pivot_wider(names_from=c("D_or_R","Stat"),values_from="value")%>%
  mutate(MPD_FC=Recipient_MPD/Donor_MPD,MNTD_FC=Recipient_MNTD/Donor_MNTD)%>%
  dplyr::select(Pair_new,MPD_FC,MNTD_FC)%>%
  gather(-Pair_new,key="Stat",value="FC")%>%
  mutate(Stat=ifelse(Stat=="MNTD_FC","Mean Nearest\nTaxon Distance","Mean Pairwise\nDistance"))%>%
  ggplot(aes(x=Pair_new,xend=Pair_new,y=1,yend=FC,col=factor(ifelse(FC>1,1,0))))+
  geom_segment(arrow=arrow(angle=30,length=unit(1,"mm")),linewidth = 0.4)+
  scale_color_manual(values=c("#11a0aa", "#c82565"),guide="none")+
  scale_y_log10()+
  facet_grid(~Stat,scales="free")+
  geom_hline(yintercept=1,linetype=2)+
  labs(x="",y="Diversity measure fold change\n(Recipent / Donor)")+
  theme_bw()+
  my_theme+
  theme(axis.line.y = element_blank(),axis.ticks.y=element_blank(),strip.text.x = element_text(size=7))+
  coord_flip()

ggsave(filename = paste0(plots_dir,"Fig3b.pdf"),plot.MNTD.MPD.FC,width =2.5,height=2)

