#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest")
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

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

#Generate information regarding loss-of-Y in male samples from X and Y coverage data
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#Read in mutational signature extraction data
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)


#========================================#
# CALCULATE A RANGE OF PHYLOGENETIC DIVERSITY MEASURES TO COMPARE THE DONOR & RECIPIENT TREES ####
#========================================#

#Compare the 'sharedness' stats of each D & R trees
sharedness_stats<-dplyr::bind_rows(lapply(all.trees.ultra,function(tree) {
  #Extract the donor/ recipient Ids and the donor/ recipient trees
  donor_ID=get_DR_ids(tree)['donor_ID']
  recip_ID=get_DR_ids(tree)['recip_ID']
  D_tree<-keep.tip(tree,grep(donor_ID,tree$tip.label,value=T))
  R_tree<-keep.tip(tree,grep(recip_ID,tree$tip.label,value=T))
  
  #The sharedness statistic is sensitive to tree size, therefore subsample the larger tree so that it is the same size as the smaller tree
  min_samples<-min(sum(grepl(donor_ID,tree$tip.label)),sum(grepl(recip_ID,tree$tip.label)))
  if(length(D_tree$tip.label)>min_samples){D_tree<-keep.tip(D_tree,sample(D_tree$tip.label,size=min_samples))}
  if(length(R_tree$tip.label)>min_samples){R_tree<-keep.tip(R_tree,sample(R_tree$tip.label,size=min_samples))}
  
  #Calculate the stats
  comb_stat<-calculate_sharedness_stat(drop.tip(tree,"Ancestral"))
  D_stat<-calculate_sharedness_stat(D_tree)
  R_stat<-calculate_sharedness_stat(R_tree)
  
  return(data.frame(comb_stat=comb_stat,D_stat=D_stat,R_stat=R_stat))
}))

#Then view as a log2(FC) from R~D
sharedness_DvsR_FC_plot<-cbind(Pair_metadata$Pair_new,sharedness_stats)%>%
  dplyr::rename(Pair_new="Pair_metadata$Pair_new")%>%
  mutate(FC=R_stat/D_stat)%>%
  gather(key="D_or_R",value="Sharedness",-FC,-Pair_new,-comb_stat)%>%
  ggplot(aes(x=Pair_new,xend=Pair_new,y=0,yend=log2(FC),col=factor(ifelse(FC>1,1,0))))+
  geom_segment(arrow=arrow(angle=30,length=unit(1,"mm")))+
  scale_color_manual(values=c("#11a0aa", "#c82565"),guide="none")+
  scale_y_continuous(limits=c(-0.8,0.8))+
  geom_hline(yintercept=0,linetype=2)+
  labs(x="",y="log2 of the Fold Change")+
  theme_bw()+
  my_theme+
  theme(axis.line.y = element_blank(),axis.ticks.y=element_blank())+
  coord_flip()

ggsave(filename = paste0(plots_dir,"Sharedness_DvsR_plot.pdf"),sharedness_DvsR_FC_plot,width =2.5,height=2)

#Similar to sharedness, calculate the "Mean nearest taxon distance" for donor & recipient trees
#Mean nearest taxonomic distance
calculate_MNTD=function(tree) {
  tree<-drop.tip(tree,"Ancestral")
  tree$edge.length<-tree$edge.length/nodeheight(tree,1) #normalize the tree to overall height of 1 (so that each is comparable)
  cophen=cophenetic.phylo(tree)
  NTDs=apply(cophen,1,function(x) min(x[x>0]))
  return(mean(NTDs))
}
#Faith's Phylogenetic Diversity
calculate_FaithsPD=function(tree) {
  tree<-drop.tip(tree,"Ancestral")
  tree$edge.length<-tree$edge.length/nodeheight(tree,1) #normalize the tree to overall height of 1 (so that each is comparable)
  sum(tree$edge.length)
}
#Mean pairwise distance
calculate_MPD=function(tree) {
  require(ape)
  tree<-drop.tip(tree,"Ancestral")
  tree$edge.length<-tree$edge.length/nodeheight(tree,1) #normalize the tree to overall height of 1 (so that each is comparable)
  return(mean(cophenetic.phylo(tree)))
}

#1-sharedness as previously defined
calculate_one_minus_sharedness= function(tree) {
  tree<-drop.tip(tree,"Ancestral")
  tree$edge.length<-tree$edge.length/nodeheight(tree,1) #normalize the tree to overall height of 1 (so that each is comparable)
  
  prop_samples<-sapply(tree$edge[,2],function(node) {
    prop_samples<-length(getTips(tree,node))/length(tree$tip.label)
    return(prop_samples)
  })
  mean_w<-weighted.mean(x=prop_samples,w=tree$edge.length)
  return(1-mean_w) #Formulated here to calculate 1-sharedness for comparability to other diversity measures
}

#Shannon diversity index - defines 'clones' as originating from ancestors at 50 mutations of molecular time (i.e. early post-embryonic period)
calculate_SDI=function(tree,height_cut_off=50) {
  tree<-drop.tip(tree,"Ancestral")
  nodeheights=nodeHeights(tree)
  n_samples=length(tree$tip.label)
  #This pulls out branches that cross the 50 mutation mark
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off]
  
  clone_proportions=sapply(nodes,function(node) length(getTips(tree = tree,node=node))/n_samples)
  SDI=-sum(clone_proportions*log(clone_proportions)) #The Shannon diversity index is: - SUM p *log p
  return(SDI)
}

#Combine all functions into a single list that can be applied over the trees
diversity_functions=list(SDI=calculate_SDI,MPD=calculate_MPD,MNTD=calculate_MNTD,FaithsPD=calculate_FaithsPD,"1-Sharedness"=calculate_one_minus_sharedness)

#Input trees should be ultrametric, but not normalized - such that diversity index can correctly label clones post embryonic period
Diversity_stats_df=dplyr::bind_rows(Map(tree=all.trees.ultra,pair=names(all.trees.ultra),function(tree,pair) {
  tree<-drop.tip(tree,"Ancestral")
  
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

plot.MNTD.difference<-Diversity_stats_df%>%
  #mutate(Pair=factor(Pair,levels=Pair_metadata%>%arrange(Age)%>%pull(Pair)))%>%
  gather(-Pair,-D_or_R,key="Stat",value="value")%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=new_pair_names))%>%
  mutate(Age=sapply(Pair_new,function(pair) Pair_metadata$Age[Pair_metadata$Pair_new==pair]))%>%
  filter(Stat=="MNTD")%>%
  mutate(Stat=ifelse(Stat=="MNTD","Clonal diversity (MNTD)"))%>%
  mutate(D_or_R=ifelse(D_or_R=="D","Donor","Recipient"))%>%
  ggplot(aes(x=Pair_new,y=value,col=D_or_R))+
  geom_line(aes(group=Pair_new),col="black",arrow = grid::arrow(angle = 30,length = unit(2,"mm")))+
  geom_point(alpha=0.5,size=2)+
  #facet_grid(cols=vars(Stat),scales="free")+
  theme_bw()+
  my_theme+
  labs(x="",y="Clonal Diversity\n(Mean Nearest Taxon Distance)",title="",col="")+
  coord_flip()

plot.MNTD.MPD.FC<-Diversity_stats_df%>%
  gather(-Pair,-D_or_R,key="Stat",value="value")%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=new_pair_names))%>%
  mutate(Age=sapply(Pair_new,function(pair) Pair_metadata$Age[Pair_metadata$Pair_new==pair]))%>%
  filter(Stat%in%c("MNTD","MPD"))%>%
  #mutate(Stat=ifelse(Stat=="MNTD","Clonal diversity (MNTD)"))%>%
  mutate(D_or_R=ifelse(D_or_R=="D","Donor","Recipient"))%>%
  pivot_wider(names_from=c("D_or_R","Stat"),values_from="value")%>%
  mutate(MPD_FC=Recipient_MPD/Donor_MPD,MNTD_FC=Recipient_MNTD/Donor_MNTD)%>%
  dplyr::select(Pair_new,MPD_FC,MNTD_FC)%>%
  gather(-Pair_new,key="Stat",value="FC")%>%
  mutate(Stat=ifelse(Stat=="MNTD_FC","Mean Nearest\nTaxon Distance","Mean Pairwise\nDistance"))%>%
  ggplot(aes(x=Pair_new,xend=Pair_new,y=1,yend=FC,col=factor(ifelse(FC>1,1,0))))+
  geom_segment(arrow=arrow(angle=30,length=unit(1,"mm")))+
  scale_color_manual(values=c("#11a0aa", "#c82565"),guide="none")+
  #scale_y_continuous(limits=c(-1,1))+
  scale_y_log10()+
  facet_grid(~Stat,scales="free")+
  geom_hline(yintercept=1,linetype=2)+
  labs(x="",y="Diversity measure fold change\n(Recipent / Donor)")+
  theme_bw()+
  my_theme+
  theme(axis.line.y = element_blank(),axis.ticks.y=element_blank(),strip.text.x = element_text(size=6))+
  coord_flip()

ggsave(filename = paste0(plots_dir,"plot.MNTD.MPD.FC.pdf"),plot.MNTD.MPD.FC,width =2.5,height=2)

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

#Look at clones that are higher/ lower in donor vs recipient
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

#Add driver status of each expansion
expanded_df_full$driver_status<-sapply(1:nrow(expanded_df_full),function(i) {
  this_node<-expanded_df_full$nodes[i]
  pair<-expanded_df_full$Pair[i]
  details<-all.muts.nodups[[pair]]
  
  this_node_and_ancestors<-c(this_node,getAncestors(tree=all.trees.ultra[[pair]],node=this_node,type="all"))

  #Add the "loss of Y" details data frame for Pair 11, as these can be considered driver events
  if(pair=="Pair11"){
    details<-bind_rows(details,Pair11_loss_of_Y_details)
  }
  
  #Get the "driver status" of all mutations in that clone
  driver_status=details%>%
    filter(node%in%this_node_and_ancestors)%>%
    pull(coding_change_chip)
  return(any(driver_status=="yes"))
})

expanded_df_full%>%
  filter(D_frac>0.02|R_frac>0.02)%>%
  group_by(driver_status)%>%
  summarise(n_clones=n())

expanded_clones_DvsR_plot<-expanded_df_full%>%
  dplyr::select(Pair,driver_status,D_frac,R_frac)%>%
  gather(-Pair,-driver_status,key="DorR",value="clonal_fraction")%>%
  mutate(DorR=gsub("_frac","",DorR))%>%
  filter(clonal_fraction>0.02)%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  arrange(Pair_new,-clonal_fraction)%>% 
  ggplot(aes(x=DorR,y=clonal_fraction,fill=driver_status,group=forcats::fct_inseq(factor(clonal_fraction))))+
  geom_bar(stat="identity",position="stack",col="black",size=0.3)+
  geom_text(data=Pair_metadata,aes(x=1.5,y=0.95,label=paste0("Age=",Age)),size=2,inherit.aes = F)+
  facet_grid(cols=vars(Pair_new),drop = F)+
  scale_fill_brewer(palette = "Set2")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
  theme_classic()+
  my_theme+
  theme(strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")))+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction of expanded clones",fill="Clone has\nknown driver?")
ggsave(filename = paste0(plots_dir,"Expanded_clones_DvsR_plot.pdf"),expanded_clones_DvsR_plot,width=6,height=2)

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
# Pathway analysis on mutated genes that appear to be most 'selected for' in recipients vs donors ####
#========================================#

R_up_genes<-unlist(lapply(1:nrow(expanded_df_full), function(i) {
  if(expanded_df_full$p.value[i]>0.5|expanded_df_full$R_frac[i]<expanded_df_full$D_frac[i]){
    return(NULL)
  } else {
    details=all.muts[[expanded_df_full$Pair[i]]]
    tree=all.trees[[expanded_df_full$Pair[i]]]
    
    genes<-details$Gene[details$node==expanded_df_full$nodes[i] & details$coding_change!="no"]
    return(genes)
  }
}))

D_up_genes<-unlist(lapply(1:nrow(expanded_df_full), function(i) {
  if(expanded_df_full$p.value[i]>0.5|expanded_df_full$R_frac[i]>expanded_df_full$D_frac[i]){
    return(NULL)
  } else {
    details=all.muts[[expanded_df_full$Pair[i]]]
    tree=all.trees[[expanded_df_full$Pair[i]]]
    
    genes<-details$Gene[details$node==expanded_df_full$nodes[i] & details$coding_change!="no"]
    return(genes)
  }
}))

# library(AnnotationHub)
# ah <- AnnotationHub()
# HumanEnsDb <- query(ah, c("EnsDb", "Homo sapiens"))[[1]]
# annotations <- genes(HumanEnsDb, return.type = "data.frame")
# colnames(annotations)
# annot <- annotations %>%
#   dplyr::filter(!is.na(entrezid))%>%
#   dplyr::select(gene_id, gene_name, entrezid)
# 
# kegg_code<-search_kegg_organism('Homo sapiens', by='scientific_name')$kegg_code
# kk_R <- enrichKEGG(gene = annot%>%dplyr::filter(gene_name%in%R_up_genes)%>%pull(entrezid)%>%unlist(), organism = kegg_code)
# kk_D <- enrichKEGG(gene = annot%>%dplyr::filter(gene_name%in%D_up_genes)%>%pull(entrezid)%>%unlist(), organism = kegg_code)

expanded_clade_comparison_plot<-expanded_df_full%>%
  dplyr::filter(D_frac>0.02|R_frac>0.02)%>%
  #filter(p.value<0.1)%>%
  dplyr::select(-clonal_fraction)%>%
  mutate(cat=ifelse(p.value>0.05,"Similar fractions",ifelse(D_frac>R_frac,"Donor higher","Recip higher")))%>%
  gather(-nodes,-cat,-D_count,-R_count,-Pair,-D_total,-R_total,-p.value,key="D_or_R",value="clonal_fraction")%>%
  mutate(D_or_R=gsub("_frac","",D_or_R))%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  mutate(cat=factor(cat,levels=c("Recip higher","Donor higher","Similar fractions")))%>%
  ggplot(aes(x=D_or_R,y=clonal_fraction,col=factor(nodes),group=factor(nodes)))+
  geom_point(alpha=0.6)+
  geom_line(arrow=arrow(length=unit(2,"mm")),linetype=2)+
  facet_grid(cat~Pair_new,scales="free")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction")

ggsave(filename = paste0(plots_dir,"Expanded_clones_trajectory_plot.pdf"),expanded_clade_comparison_plot,width =12,height=5)
