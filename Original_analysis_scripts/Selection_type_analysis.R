##START OF MAIN ANALYSIS SCRIPT
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

#Import the targeted sequencing metadata
bulk_smry_all<-readr::read_csv(file = paste0(root_dir,"/data/targeted_sequencing_data/targeted_sequencing_metadata.csv"))
##Define the pairs that have targeted sequencing results
all_pairs=unique(bulk_smry_all$Pair)

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

#Define the colour themes
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired")
names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]
names(DorR_cols)<-c("D","R")

#Read in other data objects
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
annotated_drivers_df<-read_csv(paste0(root_dir,"/data/Possible_drivers_annotated.csv"))[,c("mut_ref","node","Gene","variant_ID","Decision")]

consistency_fixed=T
if(consistency_fixed) {
  trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
  details_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))
  
  #Extract objects from these lists in a 'for' loop
  for(x in names(trees_list)) {assign(x,trees_list[[x]])}
  for(x in names(details_lists)) {assign(x,details_lists[[x]])}
  
  #Create dataframe of the LOY events
  Pair11_LOY_nodes=c(447, 468, 152, 407, 369, 56, 493, 539)
  Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                      Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                      Protein="",Type="",SO_codes="",coding_change="Coding change",
                                      coding_change_chip="yes",
                                      ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))
  
} else {
  tree_folder=paste0(root_dir,"/data/trees_no_dups")
  annotated_muts_folder=paste0(root_dir,"/data/annot_files_no_dups")
  tree_paths=list.files(tree_folder,pattern="a_j_vaf_post_mix_post_dup.tree",full.names = T)
  all.trees.cc.nodups<-lapply(all_pairs,function(PairID){read.tree(grep(PairID,tree_paths,value = T))})
  names(all.trees.cc.nodups)<-all_pairs
  
  #Create dataframe of the LOY events
  Pair11_LOY_nodes=c(48, 373,409 ,450 ,156 ,155 ,493 ,541,626)
  Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                      Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                      Protein="",Type="",SO_codes="",coding_change="Coding change",
                                      coding_change_chip="yes",
                                      ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))
}

####---------------------SELECTION TYPE ANALYSIS--------------------

#Classify 'type' of selection by looking at pattern of coalescences in D and R trees
expanded_clades_stats<-Map(comb_tree=all.trees.ultra,pair=names(all.trees.ultra),details=all.muts.nodups,function(comb_tree,pair,details){
  cat(pair,sep="\n")
  donor_id<-get_DR_ids(comb_tree)['donor_ID']
  recip_id<-get_DR_ids(comb_tree)['recip_ID']
  tree<-comb_tree
  age_of_donor_at_sampling<-Pair_metadata%>%filter(Pair==pair)%>%pull(Age)
  age_of_donor_at_HSCT<-Pair_metadata%>%filter(Pair==pair)%>%pull(Age_at_transplant)
  
  #Get a time tree to look at peri-HSCT
  tree.time<-tree
  tree.time$edge.length<-tree$edge.length*(age_of_donor_at_sampling/mean(get_mut_burden(tree)))
  D_tree.time<-keep.tip(tree.time,tip=grep(donor_id,tree$tip.label))
  R_tree.time<-keep.tip(tree.time,tip=grep(recip_id,tree$tip.label))
  
  #Find expanded driver clades
  if(pair=="Pair11") {
    Pair11_LOY_nodes=c(447, 468, 152, 407, 369, 56, 493, 541)
    Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                        Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                        Protein="",Type="",SO_codes="",coding_change="Coding change",
                                        coding_change_chip="yes",
                                        ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))
    Pair11_loss_of_Y_details$shared_coding_change_chip=sapply(Pair11_loss_of_Y_details$node,function(node) ifelse(length(getTips(comb_tree,node))>1,"yes","no"))
    details<-bind_rows(details,Pair11_loss_of_Y_details)
  }
  
  #Find the expanded clades at 100 muts of molecular time
  expanded_clades=get_expanded_clade_nodes(tree=tree,height_cut_off = 100,min_samples = 3)
  if(nrow(expanded_clades)==0) {stop(return(NULL))}
  
  #Find the driver clades
  expanded_driver_nodes<-details%>%filter(shared_coding_change_chip=="yes" & (Decision%in%c("Oncogenic","Possible")|grepl("LOY",Gene)))%>%dplyr::select(node,Gene,variant_ID)
  
  if(nrow(expanded_driver_nodes)>0) {
    expanded_driver_df<-expanded_driver_nodes%>%mutate(n_samples=sapply(node,function(this_node) length(getTips(tree,this_node))))%>%dplyr::rename("nodes"=node,"driver"=variant_ID)
    expanded_driver_df$clone_muts<-sapply(1:nrow(expanded_driver_df),function(i) {
      ancestor_nodes<-getAncestors(tree,node=expanded_driver_df$nodes[i],type="all")
      daughter_nodes<-get_all_node_children(node=expanded_driver_df$nodes[i],tree = tree)
      if(any(ancestor_nodes%in%expanded_driver_df$nodes)) {
        ancestral_drivers<-expanded_driver_df%>%filter(nodes%in%ancestor_nodes)%>%pull(driver)
        clone_drivers<-paste0(expanded_driver_df$driver[i],"/ ",paste0(ancestral_drivers,collapse="/ "))
      } else {
        clone_drivers<-expanded_driver_df$driver[i]
      }
      
      if(any(daughter_nodes%in%expanded_driver_df$nodes)) {
        daughter_drivers<-expanded_driver_df%>%filter(nodes%in%daughter_nodes)%>%pull(driver)
        clone_drivers<-paste0(clone_drivers," (",paste0(daughter_drivers,collapse=", "),")")
      }
      
      return(clone_drivers)
    })
    
    expanded_driver_df$clone_gene_muts<-sapply(1:nrow(expanded_driver_df),function(i) {
      ancestor_nodes<-getAncestors(tree,node=expanded_driver_df$nodes[i],type="all")
      daughter_nodes<-get_all_node_children(node=expanded_driver_df$nodes[i],tree = tree)
      if(any(ancestor_nodes%in%expanded_driver_df$nodes)) {
        ancestral_drivers<-expanded_driver_df%>%filter(nodes%in%ancestor_nodes)%>%pull(Gene)
        clone_drivers<-paste0(expanded_driver_df$Gene[i],"/ ",paste0(ancestral_drivers,collapse="/ "))
      } else {
        clone_drivers<-expanded_driver_df$Gene[i]
      }
      
      if(any(daughter_nodes%in%expanded_driver_df$nodes)) {
        daughter_drivers<-expanded_driver_df%>%filter(nodes%in%daughter_nodes)%>%pull(Gene)
        clone_drivers<-paste0(clone_drivers," (",paste0(daughter_drivers,collapse=", "),")")
      }
      
      return(clone_drivers)
    })
    expanded_clades<-full_join(expanded_clades,expanded_driver_df,by=c("nodes","n_samples"))
  }
  
  clade_samples=lapply(expanded_clades$nodes,function(node) getTips(tree,node=node))
  D_clade_samples=lapply(clade_samples,function(samples) {return(grep(donor_id,samples,value=T))})
  R_clade_samples=lapply(clade_samples,function(samples) {return(grep(recip_id,samples,value=T))})
  
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
                                    
                                    if(pair=="Pair11" & expanded_clades$nodes[i]==541) {
                                      if(expanded_clades$Gene[i]=="LOY") {
                                        TET2_node=expanded_clades%>%filter(driver=="TET2 p.G1275R")%>%pull(R_nodes)
                                        R_clade_nodes<-R_clade_nodes[!R_clade_nodes%in%get_all_node_children(node=TET2_node,tree=R_tree.time)] 
                                      }
                                    }
                                    
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
                                    
                                    ##Now need to perform the bootstraps
                                    n_bootstrap=500
                                    dummy_node=1000
                                    selection_stats_bootstraps<-lapply(1:n_bootstrap,function(k) {
                                      
                                      #Add a dummy node to each tree that is the theoretical "next node" added
                                      #This dummy node may be a another singleton/ a "pre-HSCT" node/ or a "peri-HSCT" node & the probabilities are weighted by existing numbers
                                      n_singletons.R<-sum(nodeHeights(R_tree.time)[,1]<5 & R_tree.time$edge[,2]%in%1:length(R_tree.time$tip.label))
                                      n_singletons.D<-sum(nodeHeights(D_tree.time)[,1]<5 & D_tree.time$edge[,2]%in%1:length(D_tree.time$tip.label))
                                      
                                      dummy_type.R=sample(size=1,x=c("developmental","pre","peri"),prob=c(max(n_singletons.R,1),max(n_pre.R,1),max(n_peri.R,1)))
                                      dummy_type.D=sample(size=1,x=c("developmental","pre","peri"),prob=c(max(n_singletons.D,1),max(n_pre.D,1),max(n_peri.D,1)))

                                      # dummy_type.R=sample(size=1,x=c("pre","peri"),prob=c(max(n_pre.R,1),max(n_peri.R,1)))
                                      # dummy_type.D=sample(size=1,x=c("pre","peri"),prob=c(max(n_pre.D,1),max(n_peri.D,1)))
                                      # 
                                      D_PDNs_resamples<-sample(c(D_all_PDNs,dummy_node),replace=T)
                                      R_PDNs_resamples<-sample(c(R_all_PDNs,dummy_node),replace=T)
                                      
                                      n_peri.D=sum(D_PDNs_resamples%in%c(D_clade_peri_nodes,ifelse(dummy_type.R=="peri",dummy_node,NA)))
                                      n_peri.R=sum(R_PDNs_resamples%in%c(R_clade_peri_nodes,ifelse(dummy_type.D=="peri",dummy_node,NA)))
                                      
                                      n_pre.D=sum(D_PDNs_resamples%in%c(D_clade_pre_nodes,ifelse(dummy_type.D=="pre",dummy_node,NA)))
                                      n_pre.R=sum(R_PDNs_resamples%in%c(R_clade_pre_nodes,ifelse(dummy_type.R=="pre",dummy_node,NA)))
                                      
                                      df<-data.frame(n_peri.D=n_peri.D,
                                                     n_peri.R=n_peri.R,
                                                     n_pre.D=n_pre.D,
                                                     n_pre.R=n_pre.R,
                                                     Type1_selection_stat=((n_pre.R+1)/(expanded_clades$R_total[i]+2))/((n_pre.D+1)/(2+expanded_clades$D_total[i])),
                                                     Type2_selection_stat=((n_peri.R+1)/(expanded_clades$R_total[i]+2))/((n_peri.D+1)/(2+expanded_clades$D_total[i])))
                                    })%>%dplyr::bind_rows()
                                    
                                    selection_stats_with_CI<-bind_cols(selection_stats_df,
                                                               Type1_lowerCI=quantile(selection_stats_bootstraps$Type1_selection_stat,0.025),
                                                               Type1_upperCI=quantile(selection_stats_bootstraps$Type1_selection_stat,0.975),
                                                               Type2_lowerCI=quantile(selection_stats_bootstraps$Type2_selection_stat,0.025),
                                                               Type2_upperCI=quantile(selection_stats_bootstraps$Type2_selection_stat,0.975))
                                    return(selection_stats_with_CI)
                                    
                                  })%>%dplyr::bind_rows())
  
  return(expanded_clades_full)
})

expanded_clades_stats_df<-expanded_clades_stats%>%
  dplyr::bind_rows()%>%
  mutate(id=paste(Pair,nodes,sep="_"))%>%
  mutate_at(c("Type1_selection_stat","Type2_selection_stat","Type1_lowerCI","Type1_upperCI","Type2_lowerCI","Type2_upperCI"),log)

confidence_interval_df<-lapply(1:nrow(expanded_clades_stats_df),function(i) {
  
  T1L<-expanded_clades_stats_df$Type1_lowerCI[i]
  T1U<-expanded_clades_stats_df$Type1_upperCI[i]
  T2L<-expanded_clades_stats_df$Type2_lowerCI[i]
  T2U<-expanded_clades_stats_df$Type2_upperCI[i]
  
  min_width=0.2
  if((T1U-T1L)<min_width){
    T1U<-mean(T1U,T1L)+(min_width/2)
    T1L<-mean(T1U,T1L)-(min_width/2)
  }
  
  if((T2U-T2L)<min_width){
    T2U<-mean(T2U,T2L)+(min_width/2)
    T2L<-mean(T2U,T2L)-(min_width/2)
  }
  
  x=c(T1L,expanded_clades_stats_df$Type1_selection_stat[i],T1U,expanded_clades_stats_df$Type1_selection_stat[i])
  y=c(expanded_clades_stats_df$Type2_selection_stat[i],T2U,expanded_clades_stats_df$Type2_selection_stat[i],T2L)
  
  return(data.frame(nodes=expanded_clades_stats_df$nodes[i],Pair=expanded_clades_stats_df$Pair[i],x=x,y=y))
})%>%bind_rows()%>%
  mutate(id=paste(Pair,nodes,sep="_"))%>%
  left_join(Pair_metadata)

include_ids<-expanded_clades_stats_df%>%
  mutate(clone_gene_muts=ifelse(clone_gene_muts=="LOY (TET2)","LOY (excluding\nTET2 subclone)",clone_gene_muts))%>%
  mutate(ratio=((1+n_R_samples)/(1+R_total))/((1+n_D_samples)/(1+D_total)))%>%
  filter((n_D_samples+n_R_samples)>9 & ratio>1.25)%>%
  mutate(id=paste(Pair,nodes,sep="_"))%>%pull(id)

plot.selection.type<-expanded_clades_stats_df%>%
  mutate(clone_gene_muts=ifelse(clone_gene_muts=="LOY (TET2)","LOY (excluding\nTET2 subclone)",clone_gene_muts))%>%
  mutate(ratio=((1+n_R_samples)/(1+R_total))/((1+n_D_samples)/(1+D_total)))%>%
  filter(id%in%include_ids)%>%
  left_join(Pair_metadata)%>%
  mutate(label=ifelse(is.na(driver),NA,driver))%>%
  ggplot(aes(x=Type1_selection_stat,y=Type2_selection_stat,col=Pair_new,label=clone_gene_muts))+
  geom_abline(slope = 1)+
  geom_point()+
  geom_polygon(data=confidence_interval_df%>%filter(id%in%include_ids),aes(x=x,y=y,fill=Pair_new,group=id),alpha=0.2,inherit.aes = F)+
  scale_x_continuous(limits=c(-1,5.1))+
  scale_y_continuous(limits=c(-1,4.6))+
  ggrepel::geom_label_repel(size=1.5,show.legend = F)+
  theme_classic()+
  my_theme+
  #geom_vline(xintercept=0)+
  scale_color_manual(values=Pair_cols)+
  scale_fill_manual(values=Pair_cols)+
  labs(col="Pair",fill="Pair",x="Log of the Type 1 selection statistic",y="Log of the Type 2 selection statistic")

ggsave(filename = paste0(plots_dir,"plot.selection.type.pdf"),plot.selection.type,width=4,height=3)

plot.bar.selection.type<-expanded_clades_stats%>%
  dplyr::bind_rows()%>%
  mutate(clone_gene_muts=ifelse(clone_gene_muts=="LOY (TET2)","LOY (excluding TET2 subclone)",clone_gene_muts))%>%
  mutate(ratio=((1+n_R_samples)/(1+R_total))/((1+n_D_samples)/(1+D_total)))%>%
  filter((n_D_samples+n_R_samples)>5 & ratio>1)%>%
  left_join(Pair_metadata)%>%
  mutate(nodes=factor(nodes))%>%
  ggplot(aes(x=nodes,fill=Pair_new))+
  geom_bar(aes(y=-Type2_selection_stat),stat="identity")+
  geom_bar(aes(y=Type1_selection_stat),stat="identity")+
  geom_text(aes(label=clone_gene_muts,y=Type1_selection_stat+15),size=1.5,angle=90)+
  scale_fill_manual(values=Pair_cols)+
  facet_grid(cols=vars(Pair_new),drop=T,scales="free_x",space="free")+
  theme_classic()+
  my_theme+
  theme(strip.text.x = element_text(angle=90),axis.text.x = element_blank())+
  geom_hline(yintercept = 0,linetype=2)+
  labs(x="Clone",y="Type 1 selection stat (positive) +\nType 2 selection stat (negative)")

ggsave(filename = paste0(plots_dir,"plot.bar.selection.type.pdf"),plot.bar.selection.type,width=4,height=3)
