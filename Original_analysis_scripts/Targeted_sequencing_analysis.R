##START OF MAIN ANALYSIS SCRIPT
library(stringr)
library(ape)
library(seqinr)
library(data.table)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)


root_dir<-ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/Clonal_dynamics_of_HSCT","/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT")
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
all_pairs=sort(unique(bulk_smry_all$Pair))

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=paste("Pair",1:nrow(Pair_metadata),sep = "_") #Useful named vector to swap from old to new names
names(new_pair_names)<-Pair_metadata%>%arrange(Age)%>%pull(Pair)

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
annotated_drivers_df<-read_csv(paste0(root_dir,"/data/Possible_drivers_annotated.csv"))[,c("mut_ref","node","Gene","variant_ID","Decision")]
remove_names=function(x) {attr(x,"names")<-NULL;return(x)}

consistency_fixed=F
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
  tree_paths=list.files(tree_folder,pattern="vaf_post_mix_post_dup.tree",full.names = T)
  all.trees.cc.nodups<-lapply(all_pairs,function(PairID){read.tree(grep(PairID,tree_paths,value = T))})
  names(all.trees.cc.nodups)<-all_pairs
  
  annotated_muts_paths=list.files(annotated_muts_folder,pattern="vaf_post_mix_post_dup",full.names = T)
  all.muts.nodups<-lapply(all_pairs,function(PairID){load(grep(PairID,annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
  names(all.muts.nodups)<-all_pairs
  
  #Create dataframe of the LOY events
  Pair11_LOY_nodes=c(48, 373,409 ,450 ,156 ,155 ,493 ,541,626)
  Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                      Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                      Protein="",Type="",SO_codes="",coding_change="Coding change",
                                      coding_change_chip="yes",
                                      ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))
}


#Import the raw targeted sequencing results - count data across all samples/ mutations sites
all_targeted_res<-readRDS(paste0(root_dir,"/data/Targeted_sequencing_data/all_targeted_res.Rds"))

#--------------COVERAGE PLOT--------------------
Targeted_coverage_plot<-lapply(all_targeted_res,function(list) as.data.frame(colMeans(list$matrices$NR)))%>%
  dplyr::bind_rows()%>%
  tibble::rownames_to_column(var="tissueID")%>%
  dplyr::rename("mean_cov"=2)%>%
  left_join(bulk_smry_all%>%dplyr::select(-sampleID)%>%filter(!duplicated(.)))%>%
  dplyr::filter(time_point==0)%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  ggplot(aes(x=cell_type,y=mean_cov,fill=cell_type))+
  geom_bar(stat="identity")+
  geom_hline(yintercept = 1000,linetype=2,col="black",size=0.4)+
  scale_fill_brewer(palette="Set2")+
  facet_grid(individual_type~factor(Pair_new,levels=paste("Pair",1:10,sep="_")))+
  scale_y_continuous(breaks=seq(1000,6000,1000))+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),
        strip.text.y = element_text(size=6,margin=unit(c(0.5,1,0.5,1),"mm")),
        strip.text.x = element_text(size=6,margin=unit(c(0.5,1,0.5,1),"mm")))+
  labs(fill="Cell type",x="Cell type",y="Mean targeted coverage")
ggsave(filename = paste0(plots_dir,"Targeted_coverage_plot.pdf"),Targeted_coverage_plot,width=7,height=3)

#----------------------GET CELL FRACTIONS FROM THE MODEL------------------------
#Import the output of the phylogeny-aware bayesian model that infers the posterior distribution of the cell fraction of each mutation
posterior_cell_fracs_file=paste0(root_dir,"/data/Targeted_sequencing_data/posterior_cell_fracs.Rds")
if(file.exists(posterior_cell_fracs_file)) {
  cat("Reading in previously saved RDS file of the posterior cell fractions",sep="\n")
  posterior_cell_fracs<-readRDS(posterior_cell_fracs_file)
} else {
  cat("Importing posterior VAFs files",sep="\n")
  setwd(paste0(root_dir,"/data/Targeted_sequencing_data/raw_model_output"))
  posterior_VAFs<-lapply(all_pairs,function(pair) {
    cat(pair,sep="\n")
    tissueIDs<-bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(tissueID)%>%unique()
    pair_posteriors<-lapply(tissueIDs,function(ID) {
      cat(ID,sep="\n")
      posterior_VAFs_file=paste0(ID,"/",ID,"_posterior_VAFs.txt")
      tissue_posteriors<-readr::read_delim(posterior_VAFs_file,delim = "\t",show_col_types = FALSE)
      return(tissue_posteriors)
    })
    names(pair_posteriors)<-tissueIDs
    return(pair_posteriors)
  })
  names(posterior_VAFs)<-all_pairs
  
  #Convert the VAF posterior to a clonal fraction posterior
  #Multiply VAFs of diploid chromosome mutations two
  cat("Converting VAFs to cell fractions",sep="\n")
  pair_sex=sapply(Pair_metadata$Pair,function(pair) {Pair_metadata$donor_sex[Pair_metadata$Pair==pair]})
  autosomes=as.character(1:22)
  posterior_cell_fracs<-Map(pair_posteriors=posterior_VAFs,sex=pair_sex,pair=all_pairs,function(pair_posteriors,sex,pair) {
    cat(pair,sep="\n")
    tissueIDs=names(pair_posteriors)
    pair_cell_fracs_posteriors<-Map(tissue_post=pair_posteriors,tissueID=tissueIDs,function(tissue_post,tissueID) {
      cat(tissueID,sep="\n")
      #Select & separate out the posterior values into individual columns
      VAF_distributions=tissue_post%>%
        dplyr::select(Posterior_VAFs)%>%
        separate(col="Posterior_VAFs",into=paste("n",1:100,sep="_"),sep=",")%>%
        mutate_all(as.numeric)
      
      if(sex=="Male"&F){ #The VAFs of XY mutations in males have now been divided by 2 in the original algorithm
        #For males, identify the XY mutations and only multiply the autosomal mutaiton VAFs by 2
        cat("Treating sample as male",sep="\n")
        Chrom=tissue_post%>%
          dplyr::mutate(Chrom=stringr::str_split(mutation_ID,pattern="-",simplify=T)[,1])%>%
          dplyr::pull(Chrom)
        cell_frac_distributions=VAF_distributions
        cell_frac_distributions[Chrom%in%autosomes,]<-2*VAF_distributions[Chrom%in%autosomes,]
      } else if(sex=="Female"|T){
        cat("Treating sample as female",sep="\n")
        cell_frac_distributions=2*VAF_distributions
      }
      cell_frac_post<-bind_cols((tissue_post%>%dplyr::select(Node_assignment,mutation_ID)),cell_frac_distributions)
      return(cell_frac_post)
    })
    return(pair_cell_fracs_posteriors)
  })
  names(posterior_cell_fracs)<-all_pairs
  
  saveRDS(posterior_cell_fracs,file = posterior_cell_fracs_file)
  
}

#Get VAFs of all clones at 100 mutations of molecular time
get_cutoff_branches=function(tree,cut_off) {
  heights=nodeHeights(tree)
  cutoff_branches=tree$edge[,2][heights[,1]<cut_off & heights[,2]>=cut_off]
  return(cutoff_branches)
}

#Get the posterior distributions for each 'clone' in each cell type defined by cutting tree at specified height
clone_cutoff=100
clone_sizes=Map(tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=all_pairs,function(tree,post,pair) {
  cat(pair,sep="\n")
  cutoff_branches=get_cutoff_branches(tree,cut_off = clone_cutoff)
  tissueIDs=names(post)
  tissue_clone_fractions<-lapply(tissueIDs,function(tissueID) {
    cat(tissueID,sep="\n")
    tissue_post=post[[tissueID]]
    clone_posteriors<-lapply(cutoff_branches,function(node) {
      branch_heights=nodeHeights(tree)[tree$edge[,2]==node,]
      min_height=branch_heights[1]
      max_height=branch_heights[2]
      prop=(clone_cutoff-min_height)/(max_height-min_height)
      n_targseq_muts=sum(tissue_post$Node_assignment==node)
      which_branch_rank=max(round(prop*n_targseq_muts),1)
      node_cell_frac_distributions=tissue_post%>%
        filter(Node_assignment==node)%>%
        dplyr::select(-Node_assignment,-mutation_ID)%>%
        as.matrix()
      
      #Sort each iteration of posterior by cell fraction & choose that which corresponds to the rank
      if(nrow(node_cell_frac_distributions)>1) {
        #Original method - sorts cell fracs of each Gibbs sample, and takes the 10th rank of each. ?significant overestimates
        #node_cell_frac_distributions_sorted<-apply(node_cell_frac_distributions,2,sort,decreasing=T)
        
        #Second method - keeps values of specific mutations linked together & just orders them by their median values
        rank=order(apply(node_cell_frac_distributions,1,median),decreasing = T)
        node_cell_frac_distributions_sorted<-node_cell_frac_distributions[rank,]
        
      } else {
        node_cell_frac_distributions_sorted<-node_cell_frac_distributions
        }
      
      clone_post=node_cell_frac_distributions_sorted[which_branch_rank,]
      return(clone_post)
    })
    names(clone_posteriors)<-paste("node",cutoff_branches,sep="_")
    return(clone_posteriors)
  })
  names(tissue_clone_fractions)<-tissueIDs
  return(tissue_clone_fractions)
})
names(clone_sizes)<-all_pairs

#Plot clone sizes across cells types
clones_df=Map(pair=clone_sizes,name=names(clone_sizes),function(pair,name) {
  tissues=names(pair)
  pair_fracs<-Map(tissue=pair,ID=tissues,function(tissue,ID) {
    data.frame(ID=ID,node=names(tissue),cell_frac=sapply(tissue,median))
  })%>%dplyr::bind_rows()%>%mutate(Pair=name)
  return(pair_fracs)
})

##PLOT OUT THE EXPANDED CLONES

colfunc=function(x,pal=c("#58b5e1", "#18786a", "#81bc31", "#c76952", "#0cc0aa", "#d10f55", "#a7af87", "#4d6c94", "#f587c8", "#0ec61c", "#f24219", "#496803", "#f79302", "#856619")) {
  set.seed(53)
  col_gen_func<-colorRampPalette(pal) #generate function to get colours based on palette
  cols<-col_gen_func(x)
  return(sample(cols)) #Mix them up so that similar colours aren't adjacent
}

clone_cols<-colfunc(length(unique(dplyr::bind_rows(clones_df)$node)))
names(clone_cols)<-unique(dplyr::bind_rows(clones_df)$node)

#Annotate clones info with whether contains a driver
clones_df<-Map(df=clones_df,tree=all.trees.cc.nodups,res=all_targeted_res,pair=names(all_targeted_res),function(df,tree,res,pair) {
  driver_details<-res$details_targ%>%filter(mut_ref%in%(annotated_drivers_df%>%filter(Decision%in%c("Oncogenic","Possible"))%>%pull(mut_ref)))
  if(pair=="Pair11") {
    driver_details<-bind_rows(driver_details,Pair11_loss_of_Y_details)
  }
  driver_nodes=driver_details%>%pull(node)
  driver_nodes_and_daughters=paste("node",unlist(lapply(driver_nodes,function(node) c(node,get_all_node_children(node = node,tree = tree)))),sep = "_")
  df$driver=df$node%in%driver_nodes_and_daughters
  return(df)
})

clones_100_all_plot<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  #mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(x=factor(cell_type),y=cell_frac,fill=node))+
  geom_bar(stat="identity",col="black",position="stack",size=0.05)+
  scale_fill_manual(values=clone_cols)+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),legend.position = "none",strip.text.x = element_text(size=6,margin = unit(c(0.7,1,0.7,1),"mm")))+
  facet_wrap(Pair_new~individual_type,nrow=1)+
  #facet_grid(Pair_new~individual_type)+
  scale_y_continuous(limits=c(0,1.1),breaks = seq(0,1,0.2))+
  labs(x="Cell type",y="Clonal fraction")
ggsave(filename = paste0(plots_dir,"clones_100_all_plot.pdf"),clones_100_all_plot,device = "pdf",width=7,height = 2.5)

#Now only include those clones with a fraction over 1%in donor or recipient
include_nodes_0.01<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  #filter(cell_type%in%c("Monocytes","Granulocytes"))%>%
  mutate(uid=paste(node,Pair,sep="_"))%>%
  filter(cell_frac>0.01)%>%
  arrange(cell_frac)%>%
  pull(uid)%>%
  unique(fromLast=T)
length(include_nodes_0.01) #114 total clones over 1% clonal fraction in at least one compartment

#What proportion harbour drivers - 17% = 19/114
sum(sapply(include_nodes_0.01,function(node) {
  this_node<-as.integer(stringr::str_split(node,pattern="_")[[1]][2])
  this_pair<-stringr::str_split(node,pattern="_")[[1]][3]
  this_tree<-all.trees.cc.nodups[[this_pair]]
  all_nodes<-c(this_node,getAncestors(this_tree,this_node,type="all"),get_all_node_children(this_node,this_tree))
  all_nodes<-all_nodes[!all_nodes%in%1:length(this_tree$tip.label)]
  if(this_pair=="Pair11") {
    details<-bind_rows(all_targeted_res[[this_pair]]$details_targ,Pair11_loss_of_Y_details)
  } else {
    details<-all_targeted_res[[this_pair]]$details_targ
  }
  node_muts<-details%>%filter(node%in%all_nodes)%>%pull(mut_ref)
  any(node_muts%in%(annotated_drivers_df%>%filter(Decision%in%c("Oncogenic","Possible"))%>%pull(mut_ref))|grepl("LOY",node_muts))
}))

clones_100_over1inDorR_plot<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  mutate(individual_type=substr(individual_type,1,5))%>%
  mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  mutate(size_vec=ifelse(driver,0.3,0.05))%>%
  filter(uid%in%include_nodes_0.01)%>%
  ggplot(aes(x=factor(cell_type),y=cell_frac,fill=node,col=driver))+
  geom_bar(lwd=0.15,stat="identity",position="stack")+
  theme_classic()+
  my_theme+
  scale_fill_manual(values=clone_cols)+
  scale_color_manual(values=c("black","red"))+
  guides(fill="none",linewidth="none")+
  theme(axis.text.x = element_text(angle=90),panel.spacing = unit(0.5, "mm"),strip.text.x = element_text(size=5,margin = unit(c(0.7,1,0.7,1),"mm")))+
  facet_wrap(Pair_new~individual_type,nrow=1)+
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  labs(x="Cell type",y="Clonal fraction",col="Clone contains\n driver?")
ggsave(filename = paste0(plots_dir,"clones_100_over1percinDorR_driver_colored_plot.pdf"),clones_100_over1inDorR_plot,device = "pdf",width=7,height = 2.5)

clones_100_over1_summed_DorR_barplot<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  filter(uid%in%include_nodes_0.01)%>%
  group_by(Pair_new,individual_type,cell_type)%>%
  dplyr::summarise(expanded_clone_cell_frac=sum(cell_frac))%>%
  ggplot(aes(x=factor(cell_type),y=expanded_clone_cell_frac,fill=cell_type,col=driver))+
  geom_bar(stat="identity",position="dodge",size=0.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),legend.position = "none",strip.text.x = element_text(size=6,margin = unit(c(0.7,1,0.7,1),"mm")))+
  facet_wrap(Pair_new~individual_type,nrow=1)+
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  labs(x="Cell type",y="Clonal fraction")

#Now only include those clones with a fraction over 0.5%in donor or recipient
include_nodes_0.005<-clones_df%>%
  dplyr::bind_rows()%>%
  mutate(uid=paste(node,Pair,sep="_"))%>%
  filter(cell_frac>0.005)%>%
  arrange(cell_frac)%>%
  pull(uid)%>%
  unique(fromLast=T)
length(include_nodes_0.005) #205 total clones over 0.5% clonal fraction in at least one compartment

clones_100_over0.5inDorR_plot<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  mutate(individual_type=substr(individual_type,1,5))%>%
  mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  filter(uid%in%include_nodes_0.005)%>%
  ggplot(aes(x=factor(cell_type),y=cell_frac,fill=node,col=driver))+
  geom_bar(stat="identity",position="stack",size=0.15)+
  theme_classic()+
  my_theme+
  scale_fill_manual(values=clone_cols)+
  scale_color_manual(values=c("black","red"))+
  guides(fill="none",linewidth="none")+
  theme(axis.text.x = element_text(angle=90),panel.spacing = unit(0.5, "mm"),strip.text.x = element_text(size=5,margin = unit(c(0.7,1,0.7,1),"mm")))+
  facet_wrap(Pair_new~individual_type,nrow=1)+
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  labs(x="Cell type",y="Clonal fraction")
ggsave(filename = paste0(plots_dir,"clones_100_over0.5inDorR_driver_colored_plot.pdf"),clones_100_over0.5inDorR_plot,device = "pdf",width=7,height = 2.5)

clones_100_over0.5_summed_DorR_barplot<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  filter(uid%in%include_nodes_0.005)%>%
  group_by(Pair_new,individual_type,cell_type)%>%
  dplyr::summarise(expanded_clone_cell_frac=sum(cell_frac))%>%
  ggplot(aes(x=factor(cell_type),y=expanded_clone_cell_frac,fill=cell_type))+
  geom_bar(stat="identity",col="black",position="dodge",size=0.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),legend.position = "none",strip.text.x = element_text(size=6,margin = unit(c(0.7,1,0.7,1),"mm")))+
  facet_wrap(Pair_new~individual_type,nrow=1)+
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  labs(x="Cell type",y="Clonal fraction\n(expanded clones)",fill="Cell type")
ggsave(filename = paste0(plots_dir,"clones_100_over0.5_summed_DorR_barplot.pdf"),clones_100_over0.5_summed_DorR_barplot,device = "pdf",width=7,height = 2)


clones_100_over0.5_summed_DorR_dotplot<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  dplyr::mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  dplyr::mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  dplyr::mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  dplyr::filter(uid%in%include_nodes_0.005)%>%
  group_by(Pair_new,individual_type,cell_type)%>%
  dplyr::summarise(expanded_clone_cell_frac=sum(cell_frac))%>%
  ggplot(aes(x=individual_type,y=expanded_clone_cell_frac,col=factor(cell_type)))+
  geom_point(alpha=0.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = unit(c(0.7,1,0.7,1),"mm")))+
  facet_wrap(~Pair_new,nrow=1)+
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  labs(x="",y="Clonal fraction (expanded clones)",col="Cell type")

clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  filter(uid%in%include_nodes_0.01)%>%
  group_by(Pair_new,individual_type,cell_type)%>%
  dplyr::summarise(expanded_clone_cell_frac=sum(cell_frac))%>%
  pivot_wider(id_cols=c("Pair_new","individual_type"),names_from=c("cell_type"),values_from="expanded_clone_cell_frac")%>%
  mutate('T_cell/ Mono'=T_cells/Monos,'B_cell/ Mono'=B_cells/Monos)%>%
  dplyr::select(-Monos,-B_cells,-T_cells,-Grans)%>%
  gather(-Pair_new,-individual_type,key="Measure",value="Fraction_of_monocytes")%>%
  ggplot(aes(x=individual_type,y=Fraction_of_monocytes,col=Measure))+
  geom_point()+
  scale_y_continuous(limits=c(0,1))+
  facet_grid(~Pair_new)+
  theme_classic()+
  my_theme+
  labs(x="",y="Fraction from clonal expansions\n(as proportion of fraction from expansions\nin monocytes)",col="")

clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  filter(uid%in%include_nodes_0.01)%>%
  group_by(Pair_new,individual_type,cell_type)%>%
  dplyr::summarise(expanded_clone_cell_frac=sum(cell_frac),.groups="drop")%>%
  pivot_wider(names_from=c("cell_type"),values_from="expanded_clone_cell_frac")%>%
  mutate('T_cell/ Mono'=T_cells/Monos,'B_cell/ Mono'=B_cells/Monos)%>%
  dplyr::select(-Monos,-B_cells,-T_cells,-Grans)%>%summarise(mean_Tfrac=mean(`T_cell/ Mono`),mean_Bfrac=mean(`B_cell/ Mono`),min_Tfrac=min(`T_cell/ Mono`),min_Bfrac=min(`B_cell/ Mono`),max_Tfrac=max(`T_cell/ Mono`),max_Bfrac=max(`B_cell/ Mono`))

#Compare clone fractions between phylogeny and targeted sequencing
phylo_clones_df<-Map(df=clones_df,tree=all.trees.cc.nodups,pair=names(clones_df),function(df,tree,pair) {
  donor_ID=get_DR_ids(tree)['donor_ID']
  recip_ID=get_DR_ids(tree)['recip_ID']
  
  total_donor=sum(grepl(donor_ID,tree$tip.label))
  total_recip=sum(grepl(recip_ID,tree$tip.label))
  
  D_count=sapply(df$node,function(this_node) {node1<-as.numeric(str_extract(this_node,"[0-9]+"));sum(grepl(donor_ID,getTips(tree,node1)))})
  R_count=sapply(df$node,function(this_node) {node1<-as.numeric(str_extract(this_node,"[0-9]+"));sum(grepl(recip_ID,getTips(tree,node1)))})
  
  D_phylo=Map(N=D_count,this_node=names(D_count),function(N,this_node) {
    res<-binom.test(N,total_donor)
    return(data.frame(ID="D_phylo",node=this_node,cell_frac=N/total_donor,Pair=pair,lower_CI=res$conf.int[1],upper_CI=res$conf.int[2],cell_type="Phylo",individual_type="Donor"))
    })%>%bind_rows()
  R_phylo=Map(N=R_count,this_node=names(R_count),function(N,this_node) {
    res<-binom.test(N,total_recip)
    return(data.frame(ID="R_phylo",node=this_node,cell_frac=N/total_recip,Pair=pair,lower_CI=res$conf.int[1],upper_CI=res$conf.int[2],cell_type="Phylo",individual_type="Recipient"))
  })%>%bind_rows()
  return(rbind(D_phylo,R_phylo))
})%>%dplyr::bind_rows()

include_nodes<-clones_df%>%
  dplyr::bind_rows()%>%
  mutate(uid=paste(node,Pair,sep="_"))%>%
  filter(cell_frac>0.05)%>%
  arrange(cell_frac)%>%
  #filter(!duplicated(uid))%>%
  pull(uid)%>%
  unique(fromLast=T)

plot_phylo_vs_targeted<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  dplyr::bind_rows(phylo_clones_df)%>%
  mutate(Pair_new=new_pair_names[Pair],uid=paste(node,Pair,sep="_"))%>%
  filter(!duplicated(.) & uid%in%include_nodes)%>%
  mutate(cell_type=factor(cell_type,levels=c("Phylo","Monocytes","Granulocytes","B_cells","T_cells")))%>%
  mutate(cell_frac=ifelse(cell_type=="Phylo",NA,cell_frac))%>%
  mutate(uid=factor(uid,levels=include_nodes))%>%
  ggplot(aes(x=uid,y=cell_frac,ymin=lower_CI,ymax=upper_CI,col=cell_type))+
  geom_errorbar(width=0.35)+
  geom_point(alpha=0.75,size=0.75)+
  facet_grid(rows=vars(individual_type))+
  scale_color_manual(values=c("gray20",brewer.pal(4,"Set2")))+
  theme_bw()+
  my_theme+
  #theme(axis.text.x=element_text(angle=90))+
  theme(axis.text.x=element_blank())+
  labs(x="Clone",y="Clonal fraction",col="Sample type")
ggsave(filename = paste0(plots_dir,"plot_phylo_vs_targeted.pdf"),plot_phylo_vs_targeted,device = "pdf",width=6,height = 4)

#Now look specifically at change in clone size between donors/ recipients
clones_logFC_Mono<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  pivot_wider(id_cols=c("node","Pair_new","cell_type"),names_from = "individual_type",values_from="cell_frac")%>%
  filter(cell_type=="Monocytes")%>%
  filter(Donor>0.01|Recipient>0.01)%>% #Include clones if expanded in the Donor OR the Recipient
  mutate(log2FC=log2(Recipient/Donor))%>%
  mutate(node=paste(node,Pair_new,sep="_"))%>%
  mutate(Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  mutate(bigger_in=ifelse(log2FC>0,"Recipient > Donor","Donor > Recipient"))%>%
  ggplot(aes(y=log2FC,x=forcats::fct_reorder(node,log2FC),fill=bigger_in))+
  geom_bar(stat="identity",col="black",position="dodge",size=0.05)+
  theme_classic()+
  my_theme+
  scale_fill_manual(values=c("#377EB8","#E41A1C"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "right",strip.text.x = element_text(size=7,angle=90),axis.text.y = element_text(size=7))+
  facet_grid(~Pair_new,scales="free_x",space="free_x")+
  labs(x="Clone",y="log2 Fold Change (Recipient/ Donor)",fill="")
ggsave(filename = paste0(plots_dir,"clones_logFC_Mono.pdf"),clones_logFC_Mono,device = "pdf",width=6,height = 2.5)

clones_logFC_B_cells<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  pivot_wider(id_cols=c("node","Pair_new","cell_type"),names_from = "individual_type",values_from="cell_frac")%>%
  filter(cell_type=="B_cells")%>%
  filter(Donor>0.01|Recipient>0.01)%>% #Include clones if expanded in the Donor OR the Recipient
  mutate(log2FC=log2(Recipient/Donor))%>%
  mutate(node=paste(node,Pair_new,sep="_"))%>%
  mutate(Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  mutate(bigger_in=ifelse(log2FC>0,"Recipient > Donor","Donor > Recipient"))%>%
  ggplot(aes(y=log2FC,x=forcats::fct_reorder(node,log2FC),fill=bigger_in))+
  geom_bar(stat="identity",col="black",position="dodge",size=0.05)+
  theme_classic()+
  my_theme+
  scale_fill_manual(values=c("#377EB8","#E41A1C"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "right",strip.text.x = element_text(size=7,angle=90),axis.text.y = element_text(size=7))+
  facet_grid(~Pair_new,scales="free_x",space="free_x")+
  labs(x="Clone",y="log2 Fold Change (Recipient/ Donor)",fill="")
ggsave(filename = paste0(plots_dir,"clones_logFC_B_cells.pdf"),clones_logFC_B_cells,device = "pdf",width=6,height = 2.5)

#Look at the SDI in different tissues
#As the total compartment is not captured (and is not the same in each cell type), need to make assumptions about non captured fraction
#Normalizing makes the assumption that the non-captured fraction has equivalent oligoclonality to the captured fraction
#In reality, this is  conservative, as the non-captured fraction is likely to more polyclonal than the captured.
SDI_df<-Map(clones=clones_df,pair=names(clones_df),function(clones,pair) {
  tissueIDs=unique(clones$ID)
  SDI=lapply(tissueIDs,function(tissueID) {
    SDI<-clones%>%filter(ID==tissueID)%>%mutate(plogp=(cell_frac*log(cell_frac)))%>%pull(plogp)%>%sum()%>%prod(-1)
    total_cell_frac=clones%>%filter(ID==tissueID)%>%pull(cell_frac)%>%sum()
    data.frame(tissueID=tissueID,SDI=SDI,total_cell_frac=total_cell_frac,SDI_normalized=SDI/total_cell_frac)
  })%>%bind_rows()%>%mutate(Pair=pair)
  return(SDI)
})%>%bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  #mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))

SDI_by_age_celltype_DoR<-SDI_df%>%
  left_join(Pair_metadata)%>%
  filter(cell_type!="Grans")%>%
  ggplot(aes(x=Age,y=SDI_normalized,col=individual_type))+
  geom_smooth(method="lm",size=0.5,alpha=0.5)+
  geom_point(alpha=0.6,size=0.6)+
  facet_grid(~cell_type)+
  theme_classic()+
  scale_color_manual(values=remove_names(DorR_cols))+
  my_theme+
  theme(legend.position = "bottom")+
  labs(x="Age",y="Shannon Diversity Index\n(Normalized)",col="")

ggsave(filename = paste0(plots_dir,"SDI_by_age_celltype_DoR.pdf"),SDI_by_age_celltype_DoR,device = "pdf",width=4,height = 2.5)

lme.SDI<-lmerTest::lmer(SDI_normalized~Age+cell_type+individual_type+(1|Pair_new),data=SDI_df%>%filter(cell_type!="Grans")%>%left_join(Pair_metadata))
summary(lme.SDI)
confint(lme.SDI)

SDI_normalized_plot<-SDI_df%>%
  ggplot(aes(x=cell_type,y=SDI_normalized,col=individual_type))+
  geom_point(size=0.5)+
  facet_grid(~Pair_new)+
  theme_classic()+
  my_theme+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        legend.key.height = unit(3,"mm"),
        axis.text.x = element_text(angle=90,size=5),
        #strip.background = element_rect(),
        strip.text.x = element_text(margin=unit(c(1,1,1,1),"mm"),size=6))+
  labs(x="Cell type",y="Shannon Diversity Index\n(Normalized)",col="")

ggsave(filename = paste0(plots_dir,"SDI_normalized.pdf"),SDI_normalized_plot,device = "pdf",width=4,height = 2)

SDI_normalized_change<-SDI_df%>%
  pivot_wider(id_cols=c("Pair_new","cell_type"),names_from="individual_type",values_from="SDI_normalized")%>%
  mutate(SDI_difference=Recipient-Donor)%>%
  ggplot(aes(x=Pair_new,y=SDI_difference,fill=cell_type))+
  geom_bar(position="dodge",stat="identity",size=0.25,col="black",width=0.6)+
  geom_hline(yintercept=0)+
  scale_y_continuous(limits=c(-1.9,1.9))+
  #facet_grid(~Pair)+
  theme_classic()+
  theme(legend.text = element_text(size=6),
        legend.title = element_text(size=7),
        legend.key.size = unit(3,"mm"),
        legend.position = "bottom",
        axis.text.x = element_text(angle=90,size=5),
        strip.text.x = element_text(size=7),axis.text.y = element_text(size=5),axis.title = element_text(size=5))+
  labs(x="Pair",y="Change in Shannon Diversity Index\n(Recipient - Donor)",fill="Cell type")

ggsave(filename = paste0(plots_dir,"SDI_normalized_change.pdf"),SDI_normalized_change,device = "pdf",width=2.5,height = 2.5)

##Look at relative shifts of DRIVER MUTATION clonal fractions
driver_FC_df<-Map(res=all_targeted_res,tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=names(all_targeted_res),function(res,tree,post,pair) {
  cat(pair,sep="\n")
  details_targ<-res$details_targ
  driver_mut_df<-details_targ%>%filter(coding_change_chip=="Coding change mutation in driver")
  if(pair=="Pair11") {
    driver_mut_df<-bind_rows(driver_mut_df,Pair11_loss_of_Y_details)
  }
  driver_mut_df<-driver_mut_df%>%
    filter(mut_ref%in%(annotated_drivers_df%>%filter(Decision%in%c("Oncogenic","Possible"))%>%pull(mut_ref))|grepl("LOY",mut_ref))
  
  #Convert each single driver into the overall clone
  driver_mut_df$clone_muts<-sapply(1:nrow(driver_mut_df),function(i) {
    ancestor_nodes<-getAncestors(tree,node=driver_mut_df$node[i],type="all")
    daughter_nodes<-get_all_node_children(node=driver_mut_df$node[i],tree = tree)
    if(any(ancestor_nodes%in%driver_mut_df$node)) {
      ancestral_drivers<-driver_mut_df%>%filter(node%in%ancestor_nodes)%>%pull(variant_ID)
      clone_drivers<-paste0(driver_mut_df$variant_ID[i],"/ ",paste0(ancestral_drivers,collapse="/ "))
    } else {
      clone_drivers<-driver_mut_df$variant_ID[i]
    }
    
    if(any(daughter_nodes%in%driver_mut_df$node)) {
      daughter_drivers<-driver_mut_df%>%filter(node%in%daughter_nodes)%>%pull(variant_ID)
      clone_drivers<-paste0(clone_drivers," (",paste0(daughter_drivers,collapse=", "),")")
    }
    
    return(clone_drivers)
  })
  
  driver_mut_df$clone_gene_muts<-sapply(1:nrow(driver_mut_df),function(i) {
    ancestor_nodes<-getAncestors(tree,node=driver_mut_df$node[i],type="all")
    daughter_nodes<-get_all_node_children(node=driver_mut_df$node[i],tree = tree)
    if(any(ancestor_nodes%in%driver_mut_df$node)) {
      ancestral_drivers<-driver_mut_df%>%filter(node%in%ancestor_nodes)%>%pull(Gene)
      clone_drivers<-paste0(driver_mut_df$Gene[i],"/ ",paste0(ancestral_drivers,collapse="/ "))
    } else {
      clone_drivers<-driver_mut_df$Gene[i]
    }
    
    if(any(daughter_nodes%in%driver_mut_df$node)) {
      daughter_drivers<-driver_mut_df%>%filter(node%in%daughter_nodes)%>%pull(Gene)
      clone_drivers<-paste0(clone_drivers," (",paste0(daughter_drivers,collapse=", "),")")
    }
    
    return(clone_drivers)
  })
  
  driver_mut_df$n_clonal_driver<-sapply(1:nrow(driver_mut_df),function(i) {
    ancestor_nodes<-getAncestors(tree,node=driver_mut_df$node[i],type="all")

    ancestral_drivers<-driver_mut_df%>%filter(node%in%ancestor_nodes)%>%pull(variant_ID)
    n_driver=1+length(ancestral_drivers)
    
    return(n_driver)
  })
  
  driver_mut_df$n_subclonal_driver<-sapply(1:nrow(driver_mut_df),function(i) {
    daughter_nodes<-get_all_node_children(node=driver_mut_df$node[i],tree = tree)
    
    daughter_drivers<-driver_mut_df%>%filter(node%in%daughter_nodes)%>%pull(variant_ID)
    n_subclonal_driver<-length(daughter_drivers)
    
    return(n_subclonal_driver)
  })
  
  cell_types_to_test=bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(cell_type)%>%unique()
  pair_out<-lapply(cell_types_to_test,function(test_cell_type) {
    Donor_tissueID=bulk_smry_all%>%filter(Pair==pair & cell_type==test_cell_type & time_point==0 & individual_type=="Donor")%>%pull(tissueID)%>%unique()
    Recip_tissueID=bulk_smry_all%>%filter(Pair==pair & cell_type==test_cell_type & time_point==0 & individual_type=="Recipient")%>%pull(tissueID)%>%unique()
    
    tissue_out<-lapply(1:nrow(driver_mut_df),function(i) {
      if(driver_mut_df$Mut_type[i]=="SNV"|driver_mut_df$Mut_type[i]=="INDEL") {
        mut_ref<-driver_mut_df$mut_ref[i]
        R_post=post[[Recip_tissueID]]%>%filter(mutation_ID==mut_ref)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
        D_post=post[[Donor_tissueID]]%>%filter(mutation_ID==mut_ref)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
        
      } else if(driver_mut_df$Mut_type[i]=="CNA") {
        mut_ref<-driver_mut_df$mut_ref[i]
        
        #For CNAs, take the clonal fraction of the mid-point mutation of the branch
        node<-driver_mut_df$node[i]
        nmuts_node<-post[[Recip_tissueID]]%>%filter(Node_assignment==node)%>%nrow()
        median_pos<-ceiling(nmuts_node/2)
        
        #Determine the position of the midpoint mutation on recip/ donor branches
        R.which.median<-post[[Recip_tissueID]]%>%filter(Node_assignment==node)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()%>%apply(1,mean)%>%order()%>%.[median_pos]
        D.which.median<-post[[Donor_tissueID]]%>%filter(Node_assignment==node)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()%>%apply(1,mean)%>%order()%>%.[median_pos]
        
        #Determine the actual midpoint mutation
        R.which.mut<-post[[Recip_tissueID]]%>%filter(Node_assignment==node)%>%pull(mutation_ID)%>%.[R.which.median]
        D.which.mut<-post[[Donor_tissueID]]%>%filter(Node_assignment==node)%>%pull(mutation_ID)%>%.[D.which.median]
        
      #Use this to get a posterior of the clonal fractions in recipient and donor
        R_post=post[[Recip_tissueID]]%>%filter(mutation_ID==R.which.mut)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
        D_post=post[[Donor_tissueID]]%>%filter(mutation_ID==D.which.mut)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
      }
      log2FC=log2(R_post/D_post)
      as.data.frame(quantile(log2FC,c(0.025,0.5,0.975)))%>%t()%>%
        as.data.frame()%>%
        dplyr::mutate(Pair=new_pair_names[pair],cell_type=test_cell_type,mut_ref=mut_ref,.before=1)%>%
        dplyr::rename("lowerCI"=`2.5%`,"median"=`50%`,"upperCI"=`97.5%`)%>%
        tibble::remove_rownames()
      
    })%>%dplyr::bind_rows()
    
  })%>%bind_rows()%>%
    left_join(driver_mut_df%>%dplyr::select(mut_ref,node,clone_muts,clone_gene_muts,n_clonal_driver,n_subclonal_driver,variant_ID,Gene))
  return(pair_out)
})%>%dplyr::bind_rows()%>%
  mutate(Pair=factor(Pair,levels=new_pair_names))

max_val= max(abs(driver_FC_df$median))+0.1;min_val=-0.1-max(abs(driver_FC_df$median))
driver_FC_by_pair_plot<-driver_FC_df%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(lowerCI=sapply(lowerCI,function(x) max(x,min_val)),upperCI=sapply(upperCI,function(x) min(x,max_val)))%>%
  mutate(robust_change=ifelse(lowerCI>0,"Recip > Donor",ifelse(upperCI<0,"Recip < Donor","No significant change")))%>%
  mutate(robust_change=factor(robust_change,levels=c("Recip > Donor","Recip < Donor","No significant change")))%>%
  ggplot(aes(x=str_wrap(clone_muts,30),ymin=lowerCI,y=median,ymax=upperCI,col=robust_change))+
  geom_point(size=1)+
  geom_errorbar(size=0.25,width=0)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_color_brewer(palette="Set2")+
  scale_y_continuous(limits=c(min_val,max_val),breaks = seq(floor(min_val),ceiling(max_val),1))+
  labs(x="Cell type",y="Log2 Fold Change (Recip/ Donor)")+
  facet_grid(Pair~cell_type,scales = "free",space="free")+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,size=5),strip.text.x = element_text(size=7),strip.text.y = element_text(size=7,angle=0),axis.text.y = element_text(size=6),axis.title = element_text(size=7),legend.position = "none")

ggsave(filename = paste0(plots_dir,"driver_FC_by_pair_plot.pdf"),driver_FC_by_pair_plot,device = "pdf",width=8,height = 8)

name_vec=str_wrap(driver_FC_df$variant_ID,width=20)
names(name_vec)=driver_FC_df$mut_ref
driver_FC_by_Gene_Monocytes_plot<-driver_FC_df%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(lowerCI=sapply(lowerCI,function(x) max(x,min_val)),upperCI=sapply(upperCI,function(x) min(x,max_val)))%>%
  mutate(mut_cat=ifelse(Gene%in%c("DNMT3A","TET2","TP53","CHEK2","LOY"),Gene,"Other"))%>%
  mutate(mut_cat=factor(mut_cat,levels=c("DNMT3A","TET2","TP53","CHEK2","LOY","Other")))%>%
  mutate(robust_change=ifelse(lowerCI>0,"Recip > Donor",ifelse(upperCI<0,"Recip < Donor","No significant change")))%>%
  mutate(robust_change=factor(robust_change,levels=c("Recip > Donor","Recip < Donor","No significant change")))%>%
  mutate(clone_muts=stringr::str_wrap(paste0(clone_muts," (",Pair,")"),width=30))%>%
  filter(cell_type=="Monocytes" & Gene!="LOY")%>% #Discard granulocytes for this plot as not all patients have results for this
  ggplot(aes(x=mut_ref,ymin=lowerCI,y=median,ymax=upperCI,col=robust_change))+
  geom_point(size=1)+
  geom_errorbar(size=0.25,width=0)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_color_brewer(palette="Set2")+
  #scale_color_manual(values=Pair_cols)+
  scale_x_discrete(labels=name_vec)+
  scale_y_continuous(limits=c(min_val,max_val),breaks = seq(floor(min_val),ceiling(max_val),1))+
  labs(x="Variant",y="Log2 Fold Change (Recip/ Donor)",col="Difference")+
  facet_grid(rows=vars(mut_cat),scales = "free",space="free")+
  coord_flip()+
  theme_bw()+
  my_theme+
  theme(legend.position="none",axis.text.x = element_text(angle=90,size=5),strip.text.x = element_text(size=7),strip.text.y = element_text(size=7),axis.text.y = element_text(size=6),axis.title = element_text(size=7))

ggsave(filename = paste0(plots_dir,"driver_FC_by_Gene_Monocytes_plot.pdf"),driver_FC_by_Gene_Monocytes_plot,device = "pdf",width=3,height = 5.5)

driver_FC_by_Gene_plot<-driver_FC_df%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(lowerCI=sapply(lowerCI,function(x) max(x,min_val)),upperCI=sapply(upperCI,function(x) min(x,max_val)))%>%
  mutate(mut_cat=ifelse(Gene%in%c("DNMT3A","TET2","TP53","CHEK2","LOY"),Gene,"Other"))%>%
  mutate(mut_cat=factor(mut_cat,levels=c("DNMT3A","TET2","TP53","CHEK2","LOY","Other")))%>%
  mutate(robust_change=ifelse(lowerCI>0,"Recip > Donor",ifelse(upperCI<0,"Recip < Donor","No significant change")))%>%
  mutate(robust_change=factor(robust_change,levels=c("Recip > Donor","Recip < Donor","No significant change")))%>%
  mutate(clone_muts=stringr::str_wrap(paste0(clone_muts," (",Pair,")"),width=50))%>%
  filter(cell_type!="Granulocytes")%>% #Discard granulocytes for this plot as not all patients have results for this
  ggplot(aes(x=clone_muts,ymin=lowerCI,y=median,ymax=upperCI,col=robust_change))+
  geom_point(size=1)+
  geom_errorbar(size=0.25,width=0)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_color_brewer(palette="Set2")+
  #scale_color_manual(values=Pair_cols)+
  scale_y_continuous(limits=c(min_val,max_val),breaks = seq(floor(min_val),ceiling(max_val),1))+
  labs(x="",y="Log2 Fold Change (Recip/ Donor)",col="Difference")+
  facet_grid(mut_cat~cell_type,scales = "free",space="free")+
  coord_flip()+
  theme_bw()+
  my_theme+
  theme(axis.title.y = element_blank(),axis.text.x = element_text(angle=90,size=5),strip.text.x = element_text(size=7),strip.text.y = element_text(size=7),axis.text.y = element_text(size=6),axis.title = element_text(size=7))

ggsave(filename = paste0(plots_dir,"driver_FC_by_Gene_plot.pdf"),driver_FC_by_Gene_plot,device = "pdf",width=7,height = 9)

driver_FC_by_ndriver_plot<-driver_FC_df%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(lowerCI=sapply(lowerCI,function(x) max(x,min_val)),upperCI=sapply(upperCI,function(x) min(x,max_val)))%>%
  filter(n_subclonal_driver==0 & cell_type!="Granulocytes")%>%
  mutate(n_clonal_driver=factor(n_clonal_driver,levels=c(1,2,3,4,5)))%>%
  mutate(robust_change=ifelse(lowerCI>0,"Recip > Donor",ifelse(upperCI<0,"Recip < Donor","No significant change")))%>%
  mutate(robust_change=factor(robust_change,levels=c("Recip > Donor","Recip < Donor","No significant change")))%>%
  mutate(clone_muts=stringr::str_wrap(paste0(clone_muts," (",Pair,")"),width=30))%>%
  ggplot(aes(x=clone_muts,ymin=lowerCI,y=median,ymax=upperCI,col=robust_change))+
  geom_point(size=1)+
  geom_errorbar(size=0.25,width=0)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_color_brewer(palette="Set2")+
  #scale_color_manual(values=Pair_cols)+
  scale_y_continuous(limits=c(min_val,max_val),breaks = seq(floor(min_val),ceiling(max_val),1))+
  labs(x="Variant",y="Log2 Fold Change (Recip/ Donor)",col="Difference")+
  facet_grid(n_clonal_driver~cell_type,scales = "free",space="free")+
  coord_flip()+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90,size=5),strip.text.x = element_text(size=7),strip.text.y = element_text(size=7),axis.text.y = element_text(size=6),axis.title = element_text(size=7))

driver_FC_by_ntotdriver_plot<-driver_FC_df%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(lowerCI=sapply(lowerCI,function(x) max(x,min_val)),upperCI=sapply(upperCI,function(x) min(x,max_val)))%>%
  filter(cell_type!="Granulocytes")%>%
  mutate(n_driver=factor(n_clonal_driver+n_subclonal_driver,levels=c(1,2,3,4,5)))%>%
  mutate(robust_change=ifelse(lowerCI>0,"Recip > Donor",ifelse(upperCI<0,"Recip < Donor","No significant change")))%>%
  mutate(robust_change=factor(robust_change,levels=c("Recip > Donor","Recip < Donor","No significant change")))%>%
  mutate(clone_muts=stringr::str_wrap(paste0(clone_muts," (",Pair,")"),width=30))%>%
  ggplot(aes(x=clone_muts,ymin=lowerCI,y=median,ymax=upperCI,col=robust_change))+
  geom_point(size=1)+
  geom_errorbar(size=0.25,width=0)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_color_brewer(palette="Set2")+
  #scale_color_manual(values=Pair_cols)+
  scale_y_continuous(limits=c(min_val,max_val),breaks = seq(floor(min_val),ceiling(max_val),1))+
  labs(x="Driver gene",y="Log2 Fold Change (Recip/ Donor)",col="Difference")+
  facet_grid(n_driver~cell_type,scales = "free",space="free")+
  coord_flip()+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90,size=5),strip.text.x = element_text(size=7),strip.text.y = element_text(size=7),axis.text.y = element_text(size=6),axis.title = element_text(size=7))

name_vec=str_wrap(driver_FC_df$clone_gene_muts,width=20)
names(name_vec)=driver_FC_df$mut_ref
driver_FC_by_ntotdriver_monos_plot<-driver_FC_df%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(lowerCI=sapply(lowerCI,function(x) max(x,min_val)),upperCI=sapply(upperCI,function(x) min(x,max_val)))%>%
  filter(cell_type=="Monocytes")%>%
  mutate(n_driver=factor(n_clonal_driver+n_subclonal_driver,levels=c(1,2,3,4,5)))%>%
  mutate(robust_change=ifelse(lowerCI>0,"Recip > Donor",ifelse(upperCI<0,"Recip < Donor","No significant change")))%>%
  mutate(robust_change=factor(robust_change,levels=c("Recip > Donor","Recip < Donor","No significant change")))%>%
  #mutate(clone_muts=stringr::str_wrap(paste0(clone_muts," (",Pair,")"),width=30))%>%
  arrange(median)%>%
  ggplot(aes(y=forcats::fct_reorder(factor(mut_ref),median),xmin=lowerCI,x=median,xmax=upperCI,col=robust_change))+
  geom_point(size=1)+
  geom_errorbar(size=0.25,width=0)+
  geom_vline(xintercept = 0,linetype=2)+
  scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(min_val,max_val),breaks = seq(floor(min_val),ceiling(max_val),1))+
  scale_y_discrete(labels=name_vec)+
  labs(y="Mutated genes within clone",x="Log2 Fold Change (Recip/ Donor)",col="Difference")+
  facet_grid(rows=vars(n_driver),scales = "free",space="free")+
  theme_bw()+
  my_theme+
  theme(axis.text.y = element_text(angle=0,size=5),strip.text.y = element_text(size=7,angle=0),strip.text.x = element_text(size=7),axis.text.x = element_text(size=6),axis.title = element_text(size=7))

ggsave(filename = paste0(plots_dir,"driver_FC_by_ntotdriver_monos_plot.pdf"),driver_FC_by_ntotdriver_monos_plot,device = "pdf",width=4,height = 5.5)

###Look at MYELOID/ LYMPHOID of driver mutations
# Generate specific dataframe with cell fracs in different compartments
driver_cellfracs_df<-Map(res=all_targeted_res,tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=names(all_targeted_res),function(res,tree,post,pair) {
  cat(pair,sep="\n")
  details_targ<-res$details_targ
  driver_mut_df<-details_targ%>%filter(coding_change_chip=="Coding change mutation in driver")
  if(pair=="Pair11") {
    driver_mut_df<-bind_rows(driver_mut_df,Pair11_loss_of_Y_details)
  }
  driver_mut_df<-driver_mut_df%>%
    filter(mut_ref%in%(annotated_drivers_df%>%filter(Decision%in%c("Oncogenic","Possible"))%>%pull(mut_ref))|grepl("LOY",mut_ref))
  
  tissues_to_test=bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(tissueID)%>%unique()
  pair_out<-lapply(tissues_to_test,function(test_tissueID) {
    cat(test_tissueID,sep="\n")
    tissue_out<-lapply(1:nrow(driver_mut_df),function(i) {
      #cat(i,sep="\n")
      if(driver_mut_df$Mut_type[i]=="SNV"|driver_mut_df$Mut_type[i]=="INDEL") {
        mut_ref<-driver_mut_df$mut_ref[i]
        driver_post=post[[test_tissueID]]%>%filter(mutation_ID==mut_ref)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
        
      } else if(driver_mut_df$Mut_type[i]=="CNA") {
        mut_ref<-driver_mut_df$mut_ref[i]
        
        #For CNAs, take the clonal fraction of the mid-point mutation of the branch
        node<-driver_mut_df$node[i]
        nmuts_node<-post[[test_tissueID]]%>%filter(Node_assignment==node)%>%nrow()
        median_pos<-ceiling(nmuts_node/2)
        
        #Determine the position of the midpoint mutation on recip/ donor branches
        which.median<-post[[test_tissueID]]%>%filter(Node_assignment==node)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()%>%apply(1,mean)%>%order()%>%.[median_pos]
        
        #Determine the actual midpoint mutation
        which.mut<-post[[test_tissueID]]%>%filter(Node_assignment==node)%>%pull(mutation_ID)%>%.[which.median]
        
        #Use this to get a posterior of the clonal fractions in recipient and donor
        driver_post=post[[test_tissueID]]%>%filter(mutation_ID==which.mut)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
      }
      as.data.frame(quantile(driver_post,c(0.025,0.5,0.975)))%>%t()%>%
        as.data.frame()%>%
        dplyr::mutate(Pair=new_pair_names[pair],tissueID=test_tissueID,mut_ref=mut_ref,.before=1)%>%
        dplyr::rename("lowerCI"=`2.5%`,"median"=`50%`,"upperCI"=`97.5%`)%>%
        tibble::remove_rownames()
      
    })%>%dplyr::bind_rows()
    
  })%>%bind_rows()%>%
    left_join(driver_mut_df%>%dplyr::select(mut_ref,Chrom,Pos,Ref,Alt,Mut_type,node,Gene,variant_ID),by="mut_ref")
})

#Generate df with cell fracs for different tissues on the same row
driver_cellfracs_wider_df<-driver_cellfracs_df%>%
  bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(PD_number,tissueID,individual_type,cell_type)%>%filter(!duplicated(.)))%>%
  dplyr::select(Pair,PD_number,individual_type,cell_type,median,lowerCI,upperCI,variant_ID,Gene)%>%
  pivot_wider(names_from=c("cell_type"),values_from=c("lowerCI","median","upperCI"))%>%
  filter(individual_type=="Donor" &
           (median_Monocytes>0.0005|median_Granulocytes>0.0005)&
           Gene%in%c("DNMT3A","TET2"))%>%
  mutate(Pair=factor(Pair,levels=paste("Pair",1:10,sep="_")))

#Generate a df to display confidence intervals as quadrilaterals within geom_polygon
confidence_interval_df<-lapply(1:nrow(driver_cellfracs_wider_df),function(i) {
  
  GM<-driver_cellfracs_wider_df$median_Monocytes[i]
  GL<-driver_cellfracs_wider_df$lowerCI_Monocytes[i]
  GU<-driver_cellfracs_wider_df$upperCI_Monocytes[i]
  
  BM<-driver_cellfracs_wider_df$median_Granulocytes[i]
  BL<-driver_cellfracs_wider_df$lowerCI_Granulocytes[i]
  BU<-driver_cellfracs_wider_df$upperCI_Granulocytes[i]
  
  x=c(BL,BM,BU,BM)
  y=c(GM,GU,GM,GL)
  
  return(data.frame(Pair=driver_cellfracs_wider_df$Pair[i],
                    individual_type=driver_cellfracs_wider_df$individual_type[i],
                    variant_ID=driver_cellfracs_wider_df$variant_ID[i],
                    Gene=driver_cellfracs_wider_df$Gene[i],
                    x=x,y=y))
})%>%bind_rows()%>%
  mutate(id=paste(Pair,individual_type,variant_ID,sep="_"))


driver_cellfracs_wider_df%>%
  ggplot(aes(x=median_Granulocytes,y=median_Monocytes,col=Pair))+
  geom_point(size=0.5)+
  geom_polygon(data=confidence_interval_df,aes(x=x,y=y,group=id,fill=Pair),linewidth=0.5,alpha=0.2)+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline()+
  scale_color_manual(values=Pair_cols)+
  scale_fill_manual(values=Pair_cols)+
  theme_bw()+
  my_theme+
  facet_wrap(~Gene)


#Look at TIMING OF ACQUISITION of driver mutations
driver_fracs_df<-Map(res=all_targeted_res,tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=names(all_targeted_res),function(res,tree,post,pair) {
  details_targ<-res$details_targ
  driver_mut_df<-details_targ%>%filter(coding_change_chip=="Coding change mutation in driver")
  cell_types_to_test=bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(cell_type)%>%unique()
  pair_out<-lapply(cell_types_to_test,function(test_cell_type) {
    Donor_tissueID=bulk_smry_all%>%filter(Pair==pair & cell_type==test_cell_type & time_point==0 & individual_type=="Donor")%>%pull(tissueID)%>%unique()
    Recip_tissueID=bulk_smry_all%>%filter(Pair==pair & cell_type==test_cell_type & time_point==0 & individual_type=="Recipient")%>%pull(tissueID)%>%unique()
    
    tissue_out<-lapply(driver_mut_df$mut_ref,function(mut_ref) {
      R_post=post[[Recip_tissueID]]%>%filter(mutation_ID==mut_ref)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
      D_post=post[[Donor_tissueID]]%>%filter(mutation_ID==mut_ref)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
      data.frame(Recipient=median(R_post),Donor=median(D_post))%>%
        dplyr::mutate(Pair=new_pair_names[pair],cell_type=test_cell_type,mut_ref=mut_ref,.before=1)
    })%>%dplyr::bind_rows()%>%
      gather(-Pair,-cell_type,-mut_ref,key="individual_type",value="cell_frac")
    
  })%>%bind_rows()%>%
    left_join(driver_mut_df%>%dplyr::select(mut_ref,variant_ID,Gene))
  return(pair_out)
})%>%dplyr::bind_rows()%>%
  mutate(Pair=factor(Pair,levels=new_pair_names))

driver_fracs_df%>%
  filter(cell_type=="Monocytes")%>%
  ggplot(aes(x=variant_ID,y=cell_frac,fill=individual_type))+
  #geom_bar(stat="identity",position="dodge")+
  geom_point(aes(col=individual_type),alpha=0.6)+
  facet_grid(rows=vars(Pair),scales="free",space="free")+
  scale_y_log10()+
  scale_color_manual(values = remove_names(DorR_cols))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.y = element_text(angle = 0))+
  coord_flip()+
  labs(x="Mutation",y="Clonal fraction")

driver_fracs_plot<-driver_fracs_df%>%
  filter(cell_type=="Monocytes")%>%
  ggplot(aes(y=variant_ID,x=cell_frac,fill=individual_type))+
  geom_bar(stat="identity",position="dodge")+
  #geom_point(aes(col=individual_type))+
  facet_grid(rows=vars(Pair),scales="free",space="free")+
  facet_wrap(~Pair,scales="free",nrow=2)+
  scale_fill_manual(values=c("#11a0aa","#c82565"))+
  scale_x_continuous(labels=scales::label_comma())+
  #scale_x_log10()+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text = element_text(margin = unit(c(1,0,1,0),"mm")))+
  labs(x="Clonal fraction (Monocytes)",y="Variant",fill="")

ggsave(filename = paste0(plots_dir,"driver_fracs_plot.pdf"),driver_fracs_plot,device = "pdf",width=7,height = 3)


driver_timing_df<-Map(res=all_targeted_res,post=posterior_cell_fracs,pair=names(all_targeted_res),function(res,post,pair) {
  driver_mut_df<-res$details_targ%>%filter(coding_change_chip=="Coding change mutation in driver")
  cell_types_to_test=bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(cell_type)%>%unique()
  pair_out<-lapply(cell_types_to_test,function(test_cell_type) {
    
    tissue_out<-lapply(driver_mut_df$mut_ref,function(Mut_ref) {
      max_cell_type=driver_fracs_df%>%filter(mut_ref==Mut_ref)%>%arrange(desc(cell_frac))%>%pull(cell_type)%>%.[1]
      max_individual_type=driver_fracs_df%>%filter(mut_ref==Mut_ref)%>%arrange(desc(cell_frac))%>%pull(individual_type)%>%.[1]
      
      max_tissueID=bulk_smry_all%>%filter(Pair==pair & cell_type==max_cell_type & time_point==0 & individual_type==max_individual_type)%>%pull(tissueID)%>%unique()
      
      mut_node=post[[max_tissueID]]%>%filter(mutation_ID==Mut_ref)%>%pull(Node_assignment)
      node_muts=post[[max_tissueID]]%>%filter(Node_assignment==mut_node)%>%pull(mutation_ID)
      which.is.driver=which(node_muts==Mut_ref)
      max_post=post[[max_tissueID]]%>%filter(Node_assignment==mut_node)%>%dplyr::select(-Node_assignment,-mutation_ID)%>%as.matrix()
      pos_on_branch=apply(max_post,2,function(x) order(x,decreasing=T)[which.is.driver])/length(node_muts)
      
      return(data.frame(Pair=pair,mut_ref=Mut_ref,pos_on_branch=pos_on_branch))
    })%>%bind_rows()
    
  })%>%bind_rows()%>%
    left_join(driver_mut_df%>%dplyr::select(mut_ref,variant_ID,Gene))
  return(pair_out)
})%>%dplyr::bind_rows()

#This shows that actually the targeted sequencing has little extra information on mutation timing
driver_timing_df%>%
  ggplot(aes(y=variant_ID,x=pos_on_branch,fill=Pair))+
  ggridges::geom_density_ridges(scale=3)

##Look at relative shifts of CLONAL EXPANSION clonal fractions

relative_output_df<-Map(clones=clones_df,post=posterior_cell_fracs,pair=names(all_targeted_res),function(clones,post,pair) {
  cat(pair,sep="\n")
  #Create a matrix of cell fractions across tissues for each clone
  cell_fracs_mat<-clones%>%
    pivot_wider(names_from = "ID",values_from="cell_frac")%>%
    dplyr::select(-node,-Pair)%>%as.matrix()
  
  #Include clones that contribute >2% to any cell compartment in either donor or recipient
  if(pair=="Pair41"|pair=="Pair24") {
    clones_to_include<-apply(cell_fracs_mat,1,function(x) any(x>0.0025))
    lowerlim<-1e-4;upperlim<-1e-2
  } else {
    clones_to_include<-apply(cell_fracs_mat,1,function(x) any(x>0.02))
    lowerlim<-1e-3;upperlim<-1
  }
  
  clone_numbers<-clones%>%pivot_wider(names_from = "ID",values_from="cell_frac")%>%pull(node)%>%.[clones_to_include]
  
  clones_enhanced<-clones%>%
    left_join(bulk_smry_all%>%dplyr::select(tissueID,individual_type,cell_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))
    
  
  node_colors=colfunc(24,c("#a0e3b7","#016876","#88dc40","#702cb4","#64903a","#e771dd","#2cf52b","#d5082d","#21f0b6","#812531","#4bd6fd","#424175","#bcc5eb","#0b5313","#f6932e","#604020","#e2c627","#1f84ec","#de19f7","#a3759e"))
  
  #Mono/ B cell comparison
  Mono_B_comparison_plot<-clones_enhanced%>%
    dplyr::select(-ID,-Pair)%>%
    dplyr::filter(node%in%clone_numbers & (cell_type=="Monocytes"|cell_type=="B_cells"))%>%
    pivot_wider(names_from = "cell_type",values_from="cell_frac")%>%
    mutate(node=factor(gsub("node_","",node)))%>%
    ggplot(aes(y=Monocytes,x=B_cells,col=factor(node),shape=factor(individual_type)))+
    geom_path(aes(group=factor(node)),linetype=1,alpha=0.3)+
    geom_point(size=1.5,alpha=0.5)+
    theme_bw()+
    scale_color_manual(values=node_colors)+
    scale_y_log10(limits=c(lowerlim,upperlim))+
    scale_x_log10(limits=c(lowerlim,upperlim))+
    geom_abline(slope=1,intercept=0,linetype=2)+
    my_theme+
    theme(title=element_text(size=7),legend.key.size = unit(2,"mm"),legend.title = element_text(size=6))+
    labs(col="Clone number",y="Monocyte clonal fraction",x="B cell clonal fraction",shape="Donor or\nRecipient",title = paste(new_pair_names[pair],"Relative lineage output\n(B cells & Monocytes)"))
  
  
  #Mono/ T cell comparison
  Mono_T_comparison_plot<-clones_enhanced%>%
    dplyr::select(-ID,-Pair)%>%
    dplyr::filter(node%in%clone_numbers & (cell_type=="Monocytes"|cell_type=="T_cells"))%>%
    pivot_wider(names_from = "cell_type",values_from="cell_frac")%>%
    mutate(node=factor(gsub("node_","",node)))%>%
    ggplot(aes(y=Monocytes,x=T_cells,col=factor(node),shape=factor(individual_type)))+
    geom_path(aes(group=factor(node)),linetype=1,alpha=0.3)+
    geom_point(size=1.5,alpha=0.5)+
    theme_bw()+
    scale_color_manual(values=node_colors)+
    scale_y_log10(limits=c(lowerlim,upperlim))+
    scale_x_log10(limits=c(lowerlim,upperlim))+
    geom_abline(slope=1,intercept=0,linetype=2)+
    my_theme+
    theme(title=element_text(size=7),legend.key.size = unit(2,"mm"),legend.title = element_text(size=6))+
    labs(col="Clone number",y="Monocyte clonal fraction",x="T cell clonal fraction",shape="Donor or\nRecipient",title = paste(new_pair_names[pair],"Relative lineage output\n(T cells & Monocytes)"))
  
  relative_output_alt_plot<-gridExtra::arrangeGrob(Mono_B_comparison_plot,Mono_T_comparison_plot,nrow = 1)
  ggsave(paste0(plots_dir,"relative_output_by_clone_",new_pair_names[pair],".pdf"),relative_output_alt_plot,width=7,height=3)
  
  return(clones_enhanced)
})%>%dplyr::bind_rows()

#Add the driver information
relative_output_df_enhanced<-relative_output_df%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  mutate(id=paste(Pair,node,sep="_"))%>%
  left_join(driver_FC_df%>%mutate(id=paste(Pair,"node",node,sep="_"))%>%dplyr::select(id,clone_muts,clone_gene_muts),by="id")%>%
  filter(cell_type%in%c("B_cells","T_cells","Monocytes"))%>%
  filter(!duplicated(.))%>%
  pivot_wider(id_cols = c("Pair","node","individual_type","clone_muts","clone_gene_muts"),names_from="cell_type",values_from="cell_frac")

threshold=0.01
relative_output_df_enhanced%>%
  mutate(BM_ratio=B_cells/Monocytes)%>%
  filter(B_cells>threshold|Monocytes>threshold)%>%
  mutate(bias=ifelse(BM_ratio<1,"Myeloid-bias","B-lymphoid bias"))%>%
  dplyr::count(individual_type,bias)%>%
  tidyr::complete(individual_type,bias,fill=list(n=0))

relative_output_df_enhanced%>%
  mutate(BM_ratio=B_cells/Monocytes)%>%
  filter(B_cells>threshold|Monocytes>threshold)%>%
  mutate(bias=ifelse(BM_ratio<1,"Myeloid-bias","B-lymphoid bias"))%>%
  dplyr::count(Pair,individual_type,bias)%>%
  tidyr::complete(Pair,individual_type,bias,fill=list(n=0))%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(x=Pair,y=n,fill=bias))+
  geom_bar(stat="identity",position="dodge")+
  facet_grid(~individual_type)+
  theme_classic()+
  my_theme+
  labs(x="",y="Number of clone",fill="Clone bias")

# 
# B_myeloid_output_comparison_plot<-relative_output_df_enhanced%>%
#   mutate(BM_ratio=B_cells/Monocytes)%>%
#   filter(B_cells>threshold|Monocytes>threshold)%>%
#   mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
#   ggplot(aes(x=B_cells,y=Monocytes,col=Pair))+
#   geom_point(alpha=0.5)+
#   scale_x_log10(limits=c(0.005,1))+
#   scale_y_log10(limits=c(0.005,1))+
#   geom_abline(linetype=2)+
#   facet_grid(~individual_type)+
#   scale_color_manual(values=Pair_cols,drop=T)+
#   theme_bw()+
#   my_theme+
#   labs(x="B cell clonal fraction",y="Monocyte clonal fraction")+
#   theme(legend.key.size = unit(3,"mm"))
# ggsave(paste0(plots_dir,"B_myeloid_output_comparison_plot.pdf"),B_myeloid_output_comparison_plot,width=5,height=1.8)
# 
# T_myeloid_output_comparison_plot<-relative_output_df_enhanced%>%
#   mutate(BM_ratio=T_cells/Monocytes)%>%
#   filter(T_cells>threshold|Monocytes>threshold)%>%
#   mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
#   ggplot(aes(x=T_cells,y=Monocytes,col=Pair))+
#   geom_point(alpha=0.5)+
#   scale_x_log10(limits=c(0.005,1))+
#   scale_y_log10(limits=c(0.005,1))+
#   geom_abline(linetype=2)+
#   facet_grid(~individual_type)+
#   scale_color_manual(values=Pair_cols,drop=T)+
#   theme_bw()+
#   my_theme+
#   labs(x="T cell clonal fraction",y="Monocyte clonal fraction")+
#   theme(legend.key.size = unit(3,"mm"))
# ggsave(paste0(plots_dir,"T_myeloid_output_comparison_plot.pdf"),T_myeloid_output_comparison_plot,width=5,height=1.8)

#Get plots where significance of difference is included by point size
relative_output_B_cells_df<-Map(clones=clones_df,clone_size=clone_sizes,post=posterior_cell_fracs,pair=names(all_targeted_res),function(clones,clone_size,post,pair) {
  cat(pair,sep="\n")
  
  #Include nodes that are >1% in any cell fraction in donor or recipient
  clone_differences<-lapply(c("Donor","Recipient"),function(DorR) {
    clones_enhanced<-clones%>%
      left_join(bulk_smry_all%>%dplyr::select(tissueID,individual_type,cell_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
      mutate(uid=paste(node,Pair,sep="_"))%>%
      filter(cell_type%in%c("Monocytes","B_cells") & individual_type==DorR & uid%in%include_nodes_0.01)%>%
      dplyr::select(-ID)%>%
      pivot_wider(names_from="cell_type",values_from="cell_frac")
    
    if(nrow(clones_enhanced)==0) {stop(return(NULL))}
    
    #Test if B and myeloid fractions in donor are significantly different
    clones_enhanced$signif=sapply(clones_enhanced$node,function(clone_number) {
      B_cell_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="B_cells" & time_point==0)%>%pull(tissueID)%>%.[1]
      Mono_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="Monocytes" & time_point==0)%>%pull(tissueID)%>%.[1]
      
      ratio_post=clone_size[[B_cell_donor_id]][[clone_number]]/clone_size[[Mono_donor_id]][[clone_number]]
      quants=quantile(ratio_post,c(0.025,0.5,0.975))
      signif=ifelse(quants['2.5%']>1|quants['97.5%']<1,T,F)
      return(signif)
    })
    return(clones_enhanced)
  })%>%dplyr::bind_rows()
  
  return(clone_differences)
})%>%dplyr::bind_rows()

B_myeloid_output_comparison_plot<-relative_output_B_cells_df%>%
  mutate(size=ifelse(signif,1,0))%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(x=B_cells,y=Monocytes,col=Pair,size=size))+
  geom_point(alpha=0.5)+
  scale_x_log10(limits=c(0.0009,1))+
  scale_y_log10(limits=c(0.0009,1))+
  scale_size_continuous(limits = c(0,2),range=c(0.75,2),guide = "none")+
  geom_abline(linetype=2)+
  facet_grid(~individual_type)+
  scale_color_manual(values=Pair_cols,drop=T)+
  theme_bw()+
  my_theme+
  labs(x="B cell clonal fraction",y="Monocyte clonal fraction")+
  theme(legend.key.size = unit(3,"mm"))
ggsave(paste0(plots_dir,"B_myeloid_output_comparison_plot.pdf"),B_myeloid_output_comparison_plot,width=5,height=2)

relative_output_T_cells_df<-Map(clones=clones_df,clone_size=clone_sizes,post=posterior_cell_fracs,pair=names(all_targeted_res),function(clones,clone_size,post,pair) {
  cat(pair,sep="\n")
  
  #Include nodes that are >1% in any cell fraction in donor or recipient
  clone_differences<-lapply(c("Donor","Recipient"),function(DorR) {
    clones_enhanced<-clones%>%
      left_join(bulk_smry_all%>%dplyr::select(tissueID,individual_type,cell_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
      mutate(uid=paste(node,Pair,sep="_"))%>%
      filter(cell_type%in%c("Monocytes","T_cells") & individual_type==DorR & uid%in%include_nodes_0.01)%>%
      dplyr::select(-ID)%>%
      pivot_wider(names_from="cell_type",values_from="cell_frac")
    
    if(nrow(clones_enhanced)==0) {stop(return(NULL))}
    
    #Test if B and myeloid fractions in donor are significantly different
    clones_enhanced$signif=sapply(clones_enhanced$node,function(clone_number) {
      T_cell_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="T_cells" & time_point==0)%>%pull(tissueID)%>%.[1]
      Mono_donor_id=bulk_smry_all%>%filter(Pair==pair & individual_type==DorR & cell_type=="Monocytes" & time_point==0)%>%pull(tissueID)%>%.[1]
      
      ratio_post=clone_size[[T_cell_donor_id]][[clone_number]]/clone_size[[Mono_donor_id]][[clone_number]]
      quants=quantile(ratio_post,c(0.025,0.5,0.975))
      signif=ifelse(quants['2.5%']>1|quants['97.5%']<1,T,F)
      return(signif)
    })
    return(clones_enhanced)
  })%>%dplyr::bind_rows()
  
  return(clone_differences)
})%>%dplyr::bind_rows()

T_myeloid_output_comparison_plot<-relative_output_T_cells_df%>%
  mutate(size=ifelse(signif,1,0))%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(x=T_cells,y=Monocytes,col=Pair,size=size))+
  geom_point(alpha=0.5)+
  scale_x_log10(limits=c(0.00075,1))+
  scale_y_log10(limits=c(0.00075,1))+
  scale_size_continuous(limits = c(0,2),range=c(0.75,2),guide = "none")+
  geom_abline(linetype=2)+
  facet_grid(~individual_type)+
  scale_color_manual(values=Pair_cols,drop=T)+
  theme_bw()+
  my_theme+
  labs(x="B cell clonal fraction",y="Monocyte clonal fraction")+
  theme(legend.key.size = unit(3,"mm"))
ggsave(paste0(plots_dir,"T_myeloid_output_comparison_plot.pdf"),T_myeloid_output_comparison_plot,width=5,height=2)

#How many of the 114 clones have significant differences in either B/ T cell comparisons in either donor or recipient
bind_rows(relative_output_B_cells_df,relative_output_T_cells_df)%>%filter(signif)%>%pull(uid)%>%unique()%>%length()

###CAPTURED CELL FRACTIONS THROUGH TIME

get_sum_of_frac=function(tissue_post,tree,cut_off,verbose=F) {
  
  if(verbose) {cat(paste("Clone cutoff of",cut_off,"being used"),sep="\n")}
  
  #Define the cutoff branches
  cutoff_branches=get_cutoff_branches(tree,cut_off = cut_off)
  
  if(verbose) {cat(paste(length(cutoff_branches)," branches at the cutoff level"),sep="\n")}
  
  #If the cutoff branches have no mutations, need to replace with the daughter branches of that node
  number_of_covered_muts=function(node,tissue_post) {sum(tissue_post$Node_assignment==node)} #convenience function
  
  while(any(sapply(cutoff_branches,number_of_covered_muts,tissue_post=tissue_post)==0)){
    updated_nodes<-unlist(lapply(cutoff_branches,function(node) {
      if(number_of_covered_muts(node,tissue_post=tissue_post)==0) {
        if(verbose) {cat(paste("No mutations on node",node,"=> replacing with daughter nodes."),sep="\n")}
        daughter_nodes<-tree$edge[,2][tree$edge[,1]==node]
        return(daughter_nodes)
      } else {
        return(node)
      }
    }))
    cutoff_branches<-updated_nodes
  }
  
  #Define the clone posteriors (as a list)
  clone_posteriors<-lapply(cutoff_branches,function(node) {
    if(verbose) {cat(node,sep="\n")}
    branch_heights=nodeHeights(tree)[tree$edge[,2]==node,]
    min_height=branch_heights[1]
    max_height=branch_heights[2]
    
    if(min_height>cut_off) {
      prop<-0
    } else {
      prop=(cut_off-min_height)/(max_height-min_height)
    }
    
    n_targseq_muts=sum(tissue_post$Node_assignment==node)
    
    which_branch_rank=max(round(prop*n_targseq_muts),1)
    node_cell_frac_distributions=tissue_post%>%
      filter(Node_assignment==node)%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()
    
    #Sort each iteration of posterior by cell fraction & choose that which corresponds to the rank
    if(nrow(node_cell_frac_distributions)>1) {
      #Original method - sorts cell fracs of each Gibbs sample, and takes the 10th rank of each. ?significant overestimates
      #node_cell_frac_distributions_sorted<-apply(node_cell_frac_distributions,2,sort,decreasing=T)
      
      #Second method - keeps values of specific mutations linked together & just orders them by their median values
      rank=order(apply(node_cell_frac_distributions,1,median),decreasing = T)
      node_cell_frac_distributions_sorted<-node_cell_frac_distributions[rank,]
      
    } else {
      node_cell_frac_distributions_sorted<-node_cell_frac_distributions
    }
    
    clone_post=node_cell_frac_distributions_sorted[which_branch_rank,]
    return(clone_post)
  })
  names(clone_posteriors)<-paste("node",cutoff_branches,sep="_")
  
  #Use these to get the sum of VAF
  n_iter=length(clone_posteriors[[1]])
  sum_of_frac_post<-sapply(1:100,function(j) {
    sum(sapply(clone_posteriors,function(x) x[j]))
  })
  return(sum_of_frac_post)
}

sum_of_frac_list_file=paste0(root_dir,"/data/Targeted_sequencing_data/sum_of_fracs.Rds")

if(file.exists(sum_of_frac_list_file)){
  sum_of_frac_list<-readRDS(sum_of_frac_list_file)
} else {
  #This is slow
  sum_of_frac_list<-Map(tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=names(all_targeted_res),function(tree,post,pair) {
    cat(pair,sep="\n")
    tree.ultra<-make.ultrametric.tree(tree)
    tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(drop.tip(tree,"Ancestral")))
    tree.ultra$edge.length[tree.ultra$edge[,2]==which(tree.ultra$tip.label=="Ancestral")]<-0
    
    tissueIDs=names(post)
    pair_sum_of_frac_by_cutoff<-lapply(tissueIDs,function(tissueID){
      cat(tissueID,sep="\n")
      tissue_post<-post[[tissueID]]
      tissue_sum_of_frac_by_cutoff<-lapply(seq(2,100,2),function(clone_cutoff) {
        cat(clone_cutoff,sep="\n")
        get_sum_of_frac(tissue_post = tissue_post,cut_off = clone_cutoff,tree = tree.ultra,verbose=F)
      })
      return(tissue_sum_of_frac_by_cutoff)
    })
    names(pair_sum_of_frac_by_cutoff)<-tissueIDs
    return(pair_sum_of_frac_by_cutoff)
  })
  saveRDS(sum_of_frac_list,file=sum_of_frac_list_file)
}

##Extract the median and 95% CI of the sum of frac statistic
all_sof<-Map(pair=names(all_targeted_res),pair_sum_of_frac_by_cutoff=sum_of_frac_list,f=function(pair,pair_sum_of_frac_by_cutoff) {
  pair_cell_fracs<-Map(tissue=pair_sum_of_frac_by_cutoff,tissueID=names(pair_sum_of_frac_by_cutoff),f=function(tissue,tissueID) {
    cat(tissueID,sep="\n")
    tissue_df<-Map(cutoff=seq(2,100,2),cutoff_sof=tissue,function(cutoff,cutoff_sof) {
      data.frame(cutoff=cutoff,lowerCI=quantile(cutoff_sof,0.025),median=median(cutoff_sof),upperCI=quantile(cutoff_sof,0.975))
    })%>%dplyr::bind_rows()%>%
      mutate(tissueID=tissueID,.before=1)
    return(tissue_df)
  })%>%dplyr::bind_rows()%>%
    mutate(PairID=pair,.before=1)
  return(pair_cell_fracs)
})


##Plot the results
SOF_plot<-all_sof%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,individual_type,cell_type),by="tissueID")%>%
  dplyr::mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")),
                individual_type=factor(individual_type))%>%
  dplyr::mutate(Pair_new=factor(Pair_new,levels=paste("Pair",1:10,sep="_")))%>%
  ggplot(aes(ymin=lowerCI,ymax=upperCI,y=median,x=cutoff,col=individual_type,fill=individual_type))+
  geom_hline(yintercept = 1,linetype=2)+
  geom_ribbon(alpha=0.4,col=NA)+
  geom_line(size=0.3)+
  facet_grid(Pair_new~cell_type)+
  scale_y_continuous(breaks=seq(0,1,0.5),limits=c(0,1))+
  scale_fill_manual(values=remove_names(DorR_cols))+
  scale_color_manual(values=remove_names(DorR_cols))+
  theme_bw()+
  my_theme+
  labs(x="Molecular time",y="Captured cell fraction")

ggsave(filename=paste0(plots_dir,"Sum_of_frac_through_time.pdf"),SOF_plot,width=6,height=8)

SOF_plot_Tcell_only<-all_sof%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,individual_type,cell_type),by="tissueID")%>%
  dplyr::mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")),
                individual_type=factor(individual_type))%>%
  dplyr::mutate(Pair_new=factor(Pair_new,levels=paste("Pair",1:10,sep="_")))%>%
  dplyr::filter(cell_type=="T_cells")%>%
  ggplot(aes(ymin=lowerCI,ymax=upperCI,y=median,x=cutoff,col=individual_type,fill=individual_type))+
  geom_hline(yintercept = 1,linetype=2)+
  geom_ribbon(alpha=0.4,col=NA)+
  geom_line(size=0.3)+
  facet_grid(cell_type~Pair_new)+
  scale_y_continuous(breaks=seq(0,1,0.5),limits=c(0,1))+
  scale_fill_manual(values=remove_names(DorR_cols))+
  scale_color_manual(values=remove_names(DorR_cols))+
  theme_bw()+
  my_theme+
  labs(x="Molecular time",y="Captured cell fraction",fill="Donor or\nRecipient",col="Donor or\nRecipient")

ggsave(filename=paste0(plots_dir,"Sum_of_frac_through_time_Tcell_only.pdf"),SOF_plot_Tcell_only,width=7,height=2)


all_sof%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,individual_type,cell_type),by="tissueID")%>%
  dplyr::mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")),
                individual_type=factor(individual_type))%>%
  dplyr::mutate(Pair_new=factor(Pair_new,levels=paste("Pair",1:10,sep="_")))%>%
  filter(cutoff==100)%>%
  ggplot(aes(x=individual_type,y=median,fill=cell_type))+
  geom_bar(stat="identity",position="dodge",col="black")+
  facet_grid(~Pair_new)+
  theme_bw()+
  my_theme+
  labs(fill="Cell type",x="Donor or Recipient",y="Captured cell fraction\n(at 100 mutations of molecular time)")
