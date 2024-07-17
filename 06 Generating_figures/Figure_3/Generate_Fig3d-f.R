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

#Read in the versions of the tree/ details object used for the Bayesian model - these are subtly different
all_pairs<-Pair_metadata$Pair
tree_folder=paste0(root_dir,"/data/tree_and_mutation_files/trees_no_dups")
annotated_muts_folder=paste0(root_dir,"/data/tree_and_mutation_files/annot_files_no_dups")
tree_paths=list.files(tree_folder,pattern="vaf_post_mix_post_dup.tree",full.names = T)
all.trees.cc.nodups<-lapply(all_pairs,function(PairID){read.tree(grep(PairID,tree_paths,value = T))})
names(all.trees.cc.nodups)<-all_pairs

annotated_muts_paths=list.files(annotated_muts_folder,pattern="vaf_post_mix_post_dup",full.names = T)
all.muts.nodups<-lapply(all_pairs,function(PairID){cat(PairID);load(grep(PairID,annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
names(all.muts.nodups)<-all_pairs

## Generate information regarding loss-of-Y in male samples from X and Y coverage data ----
LOY_files=list.files(path=paste0(root_dir,"/data/SV_and_CNA_data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Create dataframe of the LOY events
Pair11_LOY_nodes=c(48, 373,409 ,450 ,156 ,155 ,493 ,541,626)
Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                    Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                    Protein="",Type="",SO_codes="",coding_change="Coding change",
                                    coding_change_chip="yes",
                                    ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))


## Read in the spreadsheet listing other copy number changes ----
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

#========================================#
# Import the targeted sequencing metadata ####
#========================================#

## First read in the metadata for all the targeted sequencing samples
bulk_smry_all<-readr::read_csv(file = paste0(root_dir,"/data/targeted_sequencing_data/targeted_sequencing_metadata.csv"))

#Read in other data objects
annotated_drivers_df<-read_csv(paste0(root_dir,"/data/tree_and_mutation_files/Possible_drivers_annotated.csv"))[,c("mut_ref","node","Gene","variant_ID","Decision")]
remove_names=function(x) {attr(x,"names")<-NULL;return(x)}

#========================================#
# Import the raw targeted sequencing count data ####
#========================================#

#Import the raw targeted sequencing results - count data across all samples/ mutations sites
#This is generated using alleleCounter across all the targeted sequencing bams
all_targeted_res<-readRDS(paste0(root_dir,"/data/targeted_sequencing_data/all_targeted_res.Rds"))

#========================================#
# GET CELL FRACTIONS FROM THE MODEL ####
#========================================#

## Import the output of the phylogeny-aware Gibb's sampler that infers the posterior distribution of the cell fraction of each mutation ----
## For scripts to run the model using the count data and the phylogenies, please see separate directory 'Gibbs Sampler for Targ Seq'

posterior_cell_fracs_file=paste0(root_dir,"/data/targeted_sequencing_data/posterior_cell_fracs.Rds")

if(file.exists(posterior_cell_fracs_file)) {
  cat("Reading in previously saved RDS file of the posterior cell fractions",sep="\n")
  posterior_cell_fracs<-readRDS(posterior_cell_fracs_file)  ## THIS FILE NEEDS TO BE DOWNLOADED FROM MENDELEY DATA
} else {
  cat("Importing posterior VAFs files",sep="\n")
  setwd(paste0(root_dir,"/data/targeted_sequencing_data/raw_model_output"))  ## ALTERNATIVELY CAN DOWNLOAD THE RAW MODEL OUTPUT WHICH IS THEN PROCESSED IN THESE FUNCTIONS
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

#========================================#
# GET CLONAL CONTRIBUTIONS OF CLONAL EXPANSIONS TO THE DIFFERENT COMPARTMENTS ####
#========================================#
# How to define a 'clone' when looking at the phylogenetic tree? To avoid 'double-counting' you have to define a single time point
# when you consider the earliest time point for clonal expansions to begin, and any expansion beyond that point you consider as 'subclones'.
# If you define this point too early, you will start including normal expansion within development.
# We define it here as 100 mutations of molecular time. The developmental 'burst' of branch points is over by then, but it will capture the root of most abnormal expansions.

#Convenience function to work out the branches that traverse the defined cut off point.
get_cutoff_branches=function(tree,cut_off) {
  heights=nodeHeights(tree)
  cutoff_branches=tree$edge[,2][heights[,1]<cut_off & heights[,2]>=cut_off]
  return(cutoff_branches)
}

#Get the posterior distributions for each 'clone' in each cell type defined by cutting tree at 100 mutations of molecular time
clone_cutoff=100
clone_sizes=Map(tree=all.trees.cc.nodups,post=posterior_cell_fracs,pair=all_pairs,function(tree,post,pair) {
  cat(pair,sep="\n")
  cutoff_branches=get_cutoff_branches(tree,cut_off = clone_cutoff)
  tissueIDs=names(post)
  tissue_clone_fractions<-lapply(tissueIDs,function(tissueID) {
    cat(tissueID,sep="\n")
    tissue_post=post[[tissueID]]
    clone_posteriors<-lapply(cutoff_branches,function(node) {
      
      #For each branch, need to work out which mutation most closely corresponds corresponds to the '100 mutations of molecular time' cutoff
      
      #1. Find the molecular time of the ends of the branch that crosses the time point
      branch_heights=nodeHeights(tree)[tree$edge[,2]==node,]
      min_height=branch_heights[1]
      max_height=branch_heights[2]
      
      #2. From this, work out how far down the branch the '100 mutations of molecular time' is
      prop=(clone_cutoff-min_height)/(max_height-min_height)
      
      #3. Find out how many mutations from this branch are covered by the targeted sequencing.
      # Assuming that the order of these mutations along the branch can be deduced by ordering them by decreasing VAF (i.e. their rank), work out which rank most closely
      # corresponds to the branch position corresponding to '100 mutations o molecular time'
      n_targseq_muts=sum(tissue_post$Node_assignment==node)
      which_branch_rank=max(round(prop*n_targseq_muts),1)
      
      #4. Pull out the posterior distributions of mutations on this branch only, and convert into a matrix
      node_cell_frac_distributions=tissue_post%>%
        filter(Node_assignment==node)%>%
        dplyr::select(-Node_assignment,-mutation_ID)%>%
        as.matrix()
      
      #5. Sort the posterior cell fractions by their median values (to approximate the order) & select the posterior which corresponds to the rank defined above
      if(nrow(node_cell_frac_distributions)>1) {
        
        #Keeps values of specific mutations linked together & just orders them by their median values
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

#========================================#
# PLOT OUT THE EXPANDED CLONES ####
#========================================#

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

## Plot Fig 3d - the clonal contributions in the different cell types ----
clones_100_over1inDorR_plot<-clones_df%>%
  dplyr::bind_rows()%>%
  left_join(bulk_smry_all%>%dplyr::select(Pair_new,tissueID,cell_type,individual_type)%>%filter(!duplicated(.)),by=c("ID"="tissueID"))%>%
  mutate(cell_type=ifelse(cell_type=="Granulocytes","Grans",ifelse(cell_type=="Monocytes","Monos",cell_type)))%>%
  mutate(cell_type=factor(cell_type,levels=c("Grans","Monos","B_cells","T_cells")))%>%
  mutate(individual_type=substr(individual_type,1,5))%>%
  mutate(uid=paste(node,Pair,sep="_"),Pair_new=factor(Pair_new,levels = paste("Pair",1:10,sep="_")))%>%
  filter(uid%in%include_nodes_0.01)%>%
  ggplot(aes(x=factor(cell_type),y=cell_frac,fill=node,))+
  geom_bar(lwd=0.15,stat="identity",position="stack")+
  theme_classic()+
  my_theme+
  scale_fill_manual(values=clone_cols)+
  scale_color_manual(values=c("black","red"))+
  guides(fill="none",linewidth="none")+
  theme(axis.text.x = element_text(angle=90),panel.spacing = unit(0.5, "mm"),strip.text.x = element_text(size=5,margin = unit(c(0.7,1,0.7,1),"mm")))+
  facet_wrap(Pair_new~individual_type,nrow=1)+
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  labs(x="Cell type",y="Clonal fraction")
ggsave(filename = paste0(plots_dir,"Fig3d.pdf"),clones_100_over1inDorR_plot,device = "pdf",width=7,height = 2.5)

#========================================#
# SHANNON DIVERSITY ####
#========================================#

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

## Plot Fig 3e - the change in SDI between donor and recipient ----
SDI_normalized_change<-SDI_df%>%
  pivot_wider(id_cols=c("Pair_new","cell_type"),names_from="individual_type",values_from="SDI_normalized")%>%
  mutate(SDI_difference=Recipient-Donor)%>%
  ggplot(aes(x=Pair_new,y=SDI_difference,fill=cell_type))+
  geom_bar(position="dodge",stat="identity",size=0.25,col="black",width=0.6)+
  geom_hline(yintercept=0)+
  scale_y_continuous(limits=c(-1.9,1.9))+
  guides(fill=guide_legend(nrow=2,ncol=2,position="inside"))+
  theme_classic()+
  my_theme+
  theme(legend.key.size = unit(3,"mm"),
        axis.text.x=element_text(angle = 90),
        legend.position.inside = c(0.4,0.95))+
  labs(x="",y="Change in Shannon Diversity Index\n(Recipient - Donor)",fill="")

ggsave(filename = paste0(plots_dir,"Fig3e.pdf"),SDI_normalized_change,device = "pdf",width=2.5,height = 2.5)

## Plot Fig 3f - the change in SDI with age and by cell compartment ----
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
  labs(x="Age (years)",y="Shannon Diversity Index\n(Normalized)",col="")

#Replace the strip colours with different ones corresponding to the different cell types
g <- ggplot_gtable(ggplot_build(SDI_by_age_celltype_DoR))
fills <- c("#03AC13","#30D5C8","#B660CD")
strip_both <- which(grepl('strip-', g$layout$name))
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
SDI_by_age_celltype_DoR<-g

ggsave(filename = paste0(plots_dir,"Fig3f.pdf"),SDI_by_age_celltype_DoR,device = "pdf",width=4,height = 2.5)

#Do the linear mixed effects regression of Age/ Cell type/ Donor or Recipient----
lme.SDI<-lmerTest::lmer(SDI_normalized~Age+cell_type+individual_type+(1|Pair_new),data=SDI_df%>%filter(cell_type!="Grans")%>%left_join(Pair_metadata))
summary(lme.SDI)
confint(lme.SDI)
