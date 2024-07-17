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

## Read in the spreadsheet listing other copy number changes ----
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

#Create dataframe of the LOY events
Pair11_LOY_nodes=c(48, 373,409 ,450 ,156 ,155 ,493 ,541,626)
Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                    Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                    Protein="",Type="",SO_codes="",coding_change="Coding change",
                                    coding_change_chip="yes",
                                    ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))


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
  setwd(paste0(root_dir,"/data/Targeted_sequencing_data/raw_model_output"))  ## ALTERNATIVELY CAN DOWNLOAD THE RAW MODEL OUTPUT WHICH IS THEN PROCESSED IN THESE FUNCTIONS
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
# EXPLORE DRIVER MUTATIONS - relative shifts in different clonal fractions ####
#========================================#
## Generate a driver 'fold change' dataframe, including all the changes in driver cell fraction
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

## Generate Fig. 5a ----
max_val= max(abs(driver_FC_df$median))+0.1;min_val=-0.1-max(abs(driver_FC_df$median))
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
  geom_point(size=0.5)+
  geom_errorbar(size=0.25,width=0)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_color_manual(values=c("#CE4D99","#4CBB17","#8DA0CB"))+
  scale_x_discrete(labels=name_vec)+
  scale_y_continuous(limits=c(min_val,max_val),breaks = seq(round(min_val),round(max_val),2))+
  labs(x="",y=bquote(~Log[2]~Fold~Change~(Recip/ Donor)),col="Difference")+
  facet_grid(rows=vars(mut_cat),scales = "free",space="free")+
  coord_flip()+
  theme_bw()+
  my_theme+
  theme(legend.position="none",strip.text.y = element_text(size=7,face="bold.italic"),axis.text.y = element_text(size=6,face="italic"),axis.title = element_text(size=7))

ggsave(filename = paste0(plots_dir,"Fig5a.pdf"),driver_FC_by_Gene_Monocytes_plot,device = "pdf",width=3,height = 5.5)

## Generate Fig. 5b ----
#This plot divides clones up by how many drivers are contained within them, either clonal or subclonal
name_vec=str_wrap(driver_FC_df$clone_gene_muts,width=20)
names(name_vec)=driver_FC_df$mut_ref
cat_names=c("Clone with 1 driver mutation","2 drivers","3","5")
names(cat_names)=c(1,2,3,5)

driver_FC_by_ntotdriver_monos_plot<-driver_FC_df%>%
  mutate(cell_type=factor(cell_type,levels=c("Granulocytes","Monocytes","B_cells","T_cells")))%>%
  mutate(lowerCI=sapply(lowerCI,function(x) max(x,min_val)),upperCI=sapply(upperCI,function(x) min(x,max_val)))%>%
  filter(cell_type=="Monocytes")%>%
  mutate(n_driver=factor(n_clonal_driver+n_subclonal_driver,levels=c(1,2,3,4,5)))%>%
  mutate(robust_change=ifelse(lowerCI>0,"Recip > Donor",ifelse(upperCI<0,"Recip < Donor","No significant change")))%>%
  mutate(robust_change=factor(robust_change,levels=c("Recip > Donor","Recip < Donor","No significant change")))%>%
  mutate(n_driver=factor(cat_names[as.character(n_driver)],levels=cat_names))%>%
  arrange(median)%>%
  ggplot(aes(y=forcats::fct_reorder(factor(mut_ref),median),xmin=lowerCI,x=median,xmax=upperCI,col=robust_change))+
  geom_point(size=1)+
  geom_errorbar(size=0.25,width=0)+
  geom_vline(xintercept = 0,linetype=2)+
  scale_color_manual(values=c("#CE4D99","#4CBB17","#8DA0CB"))+
  scale_x_continuous(limits=c(min_val,max_val),breaks = seq(floor(min_val),ceiling(max_val),1))+
  scale_y_discrete(labels=name_vec)+
  labs(y="",x=bquote(~Log[2]~Fold~Change~(Recip/ Donor)),col="")+
  facet_grid(rows=vars(n_driver),scales = "free",space="free")+
  theme_bw()+
  my_theme+
  theme(legend.position="bottom",axis.text.y = element_text(angle=0,size=5,face="italic"),strip.text.y = element_text(size=7,angle=270,face="bold"),strip.text.x = element_text(size=7),axis.text.x = element_text(size=6),axis.title = element_text(size=7))

ggsave(filename = paste0(plots_dir,"Fig5b.pdf"),driver_FC_by_ntotdriver_monos_plot,device = "pdf",width=4,height = 5.5)
