#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","tidyr","readr","phytools","data.table","pheatmap")

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
                strip.background = element_rect(fill="white",linewidth = 0.4),
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
tree_folder=paste0(root_dir,"/data/trees_no_dups")
annotated_muts_folder=paste0(root_dir,"/data/annot_files_no_dups")
tree_paths=list.files(tree_folder,pattern="vaf_post_mix_post_dup.tree",full.names = T)
all.trees.cc.nodups<-lapply(all_pairs,function(PairID){read.tree(grep(PairID,tree_paths,value = T))})
names(all.trees.cc.nodups)<-all_pairs

annotated_muts_paths=list.files(annotated_muts_folder,pattern="vaf_post_mix_post_dup",full.names = T)
all.muts.nodups<-lapply(all_pairs,function(PairID){cat(PairID);load(grep(PairID,annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
names(all.muts.nodups)<-all_pairs

## Generate information regarding loss-of-Y in male samples from X and Y coverage data ----
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
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
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#========================================#
# Import the targeted sequencing metadata ####
#========================================#

## First read in the metadata for all the targeted sequencing samples
bulk_smry_all<-readr::read_csv(file = paste0(root_dir,"/data/targeted_sequencing_data/targeted_sequencing_metadata.csv"))

#Read in other data objects
annotated_drivers_df<-read_csv(paste0(root_dir,"/data/Possible_drivers_annotated.csv"))[,c("mut_ref","node","Gene","variant_ID","Decision")]
remove_names=function(x) {attr(x,"names")<-NULL;return(x)}

#========================================#
# Import the raw targeted sequencing count data ####
#========================================#

#Import the raw targeted sequencing results - count data across all samples/ mutations sites
#This is generated using alleleCounter across all the targeted sequencing bams
all_targeted_res<-readRDS(paste0(root_dir,"/data/Targeted_sequencing_data/all_targeted_res.Rds"))

#========================================#
# GET CELL FRACTIONS FROM THE MODEL ####
#========================================#

## Import the output of the phylogeny-aware Gibb's sampler that infers the posterior distribution of the cell fraction of each mutation ----
## For scripts to run the model using the count data and the phylogenies, please see separate directory 'Gibbs Sampler for Targ Seq'

posterior_cell_fracs_file=paste0(root_dir,"/data/Targeted_sequencing_data/posterior_cell_fracs.Rds")

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
# ASSESS SIMILARITY OF CLONAL COMPOSITION BETWEEN TISSUES ####
#========================================#
##One way of assessing similarity of clonal composition between tissues is using the
##soft cosine similarity of the targeted sequencing data between tissues/ individuals

##Define custom functions for this analysis
generate_median_fracs_across_tissues=function(post_ind,selected_nodes) {
  tissues=names(post_ind)
  mutation_order=post_ind[[1]]%>%filter(Node_assignment%in%selected_nodes)%>%pull(mutation_ID)
  temp<-lapply(post_ind,function(post) {
    post_emb<-post%>%filter(Node_assignment%in%selected_nodes)
    cell_frac_distributions=post_emb%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()
    median_cell_fracs<-apply(cell_frac_distributions,1,median)
    names(median_cell_fracs)<-post_emb$mutation_ID
    return(median_cell_fracs)
  })
  all_tissue_median_fracs<-Reduce(cbind,lapply(temp,function(x) x[mutation_order]))
  dimnames(all_tissue_median_fracs)[[2]]<-tissues
  return(all_tissue_median_fracs)
}

##Function to characterise the shortest phylogenetic distance between mutations
#Uses a reference object for distances between nodes
shortest_dist = function(mut_a,mut_b, tree, ref,node_comb_dist) { #the function to look up the distance
  node_a=ref$Node_assignment[ref$mutation_ID==mut_a]
  node_b=ref$Node_assignment[ref$mutation_ID==mut_b]
  genetic_dist<-node_comb_dist[[paste(node_a,node_b,sep="-")]]
  return(genetic_dist)
}

##Define functions to calculate the cosine similarity for tissue pairs
sum_tissue_scores = function(data,selected_muts,tissue_1,tissue_2,mut_pairs_sim) {
  tissue_pairs_df=as.data.frame(data.table::CJ(a=data[tissue_1,selected_muts],b=data[tissue_2,selected_muts],unique=FALSE,sorted = FALSE))
  tissue_pairs_df$sim=mut_pairs_sim
  tissue_pairs_df$numerator=apply(tissue_pairs_df[,1:3],1,prod)
  return(sum(tissue_pairs_df$numerator)/2)
}

get_soft_cosim=function(data,selected_muts,tissue_1,tissue_2,mut_pairs_sim) {
  num=sum_tissue_scores(data=data,selected_muts=selected_muts,tissue_1,tissue_2,mut_pairs_sim = mut_pairs_sim)
  denom=sum_tissue_scores(data=data,selected_muts=selected_muts,tissue_1,tissue_1,mut_pairs_sim = mut_pairs_sim)^0.5 * sum_tissue_scores(data=data,selected_muts=selected_muts,tissue_2,tissue_2,mut_pairs_sim = mut_pairs_sim)^0.5
  soft_cosim=num/denom
  return(soft_cosim)
}

## Generate the similarity matrices -----
## For this analysis we focus on the early embryonic mutations (<30 mutations of molecular time, corresponds to ~1st trimester)
# Most of the cell compartments will have most information on these early mutations
# The analysis also becomes very costly if start to have too many mutations (>1,000) - therefore if there are to many, need to reduce
# Could do this by random subsampling, but might miss key early branches, therefore reduce by gradually reducing the cut-off time until the number of mutations is <1,000

#This step takes some time so we have pre-saved the output to this
similarity_matrices_file=paste0(root_dir,"/data/similarity_matrices.Rds")
if(!file.exists(similarity_matrices_file)) {
  similarity_matrices=list()
  for(j in 1:length(posterior_cell_fracs)) {
    embryonic_mutation_cutoff=30
    tree=all.trees.cc.nodups[[j]]
    pair=names(all.trees.cc.nodups)[j]
    post_ind=posterior_cell_fracs[[j]]
    cat(pair,sep="\n")
    
    #Get vector of the tissue IDs that relate to the initial time point only
    initial_time_point_samples<-bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(tissueID)%>%unique()
    if(!all(names(post_ind)%in%initial_time_point_samples)) {
      post_ind<-post_ind[-which(!names(post_ind)%in%initial_time_point_samples)]
    }
    
    ##Select which branches are definitively 'embryonic' - i.e. all mutations are from early life
    embryonic_branches=tree$edge[,2][which(nodeHeights(tree)[,2]<embryonic_mutation_cutoff)]
    
    #Rename the tissue IDs with their donor/ recipient and cell type identity
    embryonic_median_cell_fracs<-generate_median_fracs_across_tissues(post_ind,selected_nodes=embryonic_branches)
    tissues=sapply(colnames(embryonic_median_cell_fracs),function(id) {
      bulk_smry_all%>%filter(tissueID==id)%>%mutate(new_id=paste(individual_type,cell_type,sep="_"))%>%pull(new_id)%>%.[1]
    })
    colnames(embryonic_median_cell_fracs)=tissues
    
    #(2) Get the names of the mutations that are on these embryonic branches
    embryonic_muts=rownames(embryonic_median_cell_fracs) #get a vector of mut_refs for these differing mutations
    nmut_max=1000
    cat(paste("There are",length(embryonic_muts),"embryonic mutations"),sep="\n")
    while(length(embryonic_muts)>nmut_max) {
      cat(paste("Reducing mutation set by reducing cut-off molecular time to",embryonic_mutation_cutoff-1),sep="\n")
      embryonic_mutation_cutoff<-embryonic_mutation_cutoff-1
      embryonic_branches=tree$edge[,2][which(nodeHeights(tree)[,2]<embryonic_mutation_cutoff)]
      embryonic_median_cell_fracs<-generate_median_fracs_across_tissues(post_ind,selected_nodes=embryonic_branches)
      embryonic_muts=rownames(embryonic_median_cell_fracs) #get a vector of mut_refs for these differing mutations
    }
    cat(paste("There are",length(embryonic_muts),"embryonic mutations included"),sep="\n")
    colnames(embryonic_median_cell_fracs)=tissues
    
    #(3) Calculate the phylogenetic distances for these mutations
    mut_pairs=as.data.frame(data.table::CJ(a=embryonic_muts,b=embryonic_muts,unique=FALSE,sorted=FALSE)) #create a df of all possible pairs
    mut_pairs$uid=apply(mut_pairs[,1:2],1,paste,collapse="-") #create a unique id for each pair
    cat(paste("This results in",nrow(mut_pairs),"mutation pairs"),sep="\n") #How many are included - gives an idea of how long it will take
    
    node_combs=as.data.frame(data.table::CJ(a=embryonic_branches,b=embryonic_branches,unique=TRUE))
    
    #Calculate genetic distances for all of these node pairs - save as a list
    cat("Calculating genetic distances between node pairs",sep="\n")
    nodeheights=nodeHeights(tree) #The lapply function needs the nodeheights object
    node_comb_dist=lapply(1:nrow(node_combs),function(i) {
      node_a=node_combs[i,1]
      node_b=node_combs[i,2]
      node_height_a=nodeheights[tree$edge[,2]==node_a,2]
      node_height_b=nodeheights[tree$edge[,2]==node_b,2]
      if(node_a==node_b) {stop(return(0))}
      ancestral_nodes_a=get_ancestral_nodes(node_a,tree$edge) #get all the ancestral nodes of each node & include the node itself
      ancestral_nodes_b=get_ancestral_nodes(node_b,tree$edge)
      common_ancestors=dplyr::intersect(c(node_a,ancestral_nodes_a),c(node_b,ancestral_nodes_b)) #pull out those that are common ancestors
      if(length(common_ancestors)==0) {
        genetic_dist <-(node_height_a+node_height_b) #If there are no common ancestors listed, then the closest common ancestor is the tree root & the genetic distance is the sum of the node heights
      } else {
        common_ancestors_heights=sapply(common_ancestors, function(x) {
          nodeheights[tree$edge[,2]==x,2] #otherwise find the heights of all common ancestors
        })
        genetic_dist = node_height_a+node_height_b-2*max(common_ancestors_heights) #the common ancestor with the maximum nodeheight is the most recent
      }
      return(genetic_dist)
    })
    
    #Name the list with the node pair in the format "node_1-node_2"
    node_comb_names=apply(node_combs,1,paste,collapse="-")
    names(node_comb_dist) <- node_comb_names
    
    #Now using this list as a reference, get the distance for all individual mutation pairs
    cat("Calculating the distance between mutation pairs",sep="\n")
    mut_pairs$dist<-NA
    #mut_pairs$dist[mut_pairs$uid %in% rownames(ref_dist)]=ref_dist[mut_pairs$uid[mut_pairs$uid %in% rownames(ref_dist)],"dist"]
    empties=which(is.na(mut_pairs$dist))
    for(i in empties) {
      mut_pairs$dist[i] <-shortest_dist(mut_a=mut_pairs[i,1],mut_b=mut_pairs[i,2],tree=tree,ref=post_ind[[1]][1:2],node_comb_dist=node_comb_dist)
      if(i%%1000==0) {print(i)}
    }
    
    ##(4) Convert the "distance" into a "similarity", by taking the inverse.
    mut_pairs$sim=1/(mut_pairs$dist+1) #Add one to denominator to avoid dividing by 0, and make max of 1.
    
    
    ##Apply this across all possible tissue pairs
    cat("Calculating the cosine similarity across all tissue pairs",sep="\n")
    tissue_pairs=as.data.frame(data.table::CJ(tissues,tissues))
    colnames(tissue_pairs) <- c("tissue_1","tissue_2")
    tissue_pairs$soft_cosim = NA
    for(i in 1:nrow(tissue_pairs)) {
      if(i%%10==0) {print(i)}
      tissue_pairs$soft_cosim[i] <-get_soft_cosim(data=t(embryonic_median_cell_fracs),selected_muts=embryonic_muts,tissue_1 = tissue_pairs[i,1],tissue_2=tissue_pairs[i,2],mut_pairs_sim = mut_pairs$sim)
    }
    
    #(6) Convert to wide format for the heatmap
    sim_mat<-tissue_pairs%>%
      tidyr::pivot_wider(names_from = tissue_2,values_from = soft_cosim) %>%
      dplyr::select(-1) %>%
      as.matrix()
    rownames(sim_mat) <- colnames(sim_mat)
    
    cat("Completed similarity matrix",sep="\n")
    similarity_matrices[[j]]<-sim_mat
  }
  
  saveRDS(object=similarity_matrices,file = similarity_matrices_file)
  
} else {
  similarity_matrices=readRDS(similarity_matrices_file)
}

#========================================#
# VISUALIZE THE RESULTS AS HEATMAPS OF SIMILARITY ####
#========================================#

#First bind the results into a single tidy dataframe for easy visualization with ggplot2
all_tidy_df<-Map(sim_mat=similarity_matrices,pair=names(all.trees.cc.nodups),function(sim_mat,pair) {
  cat(pair,sep="\n")
  if(!is.null(sim_mat)){
    #diag(sim_mat)<-NA
    tidy_df<-as.data.frame(sim_mat)%>%
      tibble::rownames_to_column(var = "Cell_type1")%>%
      tidyr::gather(-Cell_type1,key="Cell_type2",value="similarity")%>%
      mutate(Pair=pair,.before=1)
    return(tidy_df)
  } else {
    return(NULL)
  }
})%>%dplyr::bind_rows()

#Define the levels to order by cell type
Cell_type_by_type_levels=c("Donor Monocytes","Recipient Monocytes","Donor B cells","Recipient B cells","Donor T cells","Recipient T cells")

## Generate Extended Data Fig. 13d ----
#Plot grouped by cell type (as facets with a uniform colour scale)
sim_matrices_groupby_celltype_plot<-all_tidy_df%>%
  filter(!grepl(pattern = "Granulocytes",Cell_type1)&!grepl(pattern = "Granulocytes",Cell_type2))%>%
  mutate(Cell_type1=factor(gsub("_"," ",Cell_type1),levels=Cell_type_by_type_levels),
         Cell_type2=factor(gsub("_"," ",Cell_type2),levels=Cell_type_by_type_levels),
         pair_new=factor(new_pair_names[Pair],levels=paste0("Pair_",1:10)))%>%
  ggplot(aes(x = Cell_type1, y = Cell_type2, fill = similarity)) +
  geom_tile()+
  annotate("rect",xmin = 0.5,xmax=2.5,ymin=0.5,ymax=2.5,col="black",linewidth=1,fill=NA)+
  annotate("rect",xmin = 2.5,xmax=4.5,ymin=2.5,ymax=4.5,col="black",linewidth=1,fill=NA)+
  annotate("rect",xmin = 4.5,xmax=6.5,ymin=4.5,ymax=6.5,col="black",linewidth=1,fill=NA)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
  facet_wrap(~pair_new,nrow = 2)+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90,hjust=1))+
  labs(fill="Soft cosine\nsimilarity",x="",y="")
ggsave(filename = paste0(plots_dir,"ExtDatFig13d.pdf"),sim_matrices_groupby_celltype_plot,width=7,height=4)

