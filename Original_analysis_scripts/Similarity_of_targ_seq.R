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
bulk_smry_all<-bulk_smry_all%>%filter(tissueID!="PD45813i" &tissueID!="PD45812i")

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

##One way of assessing similarity of clonal composition between tissues is using the
##soft cosine similarity of the targeted sequencing data between tissues/ individuals

##Select objects from the individual of interest
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

####
embryonic_mutation_cutoff=30
similarity_matrices_file=paste0(root_dir,"/data/similarity_matrices.Rds")
if(!file.exists(similarity_matrices_file)) {
  similarity_matrices=list()
  for(j in 1:length(posterior_cell_fracs)) {
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
    nmut_max=600
    cat(paste("There are",length(embryonic_muts),"embryonic mutations"),sep="\n")
    if(length(embryonic_muts)>nmut_max) {
      cat(paste("Randomly subsampling to",nmut_max,"mutations"),sep="\n")
      embryonic_muts<-sample(embryonic_muts,size=nmut_max)
      embryonic_median_cell_fracs<-embryonic_median_cell_fracs[embryonic_muts,]
    }
    cat(paste("There are",length(embryonic_muts),"embryonic mutations included"),sep="\n")
    
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
      pivot_wider(names_from = tissue_2,values_from = soft_cosim) %>%
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


#VISUALIZE THE RESULTS AS HEATMAPS OF SIMILARITY
all_tidy_df<-Map(sim_mat=similarity_matrices,pair=names(all.trees.cc.nodups),function(sim_mat,pair) {
  cat(pair,sep="\n")
  if(!is.null(sim_mat)){
    #diag(sim_mat)<-NA
    tidy_df<-as.data.frame(sim_mat)%>%
      tibble::rownames_to_column(var = "Cell_type1")%>%
      gather(-Cell_type1,key="Cell_type2",value="similarity")%>%
      mutate(Pair=pair,.before=1)
    return(tidy_df)
  } else {
    return(NULL)
  }
})%>%dplyr::bind_rows()

Cell_type_levels=c("Donor_Monocytes","Donor_B_cells","Donor_T_cells","Recipient_Monocytes","Recipient_B_cells","Recipient_T_cells")
Cell_type_by_type_levels=c("Donor_Monocytes","Recipient_Monocytes","Donor_B_cells","Recipient_B_cells","Donor_T_cells","Recipient_T_cells")


#First plot as facets with a uniform colour scale
sim_matrices_plot<-all_tidy_df%>%
  filter(!grepl(pattern = "Granulocytes",Cell_type1)&!grepl(pattern = "Granulocytes",Cell_type2))%>%
  mutate(Cell_type1=factor(Cell_type1,levels=Cell_type_levels),
         Cell_type2=factor(Cell_type2,levels=Cell_type_levels),
         pair_new=factor(new_pair_names[Pair],levels=paste0("Pair_",1:10)))%>%
  ggplot(aes(x = Cell_type1, y = Cell_type2, fill = similarity)) +
  geom_tile()+
  annotate("rect",xmin = 0.5,xmax=3.5,ymin=0.5,ymax=3.5,col="black",linewidth=1,fill=NA)+
  annotate("rect",xmin = 3.5,xmax=6.5,ymin=3.5,ymax=6.5,col="black",linewidth=1,fill=NA)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
  facet_wrap(~pair_new,nrow = 2)+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90))
ggsave(filename = paste0(plots_dir,"similarity_heatmaps.pdf"),sim_matrices_plot,width=7,height=4)

#Now display as separate plots with separate scales
plot_list<-lapply(paste0("Pair_",1:10),function(pair) {
  pair_plot<-all_tidy_df%>%
    filter(!grepl(pattern = "Granulocytes",Cell_type1)&!grepl(pattern = "Granulocytes",Cell_type2))%>%
    mutate(Cell_type1=factor(Cell_type1,levels=Cell_type_levels),
           Cell_type2=factor(Cell_type2,levels=Cell_type_levels),
           pair_new=factor(new_pair_names[Pair],levels=new_pair_names))%>%
    filter(pair_new==pair)%>%
    ggplot(aes(x = Cell_type1, y = Cell_type2, fill = similarity)) +
    geom_tile()+
    annotate("rect",xmin = 0.5,xmax=3.5,ymin=0.5,ymax=3.5,col="black",linewidth=1,fill=NA)+
    annotate("rect",xmin = 3.5,xmax=6.5,ymin=3.5,ymax=6.5,col="black",linewidth=1,fill=NA)+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
    facet_wrap(~pair_new,nrow = 2)+
    theme_classic()+
    my_theme+
    theme(axis.text.x = element_text(angle=90))
  return(pair_plot)
})
ggsave(filename = paste0(plots_dir,"similarity_heatmaps_separate.pdf"),gridExtra::arrangeGrob(grobs=plot_list,ncol=5),width=15,height=6)


#Now repeat but grouped by cell type (as facets with a uniform colour scale)
sim_matrices_groupby_celltype_plot<-all_tidy_df%>%
  filter(!grepl(pattern = "Granulocytes",Cell_type1)&!grepl(pattern = "Granulocytes",Cell_type2))%>%
  mutate(Cell_type1=factor(Cell_type1,levels=Cell_type_by_type_levels),
         Cell_type2=factor(Cell_type2,levels=Cell_type_by_type_levels),
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
  theme(axis.text.x = element_text(angle=90))
ggsave(filename = paste0(plots_dir,"similarity_heatmaps_groupby_celltype.pdf"),sim_matrices_groupby_celltype_plot,width=7,height=4)

#Group by cell type (with separate colour scales)
plot_groupby_celltype_list<-lapply(paste0("Pair_",1:10),function(pair) {
  pair_plot<-all_tidy_df%>%
    filter(!grepl(pattern = "Granulocytes",Cell_type1)&!grepl(pattern = "Granulocytes",Cell_type2))%>%
    mutate(Cell_type1=factor(Cell_type1,levels=Cell_type_by_type_levels),
           Cell_type2=factor(Cell_type2,levels=Cell_type_by_type_levels),
           pair_new=factor(new_pair_names[Pair],levels=paste0("Pair_",1:10)))%>%
    filter(pair_new==pair)%>%
    ggplot(aes(x = Cell_type1, y = Cell_type2, fill = similarity)) +
    geom_tile()+
    annotate("rect",xmin = 0.5,xmax=2.5,ymin=0.5,ymax=2.5,col="black",linewidth=1,fill=NA)+
    annotate("rect",xmin = 2.5,xmax=4.5,ymin=2.5,ymax=4.5,col="black",linewidth=1,fill=NA)+
    annotate("rect",xmin = 4.5,xmax=6.5,ymin=4.5,ymax=6.5,col="black",linewidth=1,fill=NA)+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
    facet_wrap(~pair_new,nrow = 2)+
    theme_classic()+
    my_theme+
    theme(axis.text.x = element_text(angle=90))
  return(pair_plot)
})
ggsave(filename = paste0(plots_dir,"similarity_heatmaps_by_celltype_separate.pdf"),gridExtra::arrangeGrob(grobs=plot_groupby_celltype_list,ncol=5),width=15,height=6)

#Now the same but clustering cell types by similarity

library(ggplot2)
library(reshape2)
library(ggdendro)


x<-similarity_matrices[[6]][Cell_type_levels,Cell_type_levels]
dd.col <- as.dendrogram(hclust(dist(x)))
col.ord <- order.dendrogram(dd.col)

dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)

xx <- x[col.ord, row.ord]
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xx)
colnames(df) <- xx_names[[2]]
df$Cell_type1 <- xx_names[[1]]
df$Cell_type1 <- with(df, factor(Cell_type1, levels=Cell_type1, ordered=TRUE))

mdf <- melt(df, id.vars="Cell_type1",variable.name = "Cell_type2",value.name = "Similarity")

ddata_x <- dendro_data(dd.row)
ddata_y <- dendro_data(dd.col)

### Set up a blank theme
theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank()
)

### Create plot components ###    
# Heatmap
granularity=0.1
full_colscale<-colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(71)
global_min=0.65
scale_min=1+round(min(mdf$Similarity)*200)-round(global_min*200)
p1 <- ggplot(mdf, aes(x=Cell_type1, y=Cell_type2)) + 
  geom_tile(aes(fill=Similarity)) + scale_fill_gradientn(colours = full_colscale[scale_min:length(full_colscale)])+
  theme(axis.text.x=element_text(angle=90),
        axis.title=element_blank(),legend.position = "none")

# Dendrogram 1
p2 <- ggplot(segment(ddata_x)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  theme_none + theme(axis.title.x=element_blank())

# Dendrogram 2
p3 <- ggplot(segment(ddata_y)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  coord_flip() + theme_none

grid::grid.newpage()
plot(p1, vp=viewport(0.8, 0.8,x=0.5,y=0.4))
plot(p2, vp=viewport(0.48,0.2,x=0.63,y=0.9))

