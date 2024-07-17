##BENCH-MARKING OF PHYLOGENY STRUCTURE
##In addition to the bootstrapping which is in a separate script, this script
#1.  Compares the tree to that built by other tree-building algorithms
#2.  Reviews internal consistency of the raw tree genotypes by two methods: (1) testing against perfect assumptions of phylogeny mutations, (2) comparing expected genotypes from the consensus phylogeny vs the actual genotype in the input matrix

library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotrix)
library(phangorn)
library(pheatmap)
library(RColorBrewer)

#Define function required later in script for comparing phylogenies of different algorithms
comparePhylo_and_plot=function(tree1,tree2,names){
  plot_comp_tree=function(tree,comp,title,col="red",lwd=1,tree_pos){
    tree_name=deparse(substitute(tree))
    shared_clades=comp[,tree_pos]
    edge_width=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),lwd,2*lwd))
    edge_col=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),"black",col))
    plot(tree,show.tip.label=F,direction="downwards",edge.color=edge_col,edge.width=edge_width,main=title)
  }
  comp<-compare_nodes(tree1,tree2)
  par(mfrow=c(1,2))
  plot_comp_tree(tree1,comp=comp,title=names[1],tree_pos = 1)
  plot_comp_tree(tree2,comp=comp,title=names[2],tree_pos = 2)
}

compare_nodes=function(x, y) 
{
  tree1 <- deparse(substitute(x))
  tree2 <- deparse(substitute(y))
  n1 <- Ntip(x)
  n2 <- Ntip(y)
  
  key1 <- makeNodeLabel(x, "md5sum")$node.label
  key2 <- makeNodeLabel(y, "md5sum")$node.label
  mk12 <- match(key1, key2)
  mk21 <- match(key2, key1)
  if (any(tmp <- is.na(mk12))) {
    nk <- sum(tmp)
  }
  if (any(tmp <- is.na(mk21))) {
    nk <- sum(tmp)
  }
  nodes1 <- which(!is.na(mk12))
  nodes2 <- mk12[!is.na(mk12)]
  
  NODES <- data.frame(nodes1 + n1,nodes2 + n2)
  names(NODES) <- c(tree1, tree2)
  return(NODES)
}

get_RF_dist=function(tree1,tree2){
  comp<-compare_nodes(tree1,tree2)
  RF_dist=1-(length(comp[,1])/length(unique(tree1$edge[,1])))
  return(RF_dist)
}

#Set file paths
my_working_directory=paste0("/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/final_versions/")
output_dir=my_working_directory

#Set working directly & load-up required functions
setwd(my_working_directory)
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

all_pairs=paste0("Pair",c(11,13,21,24,25,28,31,38,40,41))

for(exp_ID in all_pairs) {
  cat(exp_ID,sep="\n")
  
  filtered_muts_files = list.files(paste0(my_working_directory,"annot_files/"),pattern="annot",full.names = T)
  tree_file_paths=list.files(paste0(my_working_directory,"tree_files/"),pattern=".tree",full.names = T)
  
  filtered_muts_file=grep(exp_ID,filtered_muts_files,value=T)
  tree_file_path=grep(exp_ID,tree_file_paths,value=T)
  
  #IQTree paths
  iqtree_path="/lustre/scratch126/casm/team154pc/ms56/programs/IQ-TREE/build/iqtree"
  iqtree_input_file_path=paste0("iqtree_fasta_",exp_ID,".fa")
  iqtree_output_tree_path=paste0(iqtree_input_file_path,".treefile")
  
  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  load(filtered_muts_file)
  gt<-filtered_muts$Genotype_shared_bin
  NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  details<-filtered_muts$COMB_mats.tree.build$mat
  tree<-drop.tip(read.tree(tree_file_path),"Ancestral")
  tree_mod<-di2multi(tree)
  tree_mod$edge.length[tree_mod$edge.length<10]<-10
  
  #COMPARE WITH OTHER PHYLOGENY-BUILDING ALGORITHMS
  #1. iqtree. Use the binary input model with equal likelihood for all sites.
  #iqtree uses a maximum-likelihood approach
  if(!file.exists(iqtree_output_tree_path)){
    #Write the fasta file in binary format for iqtree
    write.fasta(lapply(filtered_muts$dna_strings,function(x) gsub("W","0",gsub("V","1",x))),names=names(filtered_muts$dna_strings),iqtree_input_file_path)
    system(paste0(iqtree_path," -s ",iqtree_input_file_path," -m JC2 -bb 1000 -czb -redo")) #-czb option allows collapse to polytomies
  }

  #Read in the output tree
  tree_iq=read.tree(iqtree_output_tree_path)
  tree_iq <- drop.tip(tree_iq,"Ancestral")
  tree_iq<-di2multi(tree_iq)
  tree_iq$edge.length = rep(1, nrow(tree_iq$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work

  #Assign mutations back to the tree
  df = reconstruct_genotype_summary(tree_iq) #Define df (data frame) for treeshape
  res = assign_to_tree(NV[,df$samples],NR[,df$samples],df,error_rate = c(rep(0.01, ncol(NR))))
  tree_iq$edge.length <- res$df$df$edge_length #Assign edge lengths from the res object
  tree_iq<-di2multi(tree_iq)

  #Versions of the tree where short branches are lengthened to a 10 to allow visualization of differences
  tree_iq_mod<-tree_iq
  tree_iq_mod$edge.length[tree_iq_mod$edge.length<10]<-10

  pdf(paste0(output_dir,"/MPBoot_vs_iqtree_comparisons_",exp_ID,".pdf"),width = 15,height=6)
  comparePhylo_and_plot(di2multi(tree),tree_iq,names=c("MPBoot phylogeny","IQTree phylogeny"))
  dev.off()

  pdf(paste0(output_dir,"/MPBoot_vs_iqtree_comparisons_mod_",exp_ID,".pdf"),width = 15,height=6)
  comparePhylo_and_plot(di2multi(tree_mod),tree_iq_mod,names=c("MPBoot phylogeny","IQTree phylogeny"))
  dev.off()

  #Calculate and print the similarity metrics for the iqtree tree
  tree_iq_splits_comp<-Quartet::SplitStatus(tree_iq,di2multi(tree))
  tree_iq_quartets_comp<-Quartet::QuartetStatus(tree_iq,di2multi(tree))
  print(paste("The Robinson-Foulds similarity is",Quartet::SimilarityMetrics(tree_iq_splits_comp)[,"SymmetricDifference"]))
  print(paste("The Quartet Similarity is",Quartet::SimilarityMetrics(tree_iq_quartets_comp)[,"SymmetricDifference"]))


  #2.MEASURING INTERNAL CONSISTENCY

  #(a) display heatmap of mutation genotype to get feel for the data
  pdf(paste0(output_dir,"/Shared_mutation_genotypes_heatmap_",exp_ID,".pdf"),width = 8,height=10)
  max_mut=2500
  if(dim(gt)[1]>max_mut) {
    gt_sub<-gt[sample(1:nrow(gt),size=max_mut,replace = F),]
    pheatmap(gt_sub,legend=T,show_colnames=F,show_rownames=F,color = c("#08519C","light grey","#CB181D"))
  } else {
    pheatmap(gt,legend=T,show_colnames=F,show_rownames=F,color = c("#08519C","light grey","#CB181D"))
  }
  dev.off()
}

for(exp_ID in all_pairs) {
  cat(exp_ID,sep="\n")
  
  filtered_muts_files = list.files(paste0(my_working_directory,"annot_files/"),pattern="annot",full.names = T)
  tree_file_paths=list.files(paste0(my_working_directory,"tree_files/"),pattern=".tree",full.names = T)
  
  filtered_muts_file=grep(exp_ID,filtered_muts_files,value=T)
  tree_file_path=grep(exp_ID,tree_file_paths,value=T)
  
  #SCITE paths
  scite_input_file_path=paste0("scite_input_",exp_ID)
  scite_path="/lustre/scratch126/casm/team154pc/ms56/programs/SCITE/scite"
  scite_output_tree_path=paste0(scite_input_file_path,"_ml0.newick")
  scite_allocated_muts_tree_path=paste0("scite_tree_allocated_muts_",exp_ID,".tree")

  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  load(filtered_muts_file)
  gt<-filtered_muts$Genotype_shared_bin
  NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  details<-filtered_muts$COMB_mats.tree.build$mat
  tree<-drop.tip(read.tree(tree_file_path),"Ancestral")
  tree_mod<-di2multi(tree)
  tree_mod$edge.length[tree_mod$edge.length<10]<-10
  
  #2. SCITE
  #SCITE has a binary input format where 0 is absent, 1 is present, and 3 is "missing data".
  gt=filtered_muts$Genotype_shared_bin
  if(!file.exists(scite_output_tree_path)){
    #Write file in required format for SCITE
    gt_rows=apply(gt,1,paste,collapse=" ") #Collapse into strings
    gt_rows=sapply(gt_rows, function(x) gsub("0.5","3", x)) #Replace the 0.5's with 3's as per the SCITE input
    writeLines(gt_rows,scite_input_file_path)
    
    #Write the SCITE command
    scite_command=paste0(scite_path," -i ",scite_input_file_path," -n ",nrow(gt)," -m ",ncol(gt)," -r 1 -l 1000000 -fd 0.001 -ad 0.005 -transpose")
    system(paste0("bsub -o $PWD/log.%J -e $PWD/err.%J -q basement -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J scite_",exp_ID," ",scite_command)) #SCITE takes a long time: days to weeks
    
  } else {
    tree_scite=read.tree(text=paste0(readLines(scite_output_tree_path),";"))
    tree_scite$edge.length=rep(1,nrow(tree_scite$edge))
    tree_scite$tip.label<-sapply(tree_scite$tip.label,function(x) colnames(gt)[as.numeric(x)]) #Map back the sample names from the numbers
    
    #Drop tips that aren't in the original tree
    tree_scite<-drop.tip(tree_scite,colnames(gt)[!colnames(gt)%in%colnames(NV)])
    
    #Assign mutations back to the tree
    df = reconstruct_genotype_summary(tree_scite) #Define df (data frame) for treeshape
    res = assign_to_tree(NV[,df$samples],NR[,df$samples],df,error_rate = c(rep(0.01, ncol(NR))))
    tree_scite$edge.length <- res$df$df$edge_length #Assign edge lengths from the res object
    tree_scite<-di2multi(tree_scite)
    write.tree(tree_scite,file=scite_allocated_muts_tree_path)
    
    tree_scite_mod<-tree_scite
    tree_scite_mod$edge.length[tree_scite_mod$edge.length<10]<-10
    
    pdf(paste0(output_dir,"/MPBoot_vs_SCITE_comparisons_",exp_ID,"_",".pdf"),width = 15,height=6)
    comparePhylo_and_plot(di2multi(tree),tree_scite,names=c("MPBoot phylogeny","SCITE phylogeny"))
    dev.off()
    
    pdf(paste0(output_dir,"/MPBoot_vs_SCITE_comparisons_mod_",exp_ID,"_",".pdf"),width = 15,height=6)
    comparePhylo_and_plot(di2multi(tree_mod),tree_scite_mod,names=c("MPBoot phylogeny","SCITE phylogeny"))
    dev.off()
    
    #Calculate and print the similarity metrics for the scite tree
    tree_scite_splits_comp<-Quartet::SplitStatus(tree_scite,di2multi(tree))
    tree_scite_quartets_comp<-Quartet::QuartetStatus(tree_scite,di2multi(tree))
    print(paste("The Robinson-Foulds similarity is",round(Quartet::SimilarityMetrics(tree_scite_splits_comp)[,"SymmetricDifference"],digits=5)))
    print(paste("The Quartet Similarity is",round(Quartet::SimilarityMetrics(tree_scite_quartets_comp)[,"SymmetricDifference"],digits = 5)))
  }
}


##CALCULATION OF THE DISAGREEMENT SCORES ACROSS THE PHYLOGENIES

#Function for calculating the disagreement scores
#Note that uncertain genotypes should be recorded as '0.5', positives as '1', and negatives as '0'

#Set working directly & load-up required functions
setwd(my_working_directory)
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)


disagreement_score=function(binary_mut_mat,max_mut=2500) {
  if(dim(binary_mut_mat)[1]>max_mut) {
    binary_mut_mat<-binary_mut_mat[sample(1:nrow(binary_mut_mat),size=max_mut,replace = F),]
  }
  c=matrix(as.logical(binary_mut_mat),ncol=ncol(binary_mut_mat)) #c is a binary matrix with 0.5's rounded up to 1's
  d<-binary_mut_mat;d[d==0.5]<-0;d<-matrix(as.logical(d),ncol=ncol(d)) #d is a binary matrix with 0.5's rounded down to 0's
  
  #Do the two 'perfect phylogeny tests' using these matrices
  #1. Test for exclusivity of each mutation pair
  e_exc=(d)%*%t(!c)
  
  #2. Test for inclusivity of each mutation pair
  e_inc=(d)%*%t(d) #Need to use the version where 0.5 is cooerced to 0
  
  #All mutation pairs should meet either one of these two tests, therefore take the minimum of each
  min_exc=matrix(mapply("min",e_exc,t(e_exc)),ncol=ncol(e_exc))
  out=matrix(mapply("min",min_exc,e_inc),ncol=ncol(e_exc))
  mean(out)
}

disagreement_score_list=vector("list",length=length(all_pairs))
names(disagreement_score_list)<-all_pairs
nrand=5 #The number of 'random shuffles' to perform as a control
max_mut=10000

for(exp_ID in all_pairs) {
  cat(exp_ID,sep="\n")
  filtered_muts_files = list.files(paste0(my_working_directory,"annot_files/"),pattern="annot",full.names = T)
  tree_file_paths=list.files(paste0(my_working_directory,"tree_files/"),pattern=".tree",full.names = T)
  
  #Get the file paths
  filtered_muts_file=grep(exp_ID,filtered_muts_files,value=T); tree_file_path=grep(exp_ID,tree_file_paths,value=T)

  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  cat("Loading up mutations file",sep="\n")
  load(filtered_muts_file)
  gt<-filtered_muts$Genotype_shared_bin
  cat(paste("The shared mutation matrix includes",nrow(gt),"mutations. This will be reduced to",max_mut,"mutations"),sep="\n")
  NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  details<-filtered_muts$COMB_mats.tree.build$mat
  cat("Loading up tree file",sep="\n")
  tree<-drop.tip(read.tree(tree_file_path),"Ancestral")
  
  #Comparing the score to random shuffles at each locus
  cat("Calculating disagreement score of the data",sep="\n")
  dat=disagreement_score(filtered_muts$Genotype_shared_bin,max_mut=max_mut)

  cat(paste("Calculating disagreement score of",nrand,"random shuffles"),sep="\n")
  rand=sapply(1:nrand,function(i) {print(i);gt_shuffled<-apply(filtered_muts$Genotype_shared_bin,1,FUN = base::sample);disagreement_score(gt_shuffled,max_mut = max_mut)})
  my_theme=theme_classic(base_size=7,base_family="Helvetica")+theme(text=element_text(size=7,family="Helvetica"))
  
  disagreement_score_list[[exp_ID]]<-data.frame(type=c("data",rep("random_shuffles",nrand)),disagreement_score=c(dat,rand))
}

out_df=Map(df=disagreement_score_list,exp_ID=all_pairs,f=function(df,exp_ID) {
  return(df%>%mutate(pair=exp_ID,.before=1))
})%>%dplyr::bind_rows()

breaks <- 10^(-3:3)
minor_breaks <- c()

plot3<-out_df%>%
  filter(type=="random_shuffles")%>%
  ggplot(aes(y=pair,x=disagreement_score,col=type))+
  geom_jitter(size=0.15,width = 0.05,alpha=0.5)+
  geom_point(data=out_df%>%filter(type=="data"))+
  scale_x_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=minor_breaks)+
  theme_bw()+
  labs(x="Disagreement score",y="Pair",col="Type")+
  coord_flip()

plot3<-out_df%>%
  filter(type=="random_shuffles")%>%
  ggplot(aes(x=pair,y=disagreement_score,col=type))+
  geom_jitter(size=0.5,height = 0.05,alpha=0.5)+
  geom_point(data=out_df%>%filter(type=="data"))+
  scale_y_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))+
  labs(x="",y="Disagreement score",col="")+
  annotation_logticks(sides="l",short=unit(0.05,"cm"),mid=unit(0.05,"cm"),long=unit(0.1,"cm"))

ggsave(plot3,filename = "Disagreement_score.pdf",device="pdf",width=4,height=2.5)

# #--------------------------------------------------------------------------
# #(b) Comparing pairs of loci for perfect phylogeny assumptions (does not rely on the phylogeny)
# #Creates the "perfect" mutation matrix that would expect from the tree
# gt<-gt[rownames(gt)%in%details$mut_ref,]
# tree<-read.tree(tree_file_path)
# phylogeny_mut_mat=Reduce(rbind,lapply(rownames(gt),function(mut) colnames(gt)%in% getTips(tree,details$node[details$mut_ref==mut])))
# dimnames(phylogeny_mut_mat)=dimnames(gt)
# 
# #Sum the number of times that the mutation should be positive (according to the tree) but is negative in the genotype matrix
# unexpected_neg=sum(phylogeny_mut_mat*!as.logical(gt)) #0.5s will be coerced to "TRUE", which will then not count towards score
# 
# #Sum the number of times that the mutation should be negative (according to the tree) but is positive in the genotype matrix
# gt_min<-gt; gt_min[gt==0.5]<-0 #Create version where the 0.5s are negative
# unexpected_pos=sum((!phylogeny_mut_mat)*as.logical(gt_min))
# 

# #Review the inconsistent mutations to assess cause of these inconsistencies
# out=matrix(phylogeny_mut_mat*!as.logical(gt),ncol=ncol(phylogeny_mut_mat))
# out=matrix(((!as.logical(phylogeny_mut_mat))*as.logical(gt_min)),ncol=ncol(phylogeny_mut_mat))
# dimnames(out)=dimnames(phylogeny_mut_mat)
# muts=names(rowSums(out)[which(rowSums(out)>0)])
# mut=muts[6]
# details[details$mut_ref==mut,]
# 
# tree=plot_tree(tree,cex.label=0)
# add_annotation(tree,
#                details=details,
#                matrices=list(NV=NV, NR=NR),
#                annot_function = plot_MAV_mut,
#                mut1=mut,
#                lesion_node=details$node[details$mut_ref==mut],
#                cex=0.4)
# add_annotation(tree,details,matrices=NULL,annot_function = highlight_nodes,nodes=details$node[details$mut_ref==mut])
