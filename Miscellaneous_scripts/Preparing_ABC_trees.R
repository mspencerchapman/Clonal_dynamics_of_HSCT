#GENERATE THE TREES ON WHICH TO DO THE ABC

#1. Import the SNV trees that are corrected for sequencing coverage
#2. Reduce branch lengths by the proportion caused by non-ubiquitous mutational process i.e. APOBEC/ mutagenic exposures
#3. Reduce terminal branch lengths by the number of 'in vitro' or 'non linear' (e.g. mutations caused from increased replication through cell division during differentiation) additional mutations
#4. Then make ultrametric, scaling branches to the average mutation burden across all branches
#5. Estimate timing of coalescences during HSC life based on (1) 60 mutations at birth, (2) linear accumulation thereafter
#6. Create 'transplant time' shaded box based on time + poisson variation of mean mutation burden at that time

#Functions required for script
get_DR_ids=function(tree){
  tree=ape::drop.tip(tree,"Ancestral")
  PD_IDs=unique(substr(tree$tip.label,1,8))
  #Get number elements only
  PD_numbers=readr::parse_number(PD_IDs)
  names(PD_IDs)<-sapply(PD_numbers,function(n) ifelse(n%%2==0,"donor_ID","recip_ID")) #This relies on the fact that all donor IDs are even numbers, all recip IDs are odd
  return(PD_IDs)
}
##Now convert mutations to age
age_from_muts=function(n_muts,tree,sampling_age,birth_muts=50) {
  sampling_age_muts=mean(get_mut_burden(tree))
  if(n_muts<=birth_muts) {
    return(0)
  } else {
    age=sampling_age * (n_muts-birth_muts)/(sampling_age_muts-birth_muts)
    return(age)
  }
}
muts_from_age=function(age,tree,sampling_age,birth_muts=50) {
  sampling_age_muts=mean(get_mut_burden(drop.tip(tree,"Ancestral")))
  n_muts=birth_muts + (age*(sampling_age_muts-birth_muts)/sampling_age)
  return(n_muts)
}

#Function to import the HDP data and convert into an exposures data frame
generate_exposures_df=function(HDP_multi_chain_RDS_path,trinuc_mut_mat_path,key_table_path){
  mut_example_multi=readRDS(HDP_multi_chain_RDS_path)
  mutations=read.table(trinuc_mut_mat_path)
  key_table=read.table(key_table_path)
  
  sample_remove=rownames(mutations)[rowSums(mutations)<50]
  mutations=mutations[!rownames(mutations)%in%sample_remove,]
  key_table=key_table[!key_table$Sample%in%sample_remove,]
  freq=nrow(mutations)
  
  dp_distn <- comp_dp_distn(mut_example_multi)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),,drop=FALSE])
  colnames(exposures)=rownames(mutations)
  rownames(exposures)<-paste0("N",rownames(exposures))
  sigs=rownames(exposures)
  sig_profiles=mut_example_multi@comp_categ_distn$mean
  
  exposures_df<-as.data.frame(t(exposures),stringsAsFactors=F)%>%
    tibble::rownames_to_column("branch")%>%
    tidyr::separate(col="branch",into=c("node","exp_ID"),sep="_")%>%
    mutate(node=as.numeric(node))
  
  return(exposures_df)
}

#Function to extract average signature contributions in each sample
#The tree should be the corrected, non-ultrametric tree, and the exposures df should have unique node numbers 
get_signatures_in_samples=function(tree,signature_names,exposures_df) {
  sigs_in_samples=sapply(tree$tip.label,function(sample) {
    sample_node=which(tree$tip.label==sample)
    nodes_included=get_ancestral_nodes(node = sample_node,edge = tree$edge,exclude_root = T)
    branch_lengths=sapply(nodes_included,function(node) tree$edge.length[tree$edge[,2]==node])
    branch_sig_prop=sapply(nodes_included,function(node) {ifelse(node%in%exposures_df$node,sum(exposures_df[exposures_df$node==node,signature_names]),NA)})
    overall_contribution=weighted.mean(x=branch_sig_prop[!is.na(branch_sig_prop)],w = branch_lengths[!is.na(branch_sig_prop)])
    if(is.nan(overall_contribution)){overall_contribution<-0}
    return(overall_contribution)
  })
}

#START OF SCRIPT
my_working_directory<-getwd()
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch119/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
setwd(tree_mut_dir); source("treemut.R");setwd(my_working_directory)

#Manually enter the metadata data frame
exp_nos<-c(11,13,21,24,25,28,31,38,40,41)
Pair_metadata=data.frame(Pair=paste0("Pair",exp_nos),
                         Age=c(74.8,65.5,64.5,34.2,58.4,79.9,65.2,65.8, 51.9,42.4),
                         Age_at_transplant=c(66,36,50,18,47,63,43,35,35,30),
                         MNC_dose=c(2.66,4.17,16.28,NA,10.94,13.9,NA,4.05,2.48,15.01),
                         CD34_dose=c(1.56,NA,7.1,7.9,8.97,2.4,NA,NA,NA,4.51),
                         stem_cell_source=c("BM","BM","PBSC","PBSC","PBSC","PBSC","BM","BM","BM","PBSC"),
                         conditioning=c("MAC","MAC","RIC","MAC","RIC","RIC","MAC","MAC","MAC","MAC"))

library(hdp)
library(ape)
library(dplyr)
library(ggplot2)
library(tidyr)

#Import the signature information from HDP
# HDP_folder="/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/HDP/"
# tree_folder="/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/filtering_runs2/tree_files"
# annotated_muts_folder="/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/filtering_runs2/annotated_muts"

HDP_folder="~/R_work/Zurich_HSCT/Mutational_signature_extraction/HDP/HDP_280621"

filtering_type="vaf"
if(filtering_type=="vaf"){
  tree_folder="~/R_work/Zurich_HSCT/Data/vaf_filtering/tree_files"
  annotated_muts_folder="~/R_work/Zurich_HSCT/Data/vaf_filtering/annot_files"
}else if(filtering_type=="pval"){
  tree_folder="~/R_work/Zurich_HSCT/Data/pval_filtering/tree_files"
  annotated_muts_folder="~/R_work/Zurich_HSCT/Data/pval_filtering/annot_files"
}

LOY_files=list.files(path="~/R_work/Zurich_HSCT/Data/",pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))
CN_change_df=read.csv("~/R_work/Zurich_HSCT/Data/Copy_number_changes.csv")

exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)

##1. Import the SNV trees corrected for coverage & the details matrices
tree_paths=list.files(tree_folder,pattern=ifelse(filtering_type=="pval","pval_post_mix_post_dup.tree","vaf_post_mix_post_dup.tree"),full.names = T)
all.trees<-lapply(exp_nos,function(exp_no){read.tree(grep(paste0("Pair",exp_no),tree_paths,value = T))})
names(all.trees)<-paste0("Pair",exp_nos)

annotated_muts_paths=list.files(annotated_muts_folder,pattern=ifelse(filtering_type=="pval","pval_post_mix_post_dup","vaf_post_mix_post_dup"),full.names = T)
all.muts<-lapply(exp_nos,function(exp_no){load(grep(paste0("Pair",exp_no),annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
names(all.muts)<-paste0("Pair",exp_nos)

#Create dataframe of the LOY events
Pair11_LOY_nodes=c(447, 468, 152, 407, 369, 56, 493, 539)
Pair11_loss_of_Y_details=data.frame(Chrom="Y",Pos=NA,Ref=NA,Alt=NA,mut_ref=paste0("LOY_",1:length(Pair11_LOY_nodes)),
                                    Mut_type="CNA",node=Pair11_LOY_nodes,pval=NA,Gene="LOY",Transcript="",RNA="",CDS="",
                                    Protein="",Type="",SO_codes="",coding_change="Coding change",
                                    coding_change_chip="Coding change mutation in driver",
                                    ChromPos="",variant_ID=paste("LOY",1:length(Pair11_LOY_nodes)))

check_plots=T
if(check_plots){
  par(mfrow=c(5,2))
  temp=lapply(all.trees,plot.phylo,show.tip.label=F,direction="downwards") #Plot the current trees
}

##2. Correct the branch lengths by removing mutations assigned to APOBEC/ platinum signature
all.trees.adj1=lapply(exp_nos,function(exp_no) {
  tree<-all.trees[[paste0("Pair",exp_no)]]
  
  #Now a version with APOBEC/ platinum signatures removed using the HDP output
  remove_sigs=paste0("N",c(4,5)) #Check that these are the right
  exposures_df_pair=exposures_df%>%dplyr::filter(Pair==paste0("Pair",exp_no))
  tree_adj<-tree
  tree_adj$edge.length=sapply(tree$edge[,2],function(node) {
    branch_length_original<-tree$edge.length[tree$edge[,2]==node]
    if(node %in% exposures_df_pair$node) {
      remove_prop<-sum(exposures_df_pair[exposures_df_pair$node==node,remove_sigs])
      return(branch_length_original*(1-remove_prop))
    } else {
      return(branch_length_original)
    }
  })
  return(tree_adj)
})

if(check_plots) {
  par(mfrow=c(4,2))
  temp=lapply(all.trees.adj1,plot.phylo,show.tip.label=F,direction="downwards")
}

##3. Now correct terminal branches for in vitro/ differentiation-induced mutations by removing (in utero mutations - intercept) from linear regression 
in_vitro_mutation_burden=60
all.trees.adj2=lapply(all.trees.adj1,function(tree) {
  tree_adj<-tree
  tree_adj$edge.length<-sapply(1:nrow(tree$edge),function(i) {
    if(tree$edge[i,2]%in%1:length(tree$tip.label)) {
      if(tree$edge.length[i]<in_vitro_mutation_burden){
        print(paste("Branch length reduced to 1 for sample",tree$tip.label[tree$edge[i,2]]))
      }
      return(max(c(tree$edge.length[i] - in_vitro_mutation_burden,1)))
    } else {
      return(tree$edge.length[i])
    }
  })
  return(tree_adj)
})

if(check_plots) {
  par(mfrow=c(4,2))
  temp=lapply(all.trees.adj2,plot.phylo,show.tip.label=F,direction="downwards")
}

##4. Make the trees ultrametric
all.trees.ultra<-lapply(all.trees.adj2,function(tree) {
  mean_mutation_burden=mean(get_mut_burden(tree))
  tree.ultra<-make.ultrametric.tree(tree)
  tree.ultra$edge.length=tree.ultra$edge.length*mean_mutation_burden
  tree.ultra$edge.length[which(tree$edge[,2]==which(tree$tip.label=="Ancestral"))]<-0 #Set ancestral branch length back to 0
  return(tree.ultra)
})
names(all.trees.ultra)<-paste0("Pair",exp_nos)

#Edit the details matrices to have column indicating driver mutations on shared branches only
all.muts<-Map(function(details,tree) {
  mutate(details,shared_coding_change_chip=ifelse(coding_change_chip!="no"&!node%in%1:length(tree$tip.label),"yes","no"))
  },details=all.muts,tree=all.trees)

#Pair41 high mut sample = PD45813b_lo0013, sensitivity df says SNV sens = 0.313. Seems v low for coverage.

#Plot all the adjusted ultrametric trees showing the transplant period
plots_dir="~/R_work/Zurich_HSCT/plots/"
pdf(paste0(plots_dir,"Old_plots.pdf"),width = 20,7)
young_individuals=c("Pair24","Pair41")
middle_age_individuals=c("Pair40","Pair13","Pair21","Pair25","Pair31","Pair38")
old_individuals=c("Pair11","Pair28")
if(check_plots){
  par(mfrow=c(1,2))
  temp=lapply(old_individuals,function(pair) {
    tree<-all.trees.ultra[[pair]]
    details<-all.muts.nodups[[pair]]
    pair_age_at_transplant<-Pair_metadata%>%dplyr::filter(Pair==pair)%>%pull(Age_at_transplant)
    pair_age<-Pair_metadata%>%dplyr::filter(Pair==pair)%>%pull(Age)
    
    #Set up the axis for plotting
    binwidth=ifelse(pair_age>50,10,5)
    max_age_to_plot=binwidth*ceiling(pair_age/binwidth)
    axis_at=c(0,sapply(seq(0,max_age_to_plot,by=binwidth),muts_from_age,tree,pair_age))
    labels_at=c("Zyg.","Birth",seq(binwidth,max_age_to_plot,by=binwidth))
    
    tree=plot_tree(tree,cex.label = 0,plot_axis=F,vspace.reserve=0.1,title=Pair_metadata%>%filter(Pair==pair)%>%pull(Pair_new))
    if(pair%in%c("Pair11","Pair13")){
      hm<-matrix(nrow=2,ncol=length(tree$tip.label),dimnames = list(c("LOY","CNA"),tree$tip.label))
      hm[,"Ancestral"]<-"white"
      for(i in 1:ncol(hm)){hm[1,i]<-ifelse(!colnames(hm)[i]%in%Y_loss_df$id,"darkgrey",ifelse(Y_loss_df$loss_of_Y[Y_loss_df$id==colnames(hm)[i]]=="YES","red","lightgrey"))}
      for(i in 1:ncol(hm)){hm[2,i]<-ifelse(colnames(hm)[i]%in%CN_change_df$Sample,"purple","lightgrey")}
      tree=add_heatmap(tree,heatmap=hm,cex.label = 1)
    } else {
      hm<-matrix(nrow=1,ncol=length(tree$tip.label),dimnames = list(c("CNA"),tree$tip.label))
      hm[,"Ancestral"]<-"white"
      for(i in 1:ncol(hm)){hm[1,i]<-ifelse(colnames(hm)[i]%in%CN_change_df$Sample,"purple",ifelse(colnames(hm)[i]=="Ancestral","white","lightgrey"))}
      tree=add_heatmap(tree,heatmap=hm,cex.label = 1)
      }
    temp=add_annotation(tree,
                        annot_function=plot_sharing_info,
                        donor_ID=get_DR_ids(tree)['donor_ID'],
                        recip_ID=get_DR_ids(tree)['recip_ID'],
                        sharing_cols=c("black", "#11a0aa80", "#c8256580")
    )
    temp=plot_tree_labels(tree,
                          details = details,
                          type="line",
                          query.field = "shared_coding_change_chip", #alternative is 'coding_change_chip'
                          data.frame(value="yes",col="red",pch = 17,stringsAsFactors = FALSE), #if use 'coding_change_chip', value is 'Coding change mutation in driver'
                          label.field = "variant_ID",
                          cex.label = 0.8,
                          lty=2,
                          lwd=2)
    axis(side=4,at=mean(get_mut_burden(tree))-axis_at,labels = labels_at,las=2,cex.axis=0.7)
    transplant_time_median=muts_from_age(pair_age_at_transplant,tree,sampling_age=pair_age)
    #arrows(x0=-1,x1=1+length(all.trees.ultra[[i]]$tip.label),y0=mean(get_mut_burden(all.trees.ultra[[i]]))-transplant_time_median,y1=mean(get_mut_burden(all.trees.ultra[[i]]))-transplant_time_median,length=0,lty=2)
    CI_lower=transplant_time_median-2*sqrt(transplant_time_median)
    CI_upper=transplant_time_median+2*sqrt(transplant_time_median)
    rect(xleft = -1,
         xright=1+length(tree$tip.label),
         ybottom=mean(get_mut_burden(tree))-CI_upper,
         ytop=mean(get_mut_burden(tree))-CI_lower,
         col=rgb(0.1,0.1,0.1,alpha=0.3),border = NA)
  })
}
dev.off()

#Review contributions of individual signatures by sample at the bottom of the tree - here set to show APOBEC
par(mfrow=c(4,2))
temp=lapply(names(all.trees),function(pair){
  sig_in_samples=get_signatures_in_samples(tree=all.trees[[pair]],signature_names = c("N5"),exposures_df = exposures_df%>%dplyr::filter(Pair==pair))
  print(max(sig_in_samples))
  tree=plot_tree(all.trees[[pair]],title=pair,cex.label=0,bars=sig_in_samples)
  temp=add_annotation(tree,
                      annot_function=plot_sharing_info,
                      donor_ID=get_DR_ids(tree)['donor_ID'], recip_ID=get_DR_ids(tree)['recip_ID'],
                      sharing_cols=c("black", "#11a0aa80", "#c8256580")
  )
  text(x = 0, y=-0.05*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max contribution:",round(max(sig_in_samples),digits = 2)),pos = 4)
})

#----CALCULATE A RANGE OF PHYLOGENETIC DIVERSITY MEASURES TO COMPARE THE DONOR & RECIPIENT TREES
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

#First, view as a direct R~D comparison
cbind(Pair_metadata$Pair,sharedness_stats)%>%
  dplyr::rename(Pair="Pair_metadata$Pair")%>%
  ggplot(aes(x=D_stat,y=R_stat,col=Pair))+
  #geom_segment(aes(x=D_stat,xend=D_stat,y=D_stat,yend=R_stat),col="black",arrow=arrow(angle=15,length=unit(0.2,"cm")),linetype=2)+
  geom_point()+
  geom_abline(slope=1,linetype=1)+
  scale_x_continuous(limits=c(0,0.035))+
  scale_y_continuous(limits=c(0,0.035))+
  theme_classic()+
  labs(x="Sharedness of DONOR tree",y="Sharedness of RECIPIENT tree")+
  geom_polygon(data=data.frame(x=c(0,0,0.035),y=c(0,0.035,0.035)),aes(x=x,y=y),fill="grey",alpha=0.5,inherit.aes = F)+
  annotate(geom="text",x=0.009,y=0.033,label="Increased sharedness in recipient",fontface="italic")+
  annotate(geom="text",x=0.027,y=0.001,label="Decreased sharedness in recipient",fontface="italic")

#Then view as a log2(FC) from R~D
cbind(Pair_metadata$Pair,sharedness_stats)%>%
  dplyr::rename(Pair="Pair_metadata$Pair")%>%
  mutate(FC=R_stat/D_stat)%>%
  gather(key="D_or_R",value="Sharedness",-FC,-Pair,-comb_stat)%>%
  ggplot(aes(x=Pair,xend=Pair,y=0,yend=log2(FC),col=factor(ifelse(FC>1,1,0))))+
  geom_segment(arrow=arrow(angle=30,length=unit(1,"mm")))+
  scale_color_manual(values=c("#11a0aa", "#c82565"),guide="none")+
  scale_y_continuous(limits=c(-0.8,0.8))+
  geom_hline(yintercept=0,linetype=2)+
  labs(x="",y="log2(FC)",title="Change in tree 'sharedness': Recipient vs Donor")+
  theme_classic()+
  theme(axis.line.y = element_blank(),axis.ticks.y=element_blank())+
  coord_flip()

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

remove_names=function(x){attr(x,"names")<-NULL;return(x)}
new_pair_names=paste("Pair",1:nrow(Pair_metadata),sep = "_")
names(new_pair_names)<-Pair_metadata%>%arrange(Age)%>%pull(Pair)
Pair_metadata$Pair_new<-factor(remove_names(new_pair_names[Pair_metadata$Pair]),levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

Diversity_stats_df%>%
  #mutate(Pair=factor(Pair,levels=Pair_metadata%>%arrange(Age)%>%pull(Pair)))%>%
  gather(-Pair,-D_or_R,key="Stat",value="value")%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=new_pair_names))%>%
  mutate(Age=sapply(Pair_new,function(pair) Pair_metadata$Age[Pair_metadata$Pair_new==pair]))%>%
  filter(Stat=="MNTD")%>%
  mutate(Stat=ifelse(Stat=="MNTD","Clonal diversity (MNTD)"))%>%
  mutate(D_or_R=ifelse(D_or_R=="D","Donor","Recipient"))%>%
  ggplot(aes(x=Pair_new,y=value,col=D_or_R))+
  geom_line(aes(group=Pair_new),col="black",arrow = grid::arrow(angle = 30,length = unit(2,"mm")))+
  geom_point(alpha=0.5,size=3)+
  #facet_grid(cols=vars(Stat),scales="free")+
  theme_bw()+
  labs(x="",y="Clonal Diversity\n(Mean Nearest\nTaxon Distance)",title="",col="")

DR_change<-bind_rows(lapply(names(all.trees),function(pair) {
  pair_df<-Diversity_stats_df%>%dplyr::filter(Pair==pair)%>%dplyr::select(-Pair,-D_or_R)
  change=pair_df[2,]/pair_df[1,]
  return(change)
}))
rowMeans(DR_change)

##Do a "oligo-clonal" contribution diagram
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

expanded_clones_df<-dplyr::bind_rows(lapply(1:length(all.trees.ultra),function(i) {
tree<-all.trees.ultra[[i]]
  
  #Extract the donor/ recipient Ids and the donor/ recipient trees
  donor_ID=get_DR_ids(tree)['donor_ID']
  recip_ID=get_DR_ids(tree)['recip_ID']
  D_tree<-keep.tip(tree,grep(donor_ID,tree$tip.label,value=T))
  R_tree<-keep.tip(tree,grep(recip_ID,tree$tip.label,value=T))
  
  #Get the expanded clades from each tree
  D_nodes=get_expanded_clade_nodes(tree=D_tree,height_cut_off=50);D_df<-D_nodes%>%mutate(D_or_R="D")
  R_nodes=get_expanded_clade_nodes(tree=R_tree,height_cut_off=50);R_df<-R_nodes%>%mutate(D_or_R="R")
  
  return(rbind(D_df,R_df)%>%mutate(Pair=Pair_metadata$Pair[i]))
}))

expanded_clones_df%>%
  arrange(clonal_fraction)%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=D_or_R,y=clonal_fraction,fill=clonal_fraction))+
  geom_bar(stat="identity",position="stack",alpha=1,col="black")+
  geom_label(data=Pair_metadata,aes(x=1.5,y=1,label=paste("Age =",Age)),inherit.aes = F)+
  facet_grid(cols=vars(Pair_new),drop = F)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(4,"YlOrRd"))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
  theme_bw()+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction")
  

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

#Add the p.value from proportions test to see which are robustly different between donor & recipient
expanded_df_full$p.value=sapply(1:nrow(expanded_df_full),function(i) {
  prop.test(x=c(expanded_df_full$D_count[i],expanded_df_full$R_count[i]),n=c(expanded_df_full$D_total[i],expanded_df_full$R_total[i]))$p.value
})

#Pathway analysis on mutated genes that appear to be most 'selected for' in recipients vs donors
library(clusterProfiler)
library(tidyverse)
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

expanded_df_full%>%
  dplyr::filter(D_frac>0.02|R_frac>0.02)%>%
  #filter(p.value<0.1)%>%
  dplyr::select(-clonal_fraction)%>%
  mutate(cat=ifelse(p.value>0.05,"Similar fractions",ifelse(D_frac>R_frac,"Donor higher","Recip higher")))%>%
  gather(-nodes,-cat,-D_count,-R_count,-Pair,-D_total,-R_total,-p.value,key="D_or_R",value="clonal_fraction")%>%
  mutate(D_or_R=gsub("_frac","",D_or_R))%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=D_or_R,y=clonal_fraction,col=factor(nodes),group=factor(nodes)))+
  geom_point(alpha=0.6)+
  geom_line(arrow=arrow(length=unit(2,"mm")),linetype=2)+
  facet_grid(cat~Pair_new,scales="free")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction")
    

#Now look specifically at trajectories of clonal haem mutations between D -> R
CH_df<-Map(f=function(details,pair) {
  if(pair=="Pair11") {details<-dplyr::bind_rows(details,Pair11_loss_of_Y_details)}
  df<-details%>%
    dplyr::filter(coding_change_chip=="yes")%>%
    dplyr::select(node,Gene,variant_ID)%>%
    dplyr::rename(nodes=node)
  return(df)
},details=all.muts.nodups,pair=names(all.muts))

CH_df_full=dplyr::bind_rows(Map(f=function(pair,df,tree) {
  print(pair)
  if(nrow(df)==0) {stop(return(df))}
  #tree<-all.trees.ultra[[pair]]
  donor_ID=get_DR_ids(tree)['donor_ID']
  recip_ID=get_DR_ids(tree)['recip_ID']
  
  total_donor=sum(grepl(donor_ID,tree$tip.label))
  total_recip=sum(grepl(recip_ID,tree$tip.label))
  
  df$D_count=sapply(df$nodes,function(node) {sum(grepl(donor_ID,getTips(tree,node)))})
  df$R_count=sapply(df$nodes,function(node) {sum(grepl(recip_ID,getTips(tree,node)))})
  
  df$latest_acquisition<-sapply(df$nodes,function(node) {
    muts=nodeheight(tree,node)
    age=age_from_muts(muts,tree = tree,sampling_age = Pair_metadata%>%dplyr::filter(Pair==pair)%>%pull(Age))
    return(age)
    })
  
  return(df%>%mutate(Pair=pair,D_total=total_donor,R_total=total_recip)%>%mutate(D_frac=D_count/D_total,R_frac=R_count/R_total))
},pair=names(CH_df),df=CH_df,tree=all.trees.ultra))

CH_df_full%>%
  dplyr::filter(D_frac>0.02|R_frac>0.02)%>%
  gather(-nodes,-Gene,-latest_acquisition,-variant_ID,-D_count,-R_count,-Pair,-D_total,-R_total,key="D_or_R",value="clonal_fraction")%>%
  mutate(D_or_R=gsub("_frac","",D_or_R))%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=D_or_R,y=clonal_fraction,col=Pair_new,group=factor(nodes)))+
  geom_point(alpha=0.6)+
  geom_line(arrow=arrow(length=unit(2,"mm")),linetype=2)+
  facet_grid(~Gene,scales="free")+
  theme_bw()+
  #theme(legend.position = "none")+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction")

#Plot Acquisition timing of recipient-only drivers
library(ggrepel)
driver_mutation_timing_plot<-CH_df_full%>%
  dplyr::filter(D_frac<0.02 & R_frac>0.02)%>%
  left_join(Pair_metadata,by="Pair")%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=Pair,xend=Pair,y = Age_at_transplant,yend=Age))+
  #geom_segment(size=1,linetype=1,col="darkred")+
  geom_segment(arrow=arrow(angle = 10,ends = "last",length=unit(0.4,"cm")),size=1,linetype=1,col="darkred")+
  geom_segment(aes(x=Pair,xend=Pair,y=0,yend=Age_at_transplant),col="darkgrey",lineend = "round")+
  geom_point(aes(x=Pair,y=0),col="darkgrey")+
  geom_point(aes(x=Pair,y=Age_at_transplant),col="darkred")+
  geom_point(aes(x=Pair,y=latest_acquisition),col="red",alpha=0.6,size=2)+
  geom_label_repel(aes(x=Pair,y=latest_acquisition,label=variant_ID),nudge_y=5,nudge_x=-0.25)+
  scale_y_reverse(breaks = seq(0,80,10))+
  labs(x="",y="Age (years)",title="Latest acquisition of recipient-only drivers")+
  theme_classic()+
  theme(axis.line.x=element_blank(),axis.ticks.x = element_blank(),panel.grid.major.y = element_line())



#Plot the trajectory of CH driver clones, taking into account the presence of clones with multiple driver mutations
CH_df_full=dplyr::bind_rows(Map(f=function(pair,df,tree) {
  print(pair)
  if(nrow(df)==0) {stop(return(df))}

  donor_ID=get_DR_ids(tree)['donor_ID']; total_donor=sum(grepl(donor_ID,tree$tip.label))
  recip_ID=get_DR_ids(tree)['recip_ID']; total_recip=sum(grepl(recip_ID,tree$tip.label))
  
  sample_genotypes=sapply(tree$tip.label,function(SampleID) {
    included_nodes=get_ancestral_nodes(which(tree$tip.label==SampleID),edge = tree$edge)
    CH_variants<-df%>%dplyr::filter(nodes%in%included_nodes)%>%pull(variant_ID)%>%paste0(collapse=", ")
    return(CH_variants)
  })
  sample_genotypes<-sample_genotypes[sample_genotypes!=""]
  D_genotype_counts=sapply(unique(sample_genotypes),function(genotype) {sum(grepl(donor_ID,names(sample_genotypes))&sample_genotypes==genotype)})
  R_genotype_counts=sapply(unique(sample_genotypes),function(genotype) {sum(grepl(recip_ID,names(sample_genotypes))&sample_genotypes==genotype)})

  return(data.frame(Pair=pair,clone=unique(sample_genotypes),D_samps=D_genotype_counts,R_samps=R_genotype_counts,D_total=total_donor,R_total=total_recip)%>%mutate(D_frac=D_samps/D_total,R_frac=R_samps/R_total)%>%tibble::remove_rownames())
},pair=names(CH_df),df=CH_df,tree=all.trees.ultra))

#Make df that only includes the multiple
CH_summary_df<-CH_df_full%>%
  mutate(Gene=ifelse(grepl(",",clone),"Multiple drivers",str_remove(clone," .*")))%>%
  dplyr::filter(D_frac>0.02|R_frac>0.02)%>%
  gather(-clone,-Gene,-D_samps,-R_samps,-Pair,-D_total,-R_total,key="D_or_R",value="clonal_fraction")%>%
  mutate(D_or_R=gsub("_frac","",D_or_R))%>%
  mutate(Gene=factor(Gene,levels = c("DNMT3A","TET2","ASXL1","BCOR","LOY","Multiple drivers")))

multiple_CH_only<-CH_summary_df%>%
  dplyr::filter(Gene=="Multiple drivers" & D_or_R=="R")

library(ggrepel)
CH_summary_df%>%
  ggplot(aes(x=D_or_R,y=clonal_fraction,col=Pair,group=factor(clone)))+
  geom_point(alpha=0.6)+
  geom_line(arrow=arrow(length=unit(2,"mm")),linetype=2)+
  geom_label_repel(data=multiple_CH_only,aes(label=stringr::str_replace_all(clone,pattern=", ",replacement="\n")),size=2.5,nudge_y = 0.1,col="black")+
  facet_grid(~Gene,scales="free")+
  theme_bw()+
  #theme(legend.position = "none")+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction")


#Now do the ABC
#Define functions needed for script
get_ltt = function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  ltt_tree = sapply(time_points, function(x) {
    sum(nodeheights[,1] < x & !nodeheights[,2] < x)
  })
  return(ltt_tree)
}

get_coalescences = function(ltt) {
  coals=sapply(2:length(ltt), function(i) {return(ltt[i]-ltt[i-1])})
  return(coals)
}

#Go through each pair in turn
#1. Extract summary statistics from the data
#2. Import all the tree & parameter files from simulations & extract/ store the information in a matrix
#3. Save all the data in a list as a .RDS file

for(k in 1:length(all.trees.ultra)) {
  tree<-all.trees.ultra[[k]]
  recip_ID<-get_DR_ids(tree)['recip_ID']
  tree.R<-keep.tip(tree,grep(recip_ID,tree$tip.label,value=T))
  age_at_transplant<-Pair_metadata$Age_at_transplant[k]
  final_age<-Pair_metadata$Age[k]
  mean_transplant_time_mutations=muts_from_age(age = age_at_transplant,tree = tree,sampling_age = final_age) #The 0.558 is from the simulations
  final_mutations=get_mut_burden(tree)[1]
  
  #Calculate time points to assess coalescences - make this every five years, or 4 years for Pair 11 which has a shorter post-transplant period
  if(Pair_metadata$Pair[k]=="Pair11") {
    years_post_transplant=c(0,4,8)
  } else {
    years_post_transplant=seq(0,5*floor((final_age-age_at_transplant)/5),by=5)
  }
  
  plot.phylo(tree.R,direction="downwards",show.tip.label=F)
  transplant_time_mutations=mean_transplant_time_mutations
  phylo_units_per_year=(final_mutations-transplant_time_mutations)/(final_age-age_at_transplant)
  mutation_height_points = transplant_time_mutations + (years_post_transplant * phylo_units_per_year) #Convert the years into x coordinates
  abline(h=(final_mutations-mutation_height_points),col="red")
  
  #Set the time points to assess LTT/ coalescences
  df=dplyr::bind_rows(lapply(1:1000,function(j) {
    
    #Define the 'time of transplant' on the tree for each iteration
    #This is a random draw from a normal distribution with: mean = calculated number of mutations at the time of transplant, variance = mean (i.e. assume poisson variation)
    transplant_time_mutations=rnorm(1,mean=mean_transplant_time_mutations,sd = sqrt(mean_transplant_time_mutations))
    phylo_units_per_year=(final_mutations-transplant_time_mutations)/(final_age-age_at_transplant)
    mutation_height_points = transplant_time_mutations + (years_post_transplant * phylo_units_per_year) #Convert the years into x coordinates
    
    #Now use these mutation height 'cut offs' to extract the summary statistics from the ultra-metric tree
    ltt_data=get_ltt(tree.R,time_points = mutation_height_points); names(ltt_data)=paste("ltt",1:length(ltt_data),sep="_")
    coalescences_data=get_coalescences(ltt_data); names(coalescences_data)=paste("coals",1:length(coalescences_data),sep="_")
    fractional_ltt_increase=coalescences_data/ltt_data[1:(length(ltt_data)-1)]; names(fractional_ltt_increase)=paste("frac_increase",1:length(fractional_ltt_increase),sep="_")
    
    #Return these as a data frame with a single row
    combined_ss<-matrix(data=c(final_mutations,ltt_data,coalescences_data,fractional_ltt_increase),nrow = 1,dimnames = list(j,c("final_mutation_burden",names(ltt_data),names(coalescences_data),names(fractional_ltt_increase))))
    return(as.data.frame(combined_ss,stringsAsFactors=F))
  }))
  
  #Is there a way to incorporate different data summary statistic values? If not, use the means of each summary stat
  stats=colMeans(df)
  
  #Import simulation statistics
  #setwd("~/Mounts/Lustre/Zur_HSCT/")
  setwd("/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT")
  simulation_dir=paste0("ABC_models/",Pair_metadata$Pair[k],"_initial")
  tree_files=list.files(paste0(simulation_dir,"/trees/"),pattern="tree")
  parameter_files=list.files(paste0(simulation_dir,"/parameters/"),pattern="param")
  
  #Set up the parameters and stats matrices
  load(paste0(simulation_dir,"/parameters/",parameter_files[2500]))
  params.sim=matrix(NA,nrow=length(tree_files),ncol=length(unlist(this_sim_params)))
  colnames(params.sim)<-names(unlist(this_sim_params))
  
  stats.sim=matrix(NA,nrow=length(tree_files),ncol=length(stats))
  colnames(stats.sim)<-names(stats)
  
  #Now fill the matrices from the simulation files
  for(i in 1:length(tree_files)) {
    tree.ABC=read.tree(paste0(simulation_dir,"/trees/",tree_files[i]))
    load(paste0(simulation_dir,"/parameters/",parameter_files[i]))
    mean_muts_sim=mean(get_mut_burden(tree.ABC));names(mean_muts_sim)="mean_muts"
    tree.ABC=make.ultrametric.tree(tree = tree.ABC)
    transplant_time=this_sim_params$mean_mut_burden_at_transplant/this_sim_params$mean_mut_burden_at_sampling
    phylo_units=(1-transplant_time)/(this_sim_params$Time_from_transplant)
    time_points=transplant_time+(years_post_transplant*phylo_units)
    
    # plot(tree.ABC,show.tip.label=FALSE)
    # abline(v=time_points,col="red")
    ltt_sim=get_ltt(tree.ABC,time_points = time_points); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
    coalescences_sim=get_coalescences(ltt_sim); names(coalescences_sim)=paste("coals",1:length(coalescences_sim),sep="_")
    fractional_ltt_increase_sim=coalescences_sim/ltt_sim[1:(length(ltt_sim)-1)]; names(fractional_ltt_increase_sim)=paste("frac_increase",1:length(fractional_ltt_increase_sim),sep="_")
    stats = c(mean_muts_sim,ltt_sim,coalescences_sim,fractional_ltt_increase_sim)
    
    #Fill up the matrix
    params.sim[i,] <- unlist(this_sim_params)
    stats.sim[i,]<-stats
    if(i%%10==0) {print(i)}
  }
  
  #Log the appropriate parameters to get back onto a log scale
  params_to_log = c("LT_engrafted_cell_no", "HSC_pop_size")
  params.sim.log = params.sim
  params.sim.log[,colnames(params.sim.log) %in% params_to_log] = log10(params.sim.log[,colnames(params.sim.log) %in% params_to_log])
  log_param_names = colnames(params.sim)
  log_param_names[colnames(params.sim) %in% params_to_log] = paste0("log_", colnames(params.sim)[colnames(params.sim) %in% params_to_log])
  colnames(params.sim.log) = log_param_names
  
  saveRDS(list(tree.data=tree.R,
               stats.data=stats,
               params.sim=params.sim,
               params.sim.log=params.sim.log,
               stats.sim=stats.sim),file = paste0("ABC_models/",Pair_metadata$Pair[k],"_tables.RDS"))
}

#Now can do the actual ABC step with this data
library(abc)
library(gridExtra)
library(ggplot2)

setwd("~/Mounts/Lustre/Zur_HSCT")

all_data<-list()
for(k in 1:length(Pair_metadata$Pair)) {
  print(Pair_metadata$Pair[k])
  #Re-import
  data<-readRDS(paste0("ABC_models/",Pair_metadata$Pair[k],"_tables.RDS"))
  
  #DECIDE WHICH SUMMARY STATS TO USE
  #Define the sets of summary stats
  all_stats = 1:ncol(data$stats.sim)
  ltt_stats = grep("ltt",colnames(data$stats.sim))
  coals_stats = grep("coals",colnames(data$stats.sim))
  frac_increase_stats = grep("frac_increase",colnames(data$stats.sim))
  
  #Define the set to use
  sumstats_to_include = c(ltt_stats,coals_stats,frac_increase_stats)
  
  #SET THE TOLERANCE - need at least 250 'accepted' results to get a sense of the distribution
  #tol = 250/nrow(sumstat)
  tol=0.05
  print(tol)
  
  #RUN THE ABC USING DIFFERENT SETTINGS: rejection method, ridge regression, neural network regression
  results=list()
  results$rej<-abc(target = data$stats.data[sumstats_to_include], param = data$params.sim.log, sumstat = data$stats.sim[,sumstats_to_include], tol = tol, method = "rejection")
  results$ridge<-abc(target = data$stats.data[sumstats_to_include], param = data$params.sim.log, sumstat = data$stats.sim[,sumstats_to_include], tol = tol, hcorr = FALSE, method = "ridge",transf = "none")
  results$neural_net<-abc(target = data$stats.data[sumstats_to_include], param = data$params.sim.log, sumstat = data$stats.sim[,sumstats_to_include], tol = tol, hcorr = TRUE, method = "neuralnet",transf = "none")
  data$results<-results #Save these to the data frame
  
  #Plot the results
  CI_range=0.8; quantiles<-c(0.5*(1-CI_range),1-(0.5*(1-CI_range)))
  data$plots<-list()
  data$summary_dfs<-list()
  for(abc_type in c("rej","ridge","neural_net")){
    
    data_and_prior<-as.data.frame(data$results[[abc_type]][[ifelse(abc_type=="rej","unadj.values","adj.values")]])%>%mutate(type="posterior")%>%
      rbind(as.data.frame(data$params.sim.log)%>%mutate(type="prior"))
    
    #Calculate and print quantile results
    result_CI=rethinking::HPDI(data$results[[abc_type]][[ifelse(abc_type=="rej","unadj.values","adj.values")]][,"log_LT_engrafted_cell_no"],prob=0.95)
    hpdi=rethinking::HPDI(data$results[[abc_type]][[ifelse(abc_type=="rej","unadj.values","adj.values")]][,"log_LT_engrafted_cell_no"],prob=0.02)
    
    p1<-data_and_prior %>%
      ggplot(aes(x=log_LT_engrafted_cell_no, col=factor(type))) + 
      geom_density() +
      theme_classic() +
      theme(axis.title = element_text(size=10))+
      labs(x="Log10 long-term engrafting HSC number",y="Density")+
      annotate(geom="label",x=2.7,y=0.6,vjust=+0.2,label=paste("Maximum posterior density:",round(mean(10^hpdi))))+
      annotate(geom="label",x=2.7,y=0.3,vjust=+0.1,label=paste0(paste0(round(10^result_CI),collapse="-")," (",100*CI_range,"% confidence interval)"))
    
    p2<-data_and_prior %>%
      ggplot(aes(x=log_HSC_pop_size, col=factor(type))) + 
      geom_density() +
      theme_classic() +
      theme(axis.title = element_text(size=10))+
      labs(x="Log10 HSC population",y="Density")
    
    p3<-data_and_prior %>%
      filter(type=="posterior")%>%
      ggplot(aes(y=log_HSC_pop_size,x=log_LT_engrafted_cell_no)) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
      scale_y_continuous(limits = c(4.7,5.9)) +
      scale_x_continuous(limits=c(2,4)) +
      theme(axis.title = element_text(size=10))+
      theme_classic()+
      labs(y="Log10 HSC population",x="Log10 long-term engrafting HSC number")
    
    
    data$plots[[abc_type]]<-list(p1,p2,p3)
    data$summary_dfs[[abc_type]]<-data.frame(Pair=Pair_metadata$Pair[k],abc_type=abc_type,max_dens=round(mean(10^hpdi)),cilow=round(10^result_CI)[1],cihigh=round(10^result_CI)[2])
    
    all_data[[Pair_metadata$Pair[k]]]<-data
    #ggsave(plot=arrangeGrob(p1,p2,p3,ncol=3),filename = paste0("ABC_models/plots/ABC_",Pair_metadata$Pair[k],"_",abc_type,".pdf"))
  }
}


LT_engrafting_cell_no_df<-dplyr::bind_rows(lapply(all_data,function(list) {dplyr::bind_rows(list$summary_dfs)}))

LT_engrafting_cell_no_df%>%
  mutate(Pair=factor(Pair,levels=Pair_metadata%>%arrange(Age)%>%pull(Pair)))%>%
  ggplot(aes(x=Pair,y=max_dens,ymin=cilow,ymax=cihigh))+
  geom_point()+
  geom_errorbar()+
  facet_grid(rows=vars(abc_type))+
  coord_flip()

all_data$Pair11$results$neural_net$adj.values

library(ggridges)
abc_results<-dplyr::bind_rows(Map(function(list,Pair){return(as.data.frame(list$results$neural_net$adj.values,stringsAsFactors=F)%>%mutate(Pair=Pair))},list=all_data,Pair=names(all_data)))
HSC_posterior_ridge_plot<-abc_results%>%
  mutate(Pair=new_pair_names[Pair])%>%
  #mutate(Pair=factor(Pair,levels =Pair_metadata%>%arrange(Age)%>%pull(Pair) ))%>%
  ggplot(aes(x=10^log_LT_engrafted_cell_no,y=Pair))+
  geom_density_ridges(aes(fill=Pair)) +
  scale_x_log10()+
  scale_fill_brewer(palette="Set2")+
  #scale_fill_gradientn(olours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),name = "Log10")+
  theme_classic()+
  labs(x="Long-term engrafted HSCs",y="")
ggsave(filename="~/R_work/Zurich_HSCT/plots/HSC_posterior_ridge_plot.pdf",HSC_posterior_ridge_plot,width=6,height=3)  

