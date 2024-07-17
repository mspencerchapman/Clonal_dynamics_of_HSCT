#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest","pheatmap")

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
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

#========================================#
# Set the ggplot2 themes for plotting ####
#========================================#

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
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

## Read in Pair metadata data frame----
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

## Define colour themes for the Pairs & DorR----
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

## Read in other data objects ----
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))

## Extract objects from these lists in a 'for' loop ----
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

## Generate information regarding loss-of-Y in male samples from X and Y coverage data ----
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

## Read in the spreadsheet listing other copy number changes ----
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

## Read in mutational signature extraction data----
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)


#========================================#
# PLOT the COMBINED donor and recipient tree ####
#========================================#

#Plot the adjusted COMBINED ultrametric trees showing the transplant period

pdf(paste0(plots_dir,"ExtDatFig3.pdf"),width = 8,12)
par(mfrow=c(5,2))
temp=sapply(paste0("Pair_",1:10), function(Pair_new) {
  
  pair=Pair_metadata$Pair[Pair_metadata$Pair_new==Pair_new]
  tree<-all.trees.ultra[[pair]]
  details<-all.muts.nodups[[pair]]
  pair_age_at_transplant<-Pair_metadata%>%dplyr::filter(Pair==pair)%>%pull(Age_at_transplant)
  pair_age<-Pair_metadata%>%dplyr::filter(Pair==pair)%>%pull(Age)
  
  #Set up the axis for plotting
  binwidth=ifelse(pair_age>50,10,5)
  max_age_to_plot=binwidth*ceiling(pair_age/binwidth)
  axis_at=c(0,sapply(seq(0,max_age_to_plot,by=binwidth),muts_from_age,tree,pair_age))
  labels_at=c("Zyg.","Birth",seq(binwidth,max_age_to_plot,by=binwidth))
  
  #Define driver samples & expanded clade samples
  drivers_to_show=c("Oncogenic") #Can include possible drivers as well by making this c("Oncogenic","Possible")
  driver_samples=unlist(lapply(details%>%filter(Decision%in%drivers_to_show)%>%pull(node),function(node) getTips(tree,node)))
  expanded_clade_samples=unlist(lapply(get_expanded_clade_nodes(tree,min_clonal_fraction = 0.01)$nodes,function(node) getTips(tree,node)))
  D_expanded_clade_samples=unlist(lapply(get_DR_expanded_clades(tree,DorR = "D",min_samples = 4,min_clonal_fraction = 0.02),function(node) getTips(tree,node)))
  R_expanded_clade_samples=unlist(lapply(get_DR_expanded_clades(tree,DorR = "R",min_samples = 4,min_clonal_fraction = 0.02),function(node) getTips(tree,node)))
  
  #Plot the main tree
  tree=plot_tree(tree,cex.label = 0,plot_axis=F,vspace.reserve=0.1,title=Pair_new,lwd=0.25)
  
  #Generate the heatmaps that go beneath the tree
  #If the donor is male, include a row for 'loss of Y'
  if(pair%in%c("Pair11","Pair13","Pair25")){
    hm<-matrix(nrow=6,ncol=length(tree$tip.label),dimnames = list(c("Donor or Recip","Known Driver","Donor expanded","Recip expanded","LOY","CNA"),tree$tip.label))
    for(i in 1:ncol(hm)){hm[1,i]<-ifelse(grepl(get_DR_ids(tree)['donor_ID'],colnames(hm)[i]),"#11a0aa95","#c8256595")}
    for(i in 1:ncol(hm)){hm[2,i]<-ifelse(colnames(hm)[i]%in%driver_samples,"red","lightgrey")}
    for(i in 1:ncol(hm)){hm[3,i]<-ifelse(colnames(hm)[i]%in%D_expanded_clade_samples,"#66C2A5","lightgrey")}
    for(i in 1:ncol(hm)){hm[4,i]<-ifelse(colnames(hm)[i]%in%R_expanded_clade_samples,"#E7298A","lightgrey")}
    for(i in 1:ncol(hm)){hm[5,i]<-ifelse(!colnames(hm)[i]%in%Y_loss_df$id,"darkgrey",ifelse(Y_loss_df$loss_of_Y[Y_loss_df$id==colnames(hm)[i]]=="YES","yellow","lightgrey"))}
    for(i in 1:ncol(hm)){hm[6,i]<-ifelse(colnames(hm)[i]%in%CN_change_df$Sample,"purple","lightgrey")}
    hm[,"Ancestral"]<-"white"
    tree=add_heatmap(tree,heatmap=hm,cex.label = 0.5)
  } else {
    hm<-matrix(nrow=5,ncol=length(tree$tip.label),dimnames = list(c("Donor or Recip","Known Driver","Donor expanded","Recip expanded","CNA"),tree$tip.label))
    for(i in 1:ncol(hm)){hm[1,i]<-ifelse(grepl(get_DR_ids(tree)['donor_ID'],colnames(hm)[i]),"#11a0aa95","#c8256595")}
    for(i in 1:ncol(hm)){hm[2,i]<-ifelse(colnames(hm)[i]%in%driver_samples,"red","lightgrey")}
    for(i in 1:ncol(hm)){hm[3,i]<-ifelse(colnames(hm)[i]%in%D_expanded_clade_samples,"#66C2A5","lightgrey")}
    for(i in 1:ncol(hm)){hm[4,i]<-ifelse(colnames(hm)[i]%in%R_expanded_clade_samples,"#E7298A","lightgrey")}
    for(i in 1:ncol(hm)){hm[5,i]<-ifelse(colnames(hm)[i]%in%CN_change_df$Sample,"purple","lightgrey")}
    hm[,"Ancestral"]<-"white"
    tree=add_heatmap(tree,heatmap=hm,cex.label = 0.5)
  }
  temp=add_annotation(tree,
                      annot_function=plot_sharing_info,
                      donor_ID=get_DR_ids(tree)['donor_ID'],
                      recip_ID=get_DR_ids(tree)['recip_ID'],
                      sharing_cols=c("black", "#11a0aa80", "#c8256580")
  )
  temp=plot_tree_labels(tree,
                        details = all.muts.nodups[[pair]],
                        type="line",
                        query.field = "Decision", #alternative is 'coding_change_chip' or shared_coding_change_chip
                        data.frame(value=drivers_to_show,col=c("red","red"),pch = c(17,17),stringsAsFactors = FALSE), #if use 'coding_change_chip', value is 'Coding change mutation in driver'
                        label.field = "variant_ID",
                        cex.label = 1,
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
dev.off()
