##Script to look at all the mutations in known driver genes and analyze their features.
##Specifically
#(1) The numbers & types of mutations in different genes and their distribution across individuals
#(2) The timing of driver mutation acquisition
#(3) The relative sizes of driver mutation clones in donors vs recipients

#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest")
bioconductor_packages=c()

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

#----------------------------------
# Set the ggplot2 theme for plotting
#----------------------------------

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#----------------------------------
# Set the root directory and read in the necessary files
#----------------------------------

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
source(paste0(root_dir,"/data/HSCT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_list)) {assign(x,details_list[[x]])}

################# START ANALYSIS ################# 

#----------------------------------
# Create a driver dataframe with key driver info
#----------------------------------

driver_df<-Map(details=all.muts,pair=names(all.muts),function(details,pair) return(details%>%filter(Decision%in%c("Possible","Oncogenic"))%>%mutate(Pair=pair)))%>%bind_rows()
gene_order<-driver_df%>%group_by(Gene)%>%summarise(n=n())%>%arrange(desc(n))%>%pull(Gene)
driver_gene_numbers<-driver_df%>%
  mutate(Gene=factor(Gene,levels=gene_order))%>%
  mutate(mut_type=ifelse(grepl("splice",Type),"Splice variant",ifelse(grepl("stop_gained|frameshift",Type),"Truncating variant",ifelse(grepl("inframe",Type),"Inframe deletion","Missense variant"))))%>%
  dplyr::count(Gene,mut_type)%>%
  ggplot(aes(y=n,x=Gene,fill=mut_type))+
  geom_bar(stat="identity",col="black",linewidth=0.25)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),legend.key.size = unit(3,"mm"))+
  labs(x="Gene",y="# of unique mutations",fill="")
ggsave(filename=paste0(plots_dir,"driver_gene_numbers.pdf"),driver_gene_numbers,width=3,height=2)

driver_gene_heatmap<-driver_df%>%
  dplyr::count(Pair,Gene)%>%
  tidyr::complete(Pair,Gene,fill=list(n=0))%>%
  mutate(label=ifelse(n>0,n,""))%>%
  left_join(Pair_metadata,by="Pair")%>%
  bind_rows(driver_df%>%
              mutate(Gene=factor(Gene,levels=gene_order))%>%
              dplyr::count(Pair,Gene)%>%
              tidyr::complete(Pair,Gene,fill=list(n=0))%>%
              mutate(label=ifelse(n>0,n,""))%>%
              left_join(Pair_metadata,by="Pair")%>%group_by(Pair_new)%>%summarise(Gene="Total",n=sum(n),label=as.character(sum(n))))%>%
  mutate(Pair_new=factor(Pair_new,levels=Pair_metadata%>%arrange(Age)%>%pull(Pair_new)))%>%
  mutate(Gene=factor(Gene,levels=c(gene_order,"Total")))%>%
  ggplot(aes(x=Gene,y=Pair_new,fill=n,label=label))+
  geom_tile()+
  scale_fill_gradientn(colours=c("white",brewer.pal(8,name="YlOrRd")))+
  geom_text(size=1.5)+
  scale_x_discrete(position="top")+
  my_theme+
  theme(axis.text.x=element_text(angle=90),legend.position = "none")+
  labs(x="Gene",y="Pair")
ggsave(filename=paste0(plots_dir,"driver_gene_heatmap.pdf"),driver_gene_heatmap,width=3,height=2)

#----------------------------------
# Create a driver dataframe with additional info about molecular time of acquisition
#----------------------------------

driver_timing_df<-Map(tree=all.trees.ultra,details=all.muts.nodups,pair=names(all.trees.ultra),function(tree,details,pair) {
  drivers_df<-details%>%filter(Decision%in%c("Oncogenic","Possible"))
  tree=plot_tree(tree,cex.label = 0)
  temp=add_annotation(tree=tree,details=details,annot_function=highlight_nodes,nodes=drivers_df$node)
  
  timing_df=lapply(1:nrow(drivers_df),function(i) {
    node=drivers_df$node[i]
    heights=nodeHeights(tree = tree)[which(tree$edge[,2]==node),]
    return(data.frame(mut_ref=drivers_df$mut_ref[i],min_molecular_time=heights[1],max_molecular_time=heights[2]))
  })%>%bind_rows()%>%
    right_join(drivers_df)%>%
    mutate(PairID=pair,.before=1)
  return(timing_df)
})%>%bind_rows()

gene_cols=c("#aee39a", "#197959", "#4ad9e1", "#2f5672", "#24a5f7", "#6457d9", "#f3c5fa", "#8d2973", "#fd92fa", "#652cf6", "#e0079b", "#34f199", "#19a71f", "#bce333", "#6a9012", "#2cf52b", "#683d0d", "#fab899", "#ae301f", "#f47d0d", "#a27c59", "#f9bd3a", "#ff1c5d", "#f87574", "#ac82b4")
driver_mut_timing_plot<-driver_timing_df%>%
  mutate(Type=ifelse(Gene%in%gene_order[1:5],Gene,"Other"))%>%
  mutate(Type=factor(Type,levels=c(gene_order[1:5],"Other")))%>%
  mutate(Gene=factor(Gene,levels=gene_order))%>%
  arrange(Gene,max_molecular_time)%>%
  mutate(pos=1:nrow(.))%>%
  ggplot(aes(ymin=min_molecular_time,ymax=max_molecular_time,xmax=pos+0.4,xmin=pos-0.4,fill=factor(Gene)))+
  geom_rect()+
  theme_classic()+
  scale_fill_manual(values=gene_cols)+
  scale_y_reverse()+
  facet_grid(cols=vars(Type),scales="free",space="free")+
  my_theme+
  theme(strip.text.x = element_text(size=5),axis.line.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.key.size = unit(2.5,"mm"))+
  labs(y="Molecular time of acquisition",fill="Gene")
ggsave(filename = paste0(plots_dir,"driver_mut_timing_plot.pdf"),driver_mut_timing_plot,width=7,height = 2.5)


#----------------------------------
# CALCULATE PROPORTIONS OF CELLS HARBOURING DRIVERS
#----------------------------------

#Display the proportions of cells carrying at least one driver (with confidence intervals)
plot.proportion.with.drivers.CI<-Map(tree=all.trees.cc.nodups,details=all.muts.nodups,Pair=names(all.trees.cc.nodups),f=function(tree,details,Pair) {
  driver_nodes<-details%>%filter(Decision=="Oncogenic"|Decision=="Possible")%>%pull(node)
  driver_samples<-unique(unlist(lapply(driver_nodes,function(node) getTips(tree,node))))
  n_donor_drivers=sum(grepl(get_DR_ids(tree)['donor_ID'],driver_samples))
  n_recip_drivers=sum(grepl(get_DR_ids(tree)['recip_ID'],driver_samples))
  n_total_donor=sum(grepl(get_DR_ids(tree)['donor_ID'],tree$tip.label))
  n_total_recip=sum(grepl(get_DR_ids(tree)['recip_ID'],tree$tip.label))
  
  donor.binom=binom.test(n_donor_drivers,n_total_donor)
  recip.binom=binom.test(n_recip_drivers,n_total_recip)
  
  data.frame(Pair=Pair,
             DorR=c("D","R"),
             driver_prop=c(donor.binom$estimate,recip.binom$estimate),
             inv_driver_prop=1/c(donor.binom$estimate,recip.binom$estimate),
             driver_prop_CIlow=c(donor.binom$conf.int[1],recip.binom$conf.int[1]),
             driver_prop_CIhigh=c(donor.binom$conf.int[2],recip.binom$conf.int[2]))
})%>%
  dplyr::bind_rows()%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=DorR,y=driver_prop,ymin=driver_prop_CIlow,ymax=driver_prop_CIhigh,col=DorR))+
  geom_point(size=0.5)+
  geom_errorbar(width=0.4)+
  facet_grid(cols=vars(Pair_new))+
  scale_color_manual(values=DorR_cols)+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  my_theme+
  theme(legend.position = "none")+
  labs(x="Donor (D) or Recipient (R)",y="Proportion of samples with\n1 or more known driver mutations")
ggsave(filename = paste0(plots_dir,"plot.proportion.with.drivers.CI.pdf"),plot.proportion.with.drivers.CI,width=5,height = 2)

#Display the proportions of cells carrying different numbers of drivers
n_driver_props_df<-Map(tree=all.trees.cc.nodups,details=all.muts.nodups,Pair=names(all.trees.cc.nodups),f=function(tree,details,Pair) {
  donor_tips=grep(get_DR_ids(tree)['donor_ID'],tree$tip.label) #Define the node numbers of the donor samples
  recip_tips=grep(get_DR_ids(tree)['recip_ID'],tree$tip.label) #Define the node numbers of the recipient samples
  
  #Iterate through each sample and define how many drivers lie on branches leading up to that sample
  donor_driver_numbers=sapply(donor_tips,function(node) {all_sample_nodes=get_ancestral_nodes(node,tree$edge);return(details%>%filter(node%in%all_sample_nodes & (Decision=="Oncogenic"|Decision=="Possible"))%>%nrow(.))})
  recip_driver_numbers=sapply(recip_tips,function(node) {all_sample_nodes=get_ancestral_nodes(node,tree$edge);return(details%>%filter(node%in%all_sample_nodes & (Decision=="Oncogenic"|Decision=="Possible"))%>%nrow(.))})
  
  data.frame(Pair=Pair,
             DorR=rep(c("D","R"),each=4),
             n_drivers=rep(c(0,1,2,3),times=2),
             proportion=c((table(donor_driver_numbers)/length(donor_driver_numbers))[c("0","1","2","3")],(table(recip_driver_numbers)/length(recip_driver_numbers))[c("0","1","2","3")]))
})%>%
  dplyr::bind_rows()

plot.proportion.with.drivers<-n_driver_props_df%>%
  replace_na(replace = list(proportion=0))%>%
  mutate(Pair_new=factor(new_pair_names[Pair],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=DorR,y=proportion,fill=factor(n_drivers,levels=c(0,1,2,3))))+
  geom_bar(stat="identity",position="stack")+
  facet_grid(cols=vars(Pair_new))+
  scale_fill_brewer(palette = "YlOrRd")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  my_theme+
  labs(x="Donor (D) or Recipient (R)",y="Clonal fraction",fill="Number of\ndrivers")
ggsave(filename = paste0(plots_dir,"plot.proportion.with.drivers.pdf"),plot.proportion.with.drivers,width=5,height = 2)
